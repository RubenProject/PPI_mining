#include <map>
#include <string.h>
#include <igraph.h>


#define EDGELIST_LOC "edgelist.txt"
#define WEIGHTLIST_LOC "weightlist.txt"
#define COMPLETE_LOC "edge_weightlist.txt"
#define MAPPING_LOC "mapping.txt"

#define NCOL_LOC "links.txt"
#define GENES_LOC "genes.txt"

/* globals */
igraph_t graph;
igraph_vector_t score;
igraph_integer_t vcount, ecount, acount;

igraph_strvector_t v_names;
igraph_strvector_t a_names;
igraph_vector_t a_idx;

//prototypes
void display_score(FILE *file);

//XXX: DEPRECATED
void create_graph_2(){
    FILE *fp_el;
    igraph_t g;
    fp_el = fopen("edgelist.txt", "r");
    if (fp_el == NULL){
        printf("Read error\n");
        exit(0);
    }
    igraph_read_graph_edgelist(&g, fp_el, 0, 1);
    fclose(fp_el);
}


//XXX: DEPRECATED
void create_edgelist(const char *file_loc){
    std::map<std::string, int> vdict;
    FILE *fp_in, *fp_out_e, *fp_out_w, *fp_out_ew, *fp_out_m;
    char *v0, *v1;
    int w, vid;
    fp_in = fopen(file_loc, "r");
    fp_out_e = fopen(EDGELIST_LOC, "w");
    fp_out_w = fopen(WEIGHTLIST_LOC, "w");
    fp_out_ew = fopen(COMPLETE_LOC, "w");
    fp_out_m = fopen(MAPPING_LOC, "w");
    if (fp_out_e == NULL || fp_out_e == NULL 
        || fp_out_w == NULL || fp_out_ew == NULL){
        printf("Read/write error\n");
        exit(0);
    }
    v0 = (char *)malloc(sizeof(char) * 50);
    v1 = (char *)malloc(sizeof(char) * 50);
    /* discard first row --lazily*/
    fscanf(fp_in, "%s%s%s", v0, v0, v0);


    vid = 1;
    while (fscanf(fp_in, "%s%s%i\n", v0, v1, &w) != EOF){
        //printf("%s %s %i\n", v0, v1, w);
        if (!vdict[v0]){
            vdict[v0] = vid;
            vid++;
        }
        if (!vdict[v1]){
            vdict[v1] = vid;
            vid++;
        }
        fprintf(fp_out_e, "%i %i\n", vdict[v0], vdict[v1]);
        fprintf(fp_out_w, "%i\n", w);
        fprintf(fp_out_ew, "%i %i %i\n", vdict[v0], vdict[v1], w);


    }
    std::map<std::string, int>::iterator it;
    for (it = vdict.begin(); it != vdict.end(); it++){
        fprintf(fp_out_m, "%s %i\n", it->first.c_str(), it->second);
    }

    fclose(fp_in);
    fclose(fp_out_e);
    fclose(fp_out_w);
    fclose(fp_out_ew);
    fclose(fp_out_m);
    free(v0);
    free(v1);
}


igraph_t create_graph(){
    FILE *fp_el;
    igraph_t g;
    fp_el = fopen(NCOL_LOC, "r");
    if (fp_el == NULL){
        printf("Read error\n Exiting...\n");
        exit(0);
    }
    /* important */
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_ncol(&g, fp_el, NULL, 1, IGRAPH_ADD_WEIGHTS_YES, IGRAPH_DIRECTED);
    return g;
}


void basic_stats(){
    int i;
    igraph_bool_t connected;
    igraph_integer_t diameter;
    igraph_vector_t degree;
    igraph_real_t centralization;
    igraph_real_t th_max;
    igraph_vector_init(&degree, vcount);

    /* check weakly connected */
    igraph_is_connected(&graph, &connected, IGRAPH_WEAK);
    printf("vertex count: %i\n", (int)vcount);
    printf("edge count: %i\n", (int)ecount);
    if (connected)
        printf("weakly connected\n");
    igraph_diameter(&graph, &diameter, 0, 0, 0, IGRAPH_UNDIRECTED, connected);
    printf("diameter: %li\n", (long int)diameter);
    igraph_centralization_degree(&graph, &degree, IGRAPH_ALL, 0, &centralization, &th_max, 0);
    printf("centralization: %f\n", (double)centralization);
    printf("theoretical max: %f\n", (double)th_max);
    for (i = 0; i < vcount; i++){
        if (VECTOR(degree)[i] > 0){
            printf("v: %i degree: %f\n", i, VECTOR(degree)[i]);
        }
    }
}


/* get the list of genes associated with menier's disease */
igraph_strvector_t get_associate_names(){
    /* There is 106 associated genes */
    const int n = 106;
    int i;
    FILE *fp;
    igraph_strvector_init(&a_names, n);
    fp = fopen(GENES_LOC, "r");
    if (fp == NULL){
        printf("Read error\n Exiting...\n");
        exit(0);
    }
    for (i = 0; i < n; i++){
        fscanf(fp, "%s\n", STR(a_names, i));
    }
    fclose(fp);
    return a_names;
}


/* cuts off skip amount of characters of the string and then converts */
int parse_name(const char *s, int skip){
    return std::stoi(std::string(&s[skip]));
}



/*find idx by vname, return -1 if not found */
int get_idx_by_vname(const char *name){
    int i, s_idx, c_idx;
    s_idx = parse_name(name, 4);
    /* linear search since list not sorted */
    for (i = 0; i < vcount; i++) {
        c_idx = parse_name(STR(v_names, i), 9);
        if (s_idx == c_idx){
            return i;
        }
    }
    return -1;
}


/* build a reverse index for the associate proteins */
igraph_vector_t get_associate_idx(){
    int i;

    igraph_vector_init(&a_idx, acount);
    for (i = 0; i < acount; i++){
        VECTOR(a_idx)[i] = get_idx_by_vname((char *)STR(a_names, i));
    }
    return a_idx;
}


/*main function that calculates score*/
void mine_graph(){
    int i, j, k, idx, nwalks, nsteps;
    igraph_vector_t neis, p_neis, history;
    igraph_integer_t vid, nid, eid, weight, hist_size; 
    double tot_weight, p, init_score, v_score, share;

    igraph_vector_init(&score, vcount);
    a_names = get_associate_names();
    a_idx = get_associate_idx();

    /* pagerank */
    hist_size = 100;
    igraph_vector_init(&neis, 0);
    igraph_vector_init(&history, hist_size);
    srand(42);
    nwalks = 10000;
    nsteps = 1000;
    init_score = 100.0;
    share = 0.0001;
    /* init associated proteins */
    for (i = 0; i < acount; i++){
        idx =(int)VECTOR(a_idx)[i]; 
        if (idx >= 0)
            VECTOR(score)[idx] = init_score;
    }
    for (i = 0; i < nwalks; i++){
        /* init random walk from random vertex*/
        vid = rand() % vcount;
        igraph_vector_fill(&history, -1);
        printf("iteration: %i\n starting at v: %i\n", i, vid);
        for (j = 0; j < nsteps; j++){
            /* select next vertex and score nieghbours*/
            igraph_neighbors(&graph, &neis, vid, IGRAPH_OUT);
            /* no neighbours */
            if (igraph_vector_size(&neis) == 0)
                break;
            v_score = VECTOR(score)[vid];
            tot_weight = 0;
            /* proportional selection */
            for (k = 0; k < igraph_vector_size(&neis); k++){
                nid = VECTOR(neis)[k]; 
                igraph_get_eid(&graph, &eid, vid, nid, 1, 0);
                weight = EAN(&graph, "weight", eid);
                VECTOR(score)[nid] += share * weight * v_score;
                if (!igraph_vector_contains(&history, nid)){
                    tot_weight += weight;
                } else {
                    VECTOR(neis)[k] = -1; 
                }
            }

            igraph_vector_init(&p_neis, igraph_vector_size(&neis));
            for (k = 0; k < igraph_vector_size(&neis); k++){
                nid = VECTOR(neis)[k]; 
                if (nid == -1){
                    if (k > 0)
                        VECTOR(p_neis)[k] = VECTOR(p_neis)[k-1];
                    else
                        VECTOR(p_neis)[k] = 0;
                } else {
                    igraph_get_eid(&graph, &eid, vid, nid, 1, 0);
                    weight = EAN(&graph, "weight", eid);
                    VECTOR(p_neis)[k] = weight / tot_weight;
                    if (k > 0){
                        VECTOR(p_neis)[k] += VECTOR(p_neis)[k-1];
                    }
                }
            }
            /* spin the wheel */
            p = (double)rand() / RAND_MAX;
            k = 0;
            while (p > VECTOR(p_neis)[k]){
                k++;
            }
            nid = VECTOR(neis)[k];
            VECTOR(history)[j % hist_size] = vid;
            vid = nid;
        }
    }
}

/* initialize basic properties */
void init(){
    int i;

    graph = create_graph();
    vcount = igraph_vcount(&graph);
    ecount = igraph_ecount(&graph);
    acount = 106; /* counted from file */
    igraph_strvector_init(&v_names, vcount);
    for (i = 0; i < vcount; i++){
        igraph_strvector_set(&v_names, i, VAS(&graph, "name", i));
    }
}

int sort_score(igraph_strvector_t names, igraph_vector_t score){
    int i, j, k, n;
    double t0;
    char *t1;
    /* 100 chars should be plenty */
    t1 = (char *)malloc(100 * sizeof(char));
    if (igraph_vector_size(&score) != igraph_strvector_size(&names)){
        printf("size mismatch\n");
        return 0;
    }
    n = igraph_vector_size(&score);
    for (i = 0; i < n - 1; i++) {
        k = i;

        /* Find element with smallest time */
        for (j = i + 1; j < n; j++) {
            if (VECTOR(score)[j] > VECTOR(score)[k]) {
                k = j;
            }
        }

        /* Swap k-th and i-th element */
        if (k != i) {
            t0 = VECTOR(score)[k];
            strcpy(t1, STR(names, k));

            VECTOR(score)[k] = VECTOR(score)[i];
            VECTOR(score)[i] = t0;

            igraph_strvector_set(&names, k, STR(names, i));
            igraph_strvector_set(&names, i, t1);
        }
    }
    return 1;
}

/*show top 100*/
void display_score(FILE *file){
    int i;
    if (!sort_score(v_names, score)){
        return;
    }
    for (i = 0; i < 100; i++){
        fprintf(file, "%s : %i\n", STR(a_names, i), (int)VECTOR(score)[i]);
    }
}


int main(void){
    init();
    //basic_stats();
    mine_graph();
    display_score(stdout);
    /* TODO not all objects are destroyed */
    igraph_destroy(&graph);
    igraph_strvector_destroy(&a_names);
    igraph_vector_destroy(&a_idx);
    return 0;
}
