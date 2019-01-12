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
    igraph_read_graph_ncol(&g, fp_el, NULL, 1, IGRAPH_ADD_WEIGHTS_YES, IGRAPH_UNDIRECTED);
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
    int i;

    a_names = get_associate_names();
    a_idx = get_associate_idx();

    for (i = 0; i < acount; i++){
        //do some action for all associate genes
        //printf("%s\n", STR(a_names, i));
    }
}

/* initialize basic properties */
void init(){
    int i;

    graph = create_graph();
    vcount = igraph_vcount(&graph);
    ecount = igraph_ecount(&graph);
    acount = 106; /* counted from file */
    igraph_vector_init(&score, vcount);
    igraph_strvector_init(&v_names, vcount);
    for (i = 0; i < vcount; i++){
        igraph_strvector_set(&v_names, i, VAS(&graph, "name", i));
    }
    /*XXX:weights could be similarly queried for easy manual use */
}

/*show top 100*/
void display_score(){
    int i;
    /* TODO: disiplay only highest scores: sort and save index*/
    for (i = 0; i < 100; i++){
        printf("%i\n", (int)VECTOR(score)[i]);
    }
}


int main(void){
    init();
    //basic_stats();
    mine_graph();
    //display_score();
    /* TODO not all objects are destroyed */
    igraph_destroy(&graph);
    igraph_strvector_destroy(&a_names);
    igraph_vector_destroy(&a_idx);
    return 0;
}
