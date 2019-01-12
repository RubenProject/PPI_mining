CXX = g++
LDFLAGS = $(shell pkg-config --libs igraph)
CFLAGS = $(shell pkg-config --cflags igraph) -Wall


all: ppi

ppi: ppi.cpp
	$(CXX) ppi.cpp $(CFLAGS) $(LDFLAGS) -o ppi
