#ifndef SWGRAPH_HPP
#define SWGRAPH_HPP

#include <random>
#include <iostream>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************** class SWGraph *******************************

// This class constructs a small-world graph, based on the number of vertices,
// N, the number of connected neighbours "to the right", KR, and the "rewiring
// probability", prwr (see Watts and Strogatz 1998); for prwr = 0 the graph is
// a regular ring lattice

template<typename FltType>
class SWGraph {
public:
    using flt = FltType;
    using v_type = std::vector<flt>;
    using i_type = std::vector<int>;
    using graph_type = 
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    using adj_iter_type =
        typename boost::graph_traits<graph_type>::adjacency_iterator;
    using adj_iter_range = std::pair<adj_iter_type, adj_iter_type>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    SWGraph(int a_N, int a_KR, flt a_prwr, rand_eng& a_eng);
    int degree(int i) const
        { return (i > -1 && i < N) ? out_degree(i, gr) : 0; }
    void WriteGraphML(const std::string& outfilename) const;
    // public data members
    int N;
    int KR;
    bool complete;
    flt prwr;
    i_type k;
    graph_type gr;
    rand_eng& eng;
};

template<typename FltType>
SWGraph<FltType>::SWGraph(int a_N, int a_KR, flt a_prwr, rand_eng& a_eng) :
    N{a_N},
    KR{a_KR},
    complete{false},
    prwr{a_prwr},
    k(N, 0),
    gr(N),
    eng{a_eng}
{
    rand_uni uni(0, 1);
    // the graph has N vertices, add the required edges
    if (2*KR < N - 1) {
        // distribution of "rewiring distance" (we know that KR+1<=N-KR-1)
        rand_int uri(KR + 1, N - KR - 1);
        // construct edges to KR neighbours to the right; because the graph is
        // undirected, this implies that there are the same number of
        // neighbours to the left; note also that there will be no duplicate
        // edges
        for (int jd = 1; jd <= KR; ++jd) {
            for (int i = 0; i < N; ++i) {
                int j = (i + jd) % N;
                if (uni(eng) < prwr) {
                    // "rewire a presumed edge" between i and j
                    if (out_degree(i, gr) < N - 2) {
                        // vertex i is not connected to all other possible
                        // vertices (which does not include j)
                        int jdp = uri(eng); // possible new neighbour distance
                        int jp = (i + jdp) % N;
                        int cnt = 0; // "give up" after N attempts to rewire
                        while (cnt < N && edge(i, jp, gr).second) {
                            // if there already is an edge between i and jp
                            // (because of previous rewiring), try new possible
                            // neighbour distance
                            jdp = uri(eng);
                            jp = (i + jdp) % N;
                            ++cnt;
                        }
                        if (cnt == N) { // there were N unsuccessful attempts
                            // add edge between i and j
                            add_edge(i, j, gr);
                        } else {
                            // add "rewired edge" between i and jp
                            add_edge(i, jp, gr);
                        }
                    } else {
                        // add edge between i and j
                        add_edge(i, j, gr);
                    }
                } else {
                    // add edge between i and j
                    add_edge(i, j, gr);
                }
            }
        }
        // save degree for each vertex in k
        for (int i = 0; i < N; ++i) {
            k[i] = out_degree(i, gr);
        }
    } else {
        // assume complete graph (i.e., all vertices connected)
        complete = true;
        // save degree for each vertex in k
        for (int i = 0; i < N; ++i) {
            k[i] = N - 1;
        }
    }

}

template<typename FltType>
void SWGraph<FltType>::WriteGraphML(const std::string& outfilename) const
{
    std::ofstream outfile(outfilename.c_str(), std::ios_base::out);
    if (!outfile) {
        std::cout << "Cannot open " << outfilename << ", cannot save data \n";
    } else {
        // create "empty" dynamic properties, to satisfy function call
        boost::dynamic_properties dp;
        // write graph to file
        write_graphml(outfile, gr, dp);
        outfile.close();
    }
}

#endif // SWGRAPH_HPP
