#ifndef linearization_h_
#define linearization_h_

#include <string>
#include <vector>

#include "uniq_edge_graph.h"

class Linearization {
public:
    Linearization(size_t K) : _K(K) {
    }
    void initialize_kmer_hash();
    void initialize_kmer_hash_for_trainning();
    void free_kmer_hash();
    void initialize_edge(const std::string&);
    void initialize_scaf();
    void map_read(string, string);
    void linearize(string, string);
    void print_edge();
    void output_para();
    void estimate_insert_size();
private:
    struct Edge {
        std::string seq;
        int copy_num;
    };
    struct Component {
        size_t length;
        size_t gap;
    };
    typedef std::vector< Component > ComponentList;

    Uniq_Edge_Graph u_e_g;
    
    std::vector<Edge> edge;    

    std::vector< vector<unsigned int> > component;
    std::vector< vector<unsigned int> > gap;

    size_t _K;
};

#endif // linearization_h_
