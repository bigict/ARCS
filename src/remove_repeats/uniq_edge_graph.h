#ifndef uniq_edge_graph_h_
#define uniq_edge_graph_h_

#include <map>
#include <string>
#include <vector>
#include <set>

#include "component.h"

class UniqEdgeGraph {
public:
    UniqEdgeGraph(size_t K, size_t edge_num, size_t max_overlap, size_t iteration)
        : _K(K), _edge_num(edge_num), _max_overlap(max_overlap), _iteration(iteration) {
    }
    virtual ~UniqEdgeGraph() {
    }

    void addEdge(size_t i, size_t j, int distance, int count, int score);
    void removeEdge(size_t i, size_t j);
    void removeNode(size_t i);
    bool hasEdge(size_t i, size_t j) const;

	bool input_edge_position();
	bool input_edge_length();
	bool input_edge_link(const std::string& edge_link_file);
	bool input_component();

	void linearize();
	void output_graph(const std::string& file) const;
	
private:
    friend class ConflictResolver;
    struct Edge {
        Edge(int distance = 0, int count = 0, int score = 0) : distance(distance), count(count), score(score) {
        }
        int distance;
        int count;
        int score;
    };
    typedef std::map< size_t, Edge > Children;
    typedef std::set< size_t > Parents;
    struct Node {
        Children children;
        Parents parents;
        operator bool() const {
            return !(children.empty() && parents.empty());
        }
    };
    typedef std::map< size_t, Node > NodeList;

    int getDistance(size_t i, size_t j) const;
	int getAncestor(size_t i, size_t j) const;
	int getDescendant(size_t i, size_t j) const;
    struct EdgeInfo {
        int from;
        int to;
        Edge edge;
        bool operator < (const EdgeInfo& o) const {
            return edge.score < o.edge.score;
        }
    };
    int EdgeScorePlus(int l, const EdgeInfo& info) const;
	EdgeInfo getMinEdge(size_t i) const;

    NodeList _nodelist;
    typedef std::vector< size_t > PositionList;
    PositionList _position_tbl;
    typedef std::vector< size_t > LengthList;
    LengthList _length_tbl;
    typedef std::vector< Component > ComponentList;
    ComponentList _component_tbl;

    size_t _K;
    size_t _edge_num;
    size_t _max_overlap;
    size_t _iteration;
};

#endif // uniq_edge_graph_h_
