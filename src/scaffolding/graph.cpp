#include "graph.h"

#include <cmath>
#include <boost/format.hpp>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Contiging.main"));

Graph::Graph(size_t k, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len): _k(k), _pair_kmer_cutoff(pair_kmer_cutoff), _pair_read_cutoff(pair_read_cutoff), _percent(percent), PAIR_KMER_NUM(0), GENOME_LEN(genome_len), INSERT_SIZE(0), DELTA(0.0) {
    BOOST_ASSERT(percent>=0.0 && percent<=1.0);
    _graph.resize(size);
}

Graph::~Graph() {
}

void Graph::setPairKmerNumAndInsertSizeAndDelta(size_t pair_kmer_num, size_t insert_size, double delta) {
    PAIR_KMER_NUM = pair_kmer_num;
    INSERT_SIZE = insert_size;
    DELTA = delta;
}

void Graph::addEdge(size_t from, size_t to, long dis, bool isread) {
	
    for(GraphEdge::iterator it=_graph[from].begin(); it!=_graph[from].end(); ++it) {
        if (it->component_id == to) {
            it->dis += dis;
            ++it->kmer_cov;
            if (isread) {
                ++it->read_cov;
            }
            return ;
        }
    }
    Edge e(to, dis, 1, (isread?1:0), 0);
    _graph[from].push_back(e);
}

void Graph::scoreAndRemoveNoise(const std::vector<Component>& com) {

    size_t graph_size = 0;
	initialize_gaussian();
	std::vector<double> s;
    size_t count = 0;
    for(GraphNode::iterator i=_graph.begin(); i!=_graph.end(); ++i, ++count) {
        GraphEdge::iterator j=i->begin();
        while(j!=i->end()) {
            if (j->kmer_cov < _pair_kmer_cutoff || j->read_cov < _pair_read_cutoff){
                j = i->erase(j);
                continue;
            }
            j->dis /= j->kmer_cov;
            size_t leni = com[count].getLen();
            size_t lenj =com[j->component_id].getLen();
            j->score = score(leni, lenj, j->dis - leni + _k, j->kmer_cov);
            s.push_back(j->score);
            ++j;
            ++graph_size;
        }
    }
    
    LOG4CXX_INFO(logger, boost::format("befor remove nosize graph_size=%d") % graph_size);
    
    sort(s.begin(), s.end());
    double threshold = s[(int)(s.size() * _percent)];
    
    graph_size = 0;
    LOG4CXX_INFO(logger, boost::format("removeNoise threshold=%f") % threshold);
    for(GraphNode::iterator i=_graph.begin(); i!=_graph.end(); ++i) {
        GraphEdge::iterator j=i->begin();
        while(j!=i->end()) {
            if (j->score < threshold) {
                j = i->erase(j);
                continue;
            }
            ++graph_size;
            ++j;
        }
    }
    LOG4CXX_INFO(logger, boost::format("after remove nosize graph_size=%d") % graph_size);

    LOG4CXX_INFO(logger, boost::format("INSERT_SIZE=%d") % INSERT_SIZE);
    
    LOG4CXX_INFO(logger, boost::format("DELTA=%f") % DELTA);

    LOG4CXX_INFO(logger, boost::format("PAIR_KMER_NUM=%d") % PAIR_KMER_NUM);

    LOG4CXX_INFO(logger, boost::format("GENOME_LEN=%d") % GENOME_LEN);

    return ;
}

void Graph::outputLP(std::ostream& os) {
    
	size_t rn = 0;
	for (size_t i = 0; i < _graph.size(); ++i) {
		os << "var x_" << i << ";" << std::endl;
	}

	for(size_t i = 0; i < _graph.size(); ++i) {
		for(GraphEdge::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << "var e_" << i << "_" << it->component_id  << ";" << std::endl; 
			os << "var E_" << i << "_" << it->component_id  << ";" << std::endl;
			++rn;
		}
	}

	os << std::endl;
	os << "minimize z:  ";

	size_t count = 0;
	for(size_t i = 0; i < _graph.size(); ++i) {
        for(GraphEdge::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
		    os << " E_" << i << "_" << it->component_id << " + ";
		    ++count;
		    if (count % 10 == 0)
			    os << std::endl;
        }
	}
	os << "0;" << std::endl << std::endl;

	size_t index = 1;
	for(size_t i = 0; i < _graph.size(); ++i) {
		for(GraphEdge::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << "s.t. con" << index ++ << " : x_" << it->component_id << " - x_" << i << " + e_" << i<< "_" << it->component_id << " = " << it->dis << ";" << std::endl;
		}
	}

	for(size_t i = 0; i < _graph.size(); ++i) {
		for(GraphEdge::iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << "s.t. con" << index ++ << " : E_" << i << "_" << it->component_id << " + e_" << i << "_" << it->component_id << " >= 0;" << std::endl;
			os << "s.t. con" << index ++ << " : E_" << i << "_" << it->component_id << " - e_" << i << "_" << it->component_id << " >= 0;" << std::endl;
		}
	}

	os << std::endl << "end;";

    return;
}

std::ostream& operator<<(std::ostream& os, const Graph& g) {
    size_t index = 0;
    for(Graph::GraphNode::const_iterator i=g._graph.begin(); i!=g._graph.end(); ++i, ++index) {
        for(Graph::GraphEdge::const_iterator j=i->begin(); j!=i->end(); ++j){
            if (j->dis > 0) {
                os << boost::format("%d\t%d\t%d\t%d\t%f") % index % j->component_id % j->dis % j->kmer_cov % j->score << std::endl;
                //os << index << "|" << j->component_id << "|" << j->dis << "|" << j->kmer_cov << "|" << j->score << std::endl; 
            }
        }
    }
    return os;
}

double Graph::score(size_t leni, size_t lenj, long d, size_t c) {
	double sum = 0;

	for (int i = 0; i < leni-_k+1; ++i) {	
		sum += get_Gaussian(leni - _k + d + lenj - _k - i + 1) - get_Gaussian(leni - _k + d - i);
	}
	double lam = PAIR_KMER_NUM * sum / GENOME_LEN;
	return log_poisson(lam, c);
}

double Graph::get_Gaussian(int i) {
    if (i < 0)
		return 0;
	else if ( i >= 2*INSERT_SIZE)
		return 1;
	else
		return Gaussian[i];
}

double Graph::log_poisson(double lam, int i) {
	if( abs(lam - 0) <= 0.000001) {
		return -200000;
	} else {	
		double sum = 0;
		for (int index = 1; index <= i; index ++) {
			sum += log(lam) - log( index);
		}
		return sum - lam;
	}
}

void Graph::initialize_gaussian() {
    std::vector<double> gaussian;

    int INSERT_SIZE2 = INSERT_SIZE;
	gaussian.resize(2*INSERT_SIZE2);
	Gaussian.resize(2*INSERT_SIZE2);
	
	double temp0, temp1, temp2;

	for(int i = 0; i < 2*INSERT_SIZE2; i++) {
		temp0 = -pow(i - INSERT_SIZE2, 2)/(2.0*DELTA*DELTA);
		temp1 = exp(temp0);
		temp2 = sqrt(2*3.141592653)*DELTA;
		gaussian[i] = temp1 / temp2;
	}

	double sum = 0;
	
	for(int i = 0; i < 2*INSERT_SIZE2; i++){
		sum += gaussian[i];
		Gaussian[i] = sum;
	}	
}

