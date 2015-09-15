#include "gapped_fragment_graph.h"

#include <cmath>

#include <boost/format.hpp>
#include <boost/math/distributions.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GappedFragmentGraph"));

class Scorer {
public:
    Scorer(size_t K, size_t insert_size, double delta, size_t pair_kmer_num, size_t genome_len) : _K(K), _pair_kmer_num(pair_kmer_num), _genome_len(genome_len) {
        _gaussian.resize(2 * insert_size);

        boost::math::normal g(insert_size, delta);
        if (!_gaussian.empty()) {
            _gaussian[0] = boost::math::pdf(g, 0);
        }
        for (size_t i = 1; i < 2 * insert_size; ++i) {
            _gaussian[i] = _gaussian[i - 1] + boost::math::pdf(g, i);
        }
    }
    double score(size_t leni, size_t lenj, long d, size_t c) const {
        double sum = 0;

        for (int i = 0; i < leni - _K + 1; ++i) {	
            sum += gaussian(leni - _K + d + lenj - _K - i + 1) - gaussian(leni - _K + d - i);
        }
        double lam = _pair_kmer_num * sum / _genome_len;
        return log_poisson(lam, c);
    }
private:
    double log_poisson(double lam, int i) const {
        if (abs(lam - 0) <= 0.000001) {
            return -200000;
        }
        double sum = 0;
        for (int index = 1; index <= i; index ++) {
            sum += log(lam) - log(index);
        }
        return sum - lam;
    }
    double gaussian(int i) const {
        if (i < 0) return 0;
        else if (i >= _gaussian.size()) return 1;
        return _gaussian[i];
    }

    std::vector< double > _gaussian;
    size_t _K;
    size_t _pair_kmer_num;
    size_t _genome_len;
};


GappedFragmentGraph::GappedFragmentGraph(size_t k, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len): _k(k), _pair_kmer_cutoff(pair_kmer_cutoff), _pair_read_cutoff(pair_read_cutoff), _percent(percent), PAIR_KMER_NUM(0), GENOME_LEN(genome_len), INSERT_SIZE(0), DELTA(0.0) {
    BOOST_ASSERT(percent>=0.0 && percent<=1.0);
    _graph.resize(size);
}

GappedFragmentGraph::~GappedFragmentGraph() {
}

void GappedFragmentGraph::setPairKmerNumAndInsertSizeAndDelta(size_t pair_kmer_num, size_t insert_size, double delta) {
    PAIR_KMER_NUM = pair_kmer_num;
    INSERT_SIZE = insert_size;
    DELTA = delta;
}

void GappedFragmentGraph::addEdge(size_t from, size_t to, long distance, size_t kmer_num, size_t read_num) {
	
    for (EdgeList::iterator it = _graph[from].begin(); it != _graph[from].end(); ++it) {
        if (it->component_id == to) {
            it->distance += distance;
            it->kmer_cov += kmer_num;
            it->read_cov += read_num;
            return ;
        }
    }
    Edge e(to, distance, kmer_num, read_num, 0);
    _graph[from].push_back(e);
}

void GappedFragmentGraph::scoreAndRemoveNoise(const ComponentList& components) {

    Scorer scorer(_k, INSERT_SIZE, DELTA, PAIR_KMER_NUM, GENOME_LEN);
	std::vector< double > s;

    // scoring
    {
        size_t graph_size = 0, count = 0;
        for (NodeList::iterator i = _graph.begin(); i != _graph.end(); ++i, ++count) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if (j->kmer_cov < _pair_kmer_cutoff || j->read_cov < _pair_read_cutoff){
                    j = i->erase(j);
                    continue;
                }
                j->distance /= j->kmer_cov;
                size_t leni = components[count].length();
                size_t lenj =components[j->component_id].length();
                j->score = scorer.score(leni, lenj, j->distance - leni + _k, j->kmer_cov);
                s.push_back(j->score);
                ++j;
                ++graph_size;
            }
        }
    
        LOG4CXX_INFO(logger, boost::format("befor remove nosize graph_size=%d") % graph_size);
    }
    std::sort(s.begin(), s.end());
    double threshold = s.empty() ? 0 : s[(int)(s.size() * _percent)];
    // removeNoise
    {
        size_t graph_size = 0;
        LOG4CXX_INFO(logger, boost::format("removeNoise threshold=%f") % threshold);
        for(NodeList::iterator i=_graph.begin(); i!=_graph.end(); ++i) {
            EdgeList::iterator j=i->begin();
            while(j!=i->end()) {
                if (j->score < threshold) {
                    j = i->erase(j);
                    continue;
                }
                ++j;
                ++graph_size;
            }
        }
        LOG4CXX_INFO(logger, boost::format("after remove nosize graph_size=%d") % graph_size);
    }
}

void GappedFragmentGraph::outputLP(std::ostream& os) {
    
	size_t rn = 0;
	for (size_t i = 0; i < _graph.size(); ++i) {
		os << boost::format("var x_%d;") % i << std::endl;
	}

	for (size_t i = 0; i < _graph.size(); ++i) {
		for (EdgeList::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << boost::format("var e_%d_%d;") % i % it->component_id << std::endl; 
			os << boost::format("var E_%d_%d;") % i % it->component_id << std::endl;
			++rn;
		}
	}

	os << std::endl;
	os << "minimize z:  ";

	size_t count = 0;
	for (size_t i = 0; i < _graph.size(); ++i) {
        for(EdgeList::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
		    os << boost::format(" E_%d_%d + ") % i % it->component_id;
		    ++count;
		    if (count % 10 == 0)
			    os << std::endl;
        }
	}
	os << "0;" << std::endl << std::endl;

	size_t index = 1;
	for (size_t i = 0; i < _graph.size(); ++i) {
		for (EdgeList::const_iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << boost::format("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;") % index++ % it->component_id % i % i % it->component_id % it->distance << std::endl;
		}
	}

	for (size_t i = 0; i < _graph.size(); ++i) {
		for (EdgeList::iterator it = _graph[i].begin(); it != _graph[i].end(); ++it) {
			os << boost::format("s.t. con%d : E_%d_%d + e_%d_%d >= 0;") % index++ % i % it->component_id % i % it->component_id << std::endl;
			os << boost::format("s.t. con%d : E_%d_%d - e_%d_%d >= 0;") % index++ % i % it->component_id % i % it->component_id << std::endl;
		}
	}

	os << std::endl;
    os << "end;";
}

std::ostream& operator<<(std::ostream& os, const GappedFragmentGraph& g) {
    size_t index = 0;
    for (GappedFragmentGraph::NodeList::const_iterator i = g._graph.begin(); i != g._graph.end(); ++i, ++index) {
        for (GappedFragmentGraph::EdgeList::const_iterator j=i->begin(); j != i->end(); ++j) {
            if (j->distance > 0) {
                os << boost::format("%d\t%d\t%d\t%d\t%f") % index % j->component_id % j->distance % j->kmer_cov % j->score << std::endl;
            }
        }
    }
    return os;
}
