#include "gapped_fragment_graph.h"
#include "gmm.h"

#include <cmath>
#include <limits>

#include <fstream>
#include <map>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/math/distributions.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GappedFragmentGraph"));

class Scorer {
public:
    Scorer(size_t K, size_t insert_size, double delta, size_t pair_kmer_num, size_t genome_len) : _K(K), _pair_kmer_num(pair_kmer_num), _genome_len(genome_len), _gaussian(insert_size, delta) {
    }
    double score(size_t leni, size_t lenj, long gap, size_t c) const {
        BOOST_ASSERT(leni >= _K - 1);
        BOOST_ASSERT(lenj >= _K - 1);

        double sum = 0;
        for (size_t i = 0; i < leni - _K + 1; ++i) {	
            sum += gaussian(leni + gap + lenj - i - _K + 1) - gaussian(leni + gap - i);
        }
        double lambda = _pair_kmer_num * sum / _genome_len;
        return log_poisson(lambda, c);
    }
private:
    double log_poisson(double lambda, size_t c) const {
        if (lambda <= std::numeric_limits< double >::epsilon()) {
            return -std::numeric_limits< double >::max();
        }
        double sum = c * log(lambda);
        for (size_t index = 1; index <= c; ++index) {
            sum -= log(index);
        }
        return sum - lambda;
    }
    double gaussian(int i) const {
        return boost::math::cdf(_gaussian, i);
    }

    boost::math::normal _gaussian;
    size_t _K;
    size_t _pair_kmer_num;
    size_t _genome_len;
};


GappedFragmentGraph::GappedFragmentGraph(size_t K, size_t pair_kmer_cutoff, size_t pair_read_cutoff, double percent, size_t size, size_t genome_len): _K(K), _pair_kmer_cutoff(pair_kmer_cutoff), _pair_read_cutoff(pair_read_cutoff), _percent(percent), PAIR_KMER_NUM(0), GENOME_LEN(genome_len), INSERT_SIZE(0), DELTA(0.0){
    BOOST_ASSERT(percent>=0.0 && percent<=1.0);
    _nodelist.resize(size);
}

GappedFragmentGraph::~GappedFragmentGraph() {
}

void GappedFragmentGraph::addEdge(size_t from, size_t to, long distance, size_t kmer_num, size_t read_num) {
	
    for (EdgeList::iterator it = _nodelist[from].begin(); it != _nodelist[from].end(); ++it) {
        if (it->component_id == to) {
            it->distances.push_back( static_cast<double> (distance));
            it->distance += distance;
            it->kmer_cov += kmer_num;
            it->read_cov += read_num;
            return;
        }
    }
    Edge e(to, distance, kmer_num, read_num, 0);
    _nodelist[from].push_back(e);
}

void GappedFragmentGraph::scoreAndRemoveNoise(const ComponentList& components) {
    // remove noise links based on possion distribution

    Scorer scorer(_K, INSERT_SIZE, DELTA, PAIR_KMER_NUM, GENOME_LEN);
	std::vector< double > score_list;

    // scoring and mark repeate
    {
        size_t graph_size = 0, count = 0, multiPeakCount = 0;
        std::map<size_t, int> repeateNum;
        std::vector< std::pair<size_t, size_t> > repeatePair;
        for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i, ++count) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if (j->kmer_cov < _pair_kmer_cutoff || j->read_cov < _pair_read_cutoff){
                    j = i->erase(j);
                    continue;
                }
                if( j->distances.size() > 100 && DELTA > 0.01) {
                    std::vector< Gaussian > gau;
                    GMM gmm;
                    gmm.gmm(j->distances, 3, gau, DELTA);
                    if(gau.size() > 1) {
                        double maxWeight = 0, meanDis = 0;
                        for(int k = 0; k < gau.size(); ++k) {
                            if(gau[k].weight > maxWeight) {
                                maxWeight = gau[k].weight;
                                meanDis = gau[k].mean;
                            }
                        }
                        j->distance = static_cast<long>(meanDis * j->kmer_cov);

                        repeatePair.push_back( std::make_pair(count, j->component_id) );
                        ++repeateNum[ count ];
                        ++repeateNum[ j->component_id ];
                        ++multiPeakCount;
                        /*std::ofstream os(boost::str(boost::format("distances.%d") % graph_size));
                        for(int i = 0; i < gau.size(); ++i) {
                            std::cout << gau[i] << std::endl;
                            os << "#" << gau[i] << std::endl;
                        }
                        for(int k = 0; k < j->distances.size(); ++k) {
                            os << j->distances[k] << std::endl;
                        }*/
                    }
                }
                //j->distance /= j->kmer_cov;
                size_t leni = components[count].length();
                size_t lenj =components[j->component_id].length();
                j->score = scorer.score(leni, lenj, j->distance / j->kmer_cov - leni, j->kmer_cov);
                score_list.push_back(j->score);
                ++j;
                ++graph_size;
            }
        }
        LOG4CXX_DEBUG(logger, boost::format("gaussian mixture model multiPeakCount=%d") % multiPeakCount);
        int repeateCount = 0;
        for(std::map<size_t, int>::iterator it = repeateNum.begin(); it != repeateNum.end(); ++it) {
            if( it->second > 1) {
                _repeateList.insert( it->first );
                ++repeateCount;
            }
        }
        LOG4CXX_DEBUG(logger, boost::format("find repeate num=%d") % repeateCount);
        LOG4CXX_DEBUG(logger, boost::format("befor remove nosize graph_size=%d") % graph_size);
    }

    std::sort(score_list.begin(), score_list.end());
    double threshold = score_list.empty() ? 0 : score_list[(size_t)(score_list.size() * _percent)];
    LOG4CXX_DEBUG(logger, boost::format("removeNoise threshold=%g") % threshold);

    // removeNoise
    {
        size_t graph_size = 0;

        for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if (j->score < threshold) {
                    j = i->erase(j);
                    continue;
                }
                ++j;
                ++graph_size;
            }
        }
        LOG4CXX_DEBUG(logger, boost::format("after remove nosize graph_size=%d") % graph_size);
    }
}

void GappedFragmentGraph::outputLP(std::ostream& os) {
    
	for (size_t i = 0; i < _nodelist.size(); ++i) {
        if( _repeateList.find(i) != _repeateList.end() ) {
            continue;
        }
		os << boost::format("var x_%d;") % i << std::endl;
		for (EdgeList::const_iterator it = _nodelist[i].begin(); it != _nodelist[i].end(); ++it) {

            if( _repeateList.find(it->component_id) != _repeateList.end() ) {
                continue;
            }
			os << boost::format("var e_%d_%d;") % i % it->component_id << std::endl; 
			os << boost::format("var E_%d_%d;") % i % it->component_id << std::endl;
		}
	}

	os << std::endl;
	os << "minimize z:  ";

	size_t count = 0;
	for (size_t i = 0; i < _nodelist.size(); ++i) {
        if( _repeateList.find(i) != _repeateList.end() ) {
            continue;
        }
        for(EdgeList::const_iterator it = _nodelist[i].begin(); it != _nodelist[i].end(); ++it) {
            if( _repeateList.find(it->component_id) != _repeateList.end() ) {
                continue;
            }
		    os << boost::format(" E_%d_%d + ") % i % it->component_id;
		    ++count;
		    if (count % 10 == 0)
			    os << std::endl;
        }
	}
	os << "0;" << std::endl << std::endl;

	size_t index = 1;
	for (size_t i = 0; i < _nodelist.size(); ++i) {
        if( _repeateList.find(i) != _repeateList.end() ) {
            continue;
        }
		for (EdgeList::const_iterator it = _nodelist[i].begin(); it != _nodelist[i].end(); ++it) {
            if( _repeateList.find(it->component_id) != _repeateList.end() ) {
                continue;
            }
			os << boost::format("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;") % index++ % it->component_id % i % i % it->component_id % (it->distance / it->kmer_cov) << std::endl;
		}
	}

	for (size_t i = 0; i < _nodelist.size(); ++i) {
        if( _repeateList.find(i) != _repeateList.end() ) {
            continue;
        }
		for (EdgeList::iterator it = _nodelist[i].begin(); it != _nodelist[i].end(); ++it) {
            if( _repeateList.find(it->component_id) != _repeateList.end() ) {
                continue;
            }
			os << boost::format("s.t. con%d : E_%d_%d + e_%d_%d >= 0;") % index++ % i % it->component_id % i % it->component_id << std::endl;
			os << boost::format("s.t. con%d : E_%d_%d - e_%d_%d >= 0;") % index++ % i % it->component_id % i % it->component_id << std::endl;
		}
	}

	os << std::endl;
    os << "end;";
}

std::ostream& operator<<(std::ostream& os, const GappedFragmentGraph& g) {
    size_t index = 0;
    for (GappedFragmentGraph::NodeList::const_iterator i = g._nodelist.begin(); i != g._nodelist.end(); ++i, ++index) {
        if( g._repeateList.find(index) != g._repeateList.end() ) {
            continue;
        }
        for (GappedFragmentGraph::EdgeList::const_iterator  j = i->begin(); j != i->end(); ++j) {
            if( g._repeateList.find(j->component_id) != g._repeateList.end() ) {
                continue;
            }
            if (j->distance > 0) {
                os << boost::format("%d\t%d\t%d\t%d\t%f") % index % j->component_id % (j->distance / j->kmer_cov) % j->kmer_cov % j->score << std::endl;
            }
        }
    }
    return os;
}
