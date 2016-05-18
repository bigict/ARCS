#include "gapped_fragment_graph.h"
#include "gmm.h"

#include <cmath>
#include <limits>

#include <fstream>
#include <algorithm>
#include <map>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/math/distributions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GappedFragmentGraph"));

class Scorer {
public:
    Scorer(size_t K, size_t insert_size, double delta, size_t pair_kmer_num, size_t genome_len, size_t read_len) : _K(K), _pair_kmer_num(pair_kmer_num), _genome_len(genome_len), _gaussian(insert_size, delta), _insert_size(insert_size), _read_len(read_len) {
    }
    double score(size_t leni, size_t lenj, long gap, size_t c) const {
        BOOST_ASSERT(leni >= _K - 1);
        BOOST_ASSERT(lenj >= _K - 1);
        leni = std::min(leni, 2*_insert_size);
        lenj = std::min(lenj, 2*_insert_size);
        double estimate_num = 0;
        for(size_t i = 0; i < leni + _read_len - 2 * _K; ++i) {
            estimate_num += gaussian(_read_len + leni + gap + lenj - 2 * _K - i) - gaussian(leni + gap - i);
            //sum += norm.cdf(read_len+l1+d+l2-2*K - i, mu, sigma) - norm.cdf(l1+d-i, mu, sigma)
        }
        return estimate_num;
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
    size_t _insert_size;
    size_t _read_len;
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
            if(kmer_num == 0) {
                it->read_cov += read_num;
                return;
            }
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

void GappedFragmentGraph::scoreAndRemoveNoise(const ComponentList& components, size_t read_len) {

    // remove low kmer_cov, read_cov edge and mark repeate using gaussion mix model
    size_t avgDistanceSize = 0, edgeNum = 0;
    for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        EdgeList::iterator j = i->begin();
        while (j != i->end()) {
            if (j->kmer_cov < _pair_read_cutoff * read_len && j->read_cov < _pair_read_cutoff){
                ++j;
                continue;
            }
            ++edgeNum;
            avgDistanceSize += j->distances.size();
            ++j;
        }
    }
    LOG4CXX_DEBUG(logger, boost::format("before remove nosize graph_size=%d") % edgeNum);
    avgDistanceSize /= edgeNum;
    LOG4CXX_DEBUG(logger, boost::format("avgDistanceSize = %d") % avgDistanceSize);

    {
        size_t graph_size = 0, count = 0, multiPeakCount = 0;
        std::map<size_t, int> repeateNum;
        std::vector< std::pair<size_t, size_t> > repeatePair;
        for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i, ++count) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if (j->kmer_cov < _pair_read_cutoff * (read_len - _K + 1) && j->read_cov < _pair_read_cutoff){
                    j = i->erase(j);
                    continue;
                }
                std::sort(j->distances.begin(), j->distances.end());
                if( j->distances.size() > std::max(100.0, 1.0*avgDistanceSize) && DELTA > 0.01) {
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
                        if(components[count].length() <= read_len) {
                            ++repeateNum[ count ];
                        } 
                        if(components[j->component_id].length() <= read_len) {
                            ++repeateNum[ j->component_id ];
                        } 
                        ++multiPeakCount;
                    }
                }
                size_t leni = components[count].length();
                size_t lenj =components[j->component_id].length();
                j->score = j->read_cov;
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
    }
    //remove tips and low quality edge
    int loop = 2;
    for(int i = 0; i < loop; ++i) {
        edgeNum = 0;
        size_t count = 0;
        ReverseMap reverseMap;
        for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i, ++count) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if (_nodelist[j->component_id].size() == 0 && deepMoreThan2(count)){
                    LOG4CXX_DEBUG(logger, boost::format("remove tail tip from %d to %d") % count % j->component_id);
                    j = i->erase(j);
                    continue;
                }
                // for large degree node, remove its low quality edges 
                if (j->read_cov < 3 * _pair_read_cutoff && std::count_if(i->begin(), i->end(), [j](const Edge& e){return e.read_cov > j->read_cov;}) >= 2) {
                    LOG4CXX_DEBUG(logger, boost::format("remove low quality edge from %d to %d") % count % j->component_id);
                    j = i->erase(j);
                    continue;
                }
                reverseMap[j->component_id].insert(std::make_pair(count, j->read_cov));
                ++j;
            }
        }
        count = 0;
        for (NodeList::iterator i = _nodelist.begin(); i != _nodelist.end(); ++i, ++count) {
            EdgeList::iterator j = i->begin();
            while (j != i->end()) {
                if ((reverseMap.find(count) == reverseMap.end() || reverseMap.find(count)->second.size() == 0) && reverseDeepMoreThan2(j->component_id, reverseMap)){
                    LOG4CXX_DEBUG(logger, boost::format("remove head tip from %d to %d") % count % j->component_id);
                    j = i->erase(j);
                    continue;
                }
                if (j->read_cov < 3 * _pair_read_cutoff && 
                        std::count_if(reverseMap[j->component_id].begin(), reverseMap[j->component_id].end(), [j](const std::pair<size_t, size_t> &e){return e.second > j->read_cov;}) >= 2) {
                    LOG4CXX_DEBUG(logger, boost::format("remove low quality edge from %d to %d") % count % j->component_id);
                    j = i->erase(j);
                    continue;
                }
                ++edgeNum;
                ++j;
            }
        }
    }
    LOG4CXX_DEBUG(logger, boost::format("after remove nosize graph_size=%d") % edgeNum);
}

bool GappedFragmentGraph::deepMoreThan2(size_t startComp) {
    for(EdgeList::iterator it = _nodelist[startComp].begin(); it != _nodelist[startComp].end(); ++it) {
        if( _nodelist[it->component_id].size() > 0)
            return true;
    }
    return false;
}

bool GappedFragmentGraph::reverseDeepMoreThan2(size_t startComp, ReverseMap &reverseMap) {
    for(typename ReverseMapEle::iterator it = reverseMap[startComp].begin(); it != reverseMap[startComp].end(); ++it) {
        if(reverseMap.find(it->first) != reverseMap.end() && reverseMap.find(it->first)->second.size() > 0)
            return true;
    }
    return false;
}


void GappedFragmentGraph::outputLP(std::ostream& os, const ComponentList& components) {
    
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
            // remove the top 10% and tail 10% then calculate average
            // if change, you must change the operator<< function too!!!
            std::pair<double, size_t> sum_num = std::accumulate(it->distances.begin() + (int)std::floor(it->distances.size() * 0.1), it->distances.begin() + (int)std::ceil(it->distances.size() * 0.9), std::make_pair(0.0, 0), [](const std::pair<double, size_t> &x, double y){return std::make_pair(x.first+y, x.second+1);});

			os << boost::format("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;") % index++ % it->component_id % i % i % it->component_id % (size_t)(sum_num.first / sum_num.second) /*it->distances[ it->distances.size() / 2 ]*/ /*(it->distance / it->kmer_cov)*/ << std::endl;
            //maintain the direction information of edge, but if circle, the solution will be all 0....
            //os << boost::format("s.t. con%d : x_%d >= x_%d + %d;") % index++ % it->component_id % i % std::max(0.0, components[i].length() - (INSERT_SIZE + 3 * DELTA)) << std::endl;
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
                // remove the top 10% and tail 10% then calculate average
                // if change, you must change the outputLP function too!!!
                std::pair<double, size_t> sum_num = std::accumulate(j->distances.begin() + (int)std::floor(j->distances.size() * 0.1), j->distances.begin() + (int)std::ceil(j->distances.size() * 0.9), std::make_pair(0.0, 0), [](const std::pair<double, size_t> &x, double y){return std::make_pair(x.first+y, x.second+1);});
                os << boost::format("%d\t%d\t%d\t%d\t%f") % index % j->component_id % (size_t)(sum_num.first / sum_num.second)/*j->distances[ j->distances.size() / 2 ]*/ /*(j->distance / j->kmer_cov)*/ % j->kmer_cov % j->score << std::endl;
            }
        }
    }
    return os;
}
