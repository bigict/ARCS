#include "gap_filler.h"
#include "condensed_debruijn_graph_reader.h"
#include "constant.h"
#include "contigs.h"
#include "utils.h"

#include <deque>
#include <fstream>
#include <iterator>
#include <numeric>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GapFiller"));

//gap filling
void GapFiller::fill() {
    LOG4CXX_DEBUG(logger, "Begin Gap Filling ... ");

    size_t num_failed_gaps = 0, num_uniq_gaps = 0, num_multi_gaps = 0, num_overlaps = 0;
	for (size_t i = 0; i < _scaffolds.size(); ++i) {
		if (_scaffolds[i].contigs.size() <= 1) { // no gaps
			continue;
		}
        Component::ContigIdList contigs = _scaffolds[i].contigs;

		for (size_t j = 1; j < contigs.size(); ++j) {
			int left_index = contigs[j - 1];
			int right_index = contigs[j];

            const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
            const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
            std::string suffix = lsequence.substr(lsequence.length() - (_K - 1), _K - 1);
            std::string prefix = rsequence.substr(0, _K - 1);

			int overlap = alignment(suffix, prefix);
            int gap = _scaffolds[i].gaps[j - 1];
			if (overlap >= _OVERLAP) {
                _gapinfo_tbl[std::make_pair(i, j)] = GapInfo(-1, -overlap);
                ++num_overlaps;
			} else if (gap > _INSERT_SIZE + 3*_DELTA) {
                _gapinfo_tbl[std::make_pair(i, j)] = GapInfo(-1, gap);
                ++num_failed_gaps;
            } else {
                try {
                    GapInfo gapinfo(-1, gap);
                    BFS(left_index, right_index, gap, &gapinfo);
                    _gapinfo_tbl[std::make_pair(i, j)] = gapinfo;
                    if (gapinfo.pathlist.empty()) {
                        ++num_failed_gaps;
                    } else if (gapinfo.pathlist.size() == 1) {
                        ++num_uniq_gaps;
                    } else {
                        ++num_multi_gaps;
                    }
                } catch(std::bad_alloc) {
                    LOG4CXX_WARN(logger, "BFS is too memory-intensive, ignoring...");
                    _gapinfo_tbl[std::make_pair(i, j)] = GapInfo(-1, -gap);
                    ++num_failed_gaps;
                }
            }
		}
	}
	LOG4CXX_DEBUG(logger, boost::format("The number of gaps = %d") % _gapinfo_tbl.size());
	LOG4CXX_DEBUG(logger, boost::format("The number of failed gaps = %d") % num_failed_gaps);
	LOG4CXX_DEBUG(logger, boost::format("The number of unique gaps = %d") % num_uniq_gaps);
	LOG4CXX_DEBUG(logger, boost::format("The number of multiple gaps = %d") % num_multi_gaps);
	LOG4CXX_DEBUG(logger, boost::format("The number of overlap = %d") % num_overlaps);
}

//Align two nerghboring contigs
size_t GapFiller::alignment(const std::string& suffix, const std::string& prefix) {
	BOOST_ASSERT(prefix.length() == suffix.length());

    std::vector< size_t > score(prefix.size() + 1, 0);

	for (size_t i = 1; i <= suffix.size(); ++i) {
		for (size_t j = prefix.size(); j > 0; --j) {
			if (prefix[j - 1] == suffix[i - 1]) {
				score[j] = score[j - 1] + 1;
			} else {
				score[j] = 0;
			}
		}
	}

	size_t ret = 0;
	for (size_t j = 0; j <= prefix.size(); ++j) {
		if (score[j] == j) {
			if (ret < score[j]) {
				ret = score[j];
			}
		}
	}
	return ret;
}

struct _ContigLengthCmp_ {
    bool operator()(const std::string& c1, const std::string& c2) const {
        return c1.length() < c2.length();
    }
};

bool GapFiller::input_scaffold(const std::string& file) {
    std::ifstream stream(file.c_str());
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
        return false;
    }

    ComponentReader reader(stream);
    Component component;
    while (reader.read(component)) {
        if (!component.gaps.empty()) {
            component.gaps.pop_back();
        }
        _scaffolds.push_back(component);
    }
    return true;
}

bool GapFiller::input_contigs(const std::string& file) {
    // cdbg_copy_number.fa
    std::ifstream stream(file.c_str());
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
        return false;
    }
    ContigReader reader(stream);
    Contig contig;
    while (reader.read(contig)) {
        _uniq_graph.addEdge(contig.seq, contig.copy_num);
    }

    LOG4CXX_DEBUG(logger, boost::format("all contig number=%d") % _uniq_graph._indexer.size());

    return true;
}

bool GapFiller::input_debruijn(const std::string& file) {
    // condensed_de_bruijn_graph_before_trimming.data
    std::ifstream stream(file.c_str());
    if (!stream) {
        LOG4CXX_WARN(logger, boost::format("%s: No such file or directory!!!") % file);
        return false;
    }
    CondensedDeBruijnGraphReader reader(stream);
    CondensedDeBruijnGraphEdge debruijn;
    while (reader.read(debruijn)) {
        _all_graph.addEdge(debruijn.seq, 0);
    }
    LOG4CXX_DEBUG(logger, boost::format("All condensed edges number=%d") % _all_graph._indexer.size());

    return true;
}

std::string GapFiller::path2seq(const CondensedDeBruijnGraph& graph, const Path& path) const {
    return path2seq(graph, path, 0, path.size());
}

std::string GapFiller::path2seq(const CondensedDeBruijnGraph& graph, const Path& path, size_t i, size_t j) const {
    std::string seq;
    if (path.size() > i) {
        seq += graph._indexer[path[i++]].seq;
        while (i < j) {
            seq += graph._indexer[path[i++]].seq.substr(_K - 1);
        }
    }
    return seq;
}

void GapFiller::BFS(const CondensedDeBruijnGraph& graph, const std::string& lseq, const std::string& rseq, int gap, size_t max_depth, size_t max_queue, PathList& pathlist) {

    std::string suffix = lseq.substr(lseq.length() - (_K - 1), _K - 1);
    std::string prefix = rseq.substr(0, _K - 1);

    CondensedDeBruijnGraph::NodeList::const_iterator i = graph._children.find(suffix), j = graph._parents.find(prefix);
    if (i != graph._children.end() && j != graph._parents.end()) {
        for (CondensedDeBruijnGraph::EdgeList::const_iterator m = i->second.begin(); m != i->second.end(); ++m) {
            if (m->first == lseq[lseq.length() - _K]) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator n = j->second.begin(); n != j->second.end(); ++n) {
                    if (n->first == rseq[_K - 1]) {
                        BFS(graph, m->second, n->second, gap, _STEP, 2000, pathlist);
                    }
                }
            }
        }
    }
}

void GapFiller::BFS(const CondensedDeBruijnGraph& graph, const size_t from, const size_t to, int gap, size_t max_depth, size_t max_queue, PathList& pathlist) {
	size_t step = 0;

    Path path = boost::assign::list_of(from);
    typedef std::pair< Path, int > Choice;
    std::deque< Choice > Q = boost::assign::list_of(std::make_pair(
                path, 0 - (_K - 1)
                ));
    while (!Q.empty()) {
        ++step;
		if (step > max_depth) {
			LOG4CXX_DEBUG(logger, boost::format("BFS step is bigger than %d...") % max_depth);
			break;
		}

        if (Q.size() > max_queue) {
			LOG4CXX_DEBUG(logger, boost::format("BFS Q size=%d is greater than %d...") % Q.size() % max_queue);
            break;
        }

        std::deque< Choice > nextQ;
        while (!Q.empty()) {
            Choice choice = Q.front();
            Q.pop_front();

            size_t node = choice.first.back();
            if (node == to) {
                pathlist.push_back(choice.first);
            }
            if (choice.second > gap - (_K - 1)) {
                continue;
            }
            const CondensedDeBruijnGraph::CondensedEdge& edge = graph._indexer[node];
            std::string suffix = edge.seq.substr(edge.seq.length() - (_K - 1), _K - 1);
            CondensedDeBruijnGraph::NodeList::const_iterator i = graph._parents.find(suffix);
            if (i != graph._parents.end()) {
                for (CondensedDeBruijnGraph::EdgeList::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
                    //size_t length = graph._indexer[j->second].seq.length();
                    size_t length = graph._indexer[j->second].seq.length() - (_K - 1);

                    Choice new_choice(choice);
                    new_choice.first.push_back(j->second);
                    new_choice.second += length;

                    nextQ.push_back(new_choice);
                }
            }
        }

        Q = nextQ;
    }
}

// Run BFS to fill current gap 
void GapFiller::BFS(size_t left_index, size_t right_index, int gap, GapInfo* gapinfo) {
	int dis = gap + (_K-1+3*_DELTA) +30; //gap constraints

    PathList pathlist;
    BFS(_uniq_graph, left_index, right_index, dis, _STEP, 1000, pathlist);
    if (pathlist.empty()) {
        const std::string& lsequence = _uniq_graph._indexer[left_index].seq;
        const std::string& rsequence = _uniq_graph._indexer[right_index].seq;
        BFS(_all_graph, lsequence, rsequence, dis, _STEP, 2000, pathlist);
    }

    if (!pathlist.empty()) {
        if (pathlist.size() <= MAX_CHOICE) {
            if (gapinfo != NULL) {
                gapinfo->graph = 0;
                gapinfo->pathlist = pathlist;
            }
        } else {
            if (gapinfo != NULL) {
                gapinfo->graph = -1;
            }
        }
    }
}

std::ostream& operator<<(std::ostream& os, const GapFiller &obj) {
	LOG4CXX_DEBUG(logger, "Output each scaffold segments begin...");

    size_t num_failed_gaps = 0;
    std::vector< size_t > scaffold_length_list;

    for (size_t i = 0; i < obj._scaffolds.size(); ++i) {
        const Component::ContigIdList& contigs = obj._scaffolds[i].contigs;
        BOOST_ASSERT(!contigs.empty());

        std::string seq = obj._uniq_graph._indexer[contigs[0]].seq;
        for (size_t j = 1; j < contigs.size(); ++j) {
            size_t left_index = contigs[j - 1];
            size_t right_index = contigs[j];

            GapFiller::GapInfoTable::const_iterator it = obj._gapinfo_tbl.find(std::make_pair(i, j));
            BOOST_ASSERT(it != obj._gapinfo_tbl.end());

            if (it->second.graph == -1 && it->second.gap < 0) {     // overlap
                BOOST_ASSERT(seq.length() >= -it->second.gap);
                seq.resize(seq.length() + it->second.gap);
            } else if (it->second.graph == -1) {                    // failed gap
                BOOST_ASSERT(it->second.gap >= 0);
                //for failed gap , the count of 'N' is determined by the distance estimation
                seq += std::string(it->second.gap, 'N');
                ++num_failed_gaps;
            } else if (it->second.graph != -1 && it->second.pathlist.size() == 1) { // unique
                const GapFiller::Path& path = it->second.pathlist[0];
                BOOST_ASSERT(path.size() >= 2);
                std::string gap = obj.path2seq(it->second.graph == 0 ? obj._uniq_graph : obj._all_graph, path, 1, path.size() - 1);
                // remove overlap kmer from two sides
                BOOST_ASSERT(gap.length() >= obj._K - 1);
                seq += gap.substr(obj._K - 1);
                BOOST_ASSERT(seq.length() >= obj._K - 1);
                seq.resize(seq.length() - (obj._K - 1));
            } else {                                                // multi_gap
                if (it->second.gap > 0) {
                    seq += std::string(it->second.gap, 'N');
                }
            }

            seq += obj._uniq_graph._indexer[contigs[j]].seq;
        }
        os << boost::format(">scaf_%d_%d\n") % i %  seq.size();
        os << seq << '\n';
        scaffold_length_list.push_back(seq.length());
    }
    BOOST_ASSERT(obj._scaffolds.size() ==  scaffold_length_list.size());

    std::sort(scaffold_length_list.begin(), scaffold_length_list.end(), std::greater< size_t >());
    size_t GENOME = std::accumulate(scaffold_length_list.begin(), scaffold_length_list.end(), 0);

    LOG4CXX_INFO(logger, boost::format("Genome estimated size : %d") % GENOME);

	long N50 = GENOME * 0.5 ;
	long N90 = GENOME * 0.9 ;

	LOG4CXX_INFO(logger, "---Scaffold---");
	LOG4CXX_INFO(logger, "Scaff length : Top 20...");

	for (size_t i = 0; i < std::min((size_t)20, scaffold_length_list.size()); ++i) {
        LOG4CXX_INFO(logger, boost::format("%d") % scaffold_length_list[i]);
	}

    std::vector< size_t > sum (scaffold_length_list.size(), 0);
	
    if (!scaffold_length_list.empty()) {
        sum[0] = scaffold_length_list[0];
        for (size_t i = 1; i < scaffold_length_list.size(); ++i) {
            sum[i] +=  sum[i - 1] + scaffold_length_list[i];
        }
    }
	long N50_val = 0;
	long N90_val = 0;

	for (size_t i = 0 ; i < sum.size(); ++i) {
		if (sum[i] >= N50 && N50_val == 0) {
			N50_val = scaffold_length_list[i];
		}
		if (sum[i] >= N90 && N90_val == 0) {
			N90_val = scaffold_length_list[i];
			break;
		}
	}

	LOG4CXX_INFO(logger, boost::format("N50 : %d") % N50_val);
	LOG4CXX_INFO(logger, boost::format("N90 : %d") % N90_val);

	LOG4CXX_INFO(logger, boost::format("All gaps number in scaffolds : %d") % obj._gapinfo_tbl.size());
	LOG4CXX_INFO(logger, boost::format("filled gap numbers : %d") % (obj._gapinfo_tbl.size() - num_failed_gaps));
	LOG4CXX_INFO(logger, boost::format("failed gap numbers : %d") % num_failed_gaps);

	LOG4CXX_DEBUG(logger, "Output each scaffold segments end...");
	return os;
}
