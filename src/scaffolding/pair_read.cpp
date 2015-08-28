#include "pair_read.h"

#include <numeric>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("scaffording.pair_read"));

PairReadSet::PairReadSet(std::istream& stream_1, std::istream& stream_2, size_t K, size_t INSERT_SIZE) : _stream_1(stream_1), _stream_2(stream_2), _k(K) , INSERT_SIZE(INSERT_SIZE) {
    DELTA = 0;
    PAIR_KMER_NUM = 0;
    if (!_stream_1 || !_stream_2) {
        LOG4CXX_INFO(logger,boost::format("can not find read files"));
    }
    readFile();
}

void PairReadSet::readFile() {
    if (_stream_1 && _stream_2) {
        DNASeqReader reader_1(_stream_1);
        DNASeqReader reader_2(_stream_2);
        DNASeq read_1;
        DNASeq read_2;
        PairRead p_r;
        while (reader_1.read(read_1) && reader_2.read(read_2)) {
            p_r.set1 = read_1.seq; 
            p_r.set2 = read_2.seq;
            _pair_reads.push_back(p_r);
        }
    }
}

void PairReadSet::buildConnectGraph(Graph& g, KmerTable& tbl, const std::vector<Component>& components) {
    std::string read1, read2;
    std::vector<int> insert_len_dis;
    LOG4CXX_INFO(logger, boost::format("in build graph: K=%lld") % _k);
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        read1 = pair_read.set1;
        read2 = pair_read.set2;
        if (read1.size() < _k || read2.size() < _k)
            continue;
        std::string read1_inv = make_complement(read1);
        std::string read2_inv = make_complement(read2);
        findLink(read1, read2_inv, g, tbl, components);
        findLink(read2, read1_inv, g, tbl, components);
    }
}


void PairReadSet::estimateInsertSize(const KmerTable &tbl) {
    std::string read1, read2;
    std::vector<int> insert_len_dis;
    BOOST_FOREACH(const PairRead& pair_read, _pair_reads) {
        if (insert_len_dis.size() >= 1000) {//trick
            break;
        }
        read1 = pair_read.set1;
        read2 = pair_read.set2;
        if (read1.size() < _k || read2.size() < _k)
            continue;
        std::string read1_inv = make_complement(read1);
        std::string read2_inv = make_complement(read2);

        LOG4CXX_TRACE(logger, boost::format("pair read %s %s") % read1 % read2_inv);
        LOG4CXX_TRACE(logger, boost::format("pair read %s %s") % read2 % read1_inv);
        estimateOnePR(read1, read2_inv, insert_len_dis, tbl);
        estimateOnePR(read2, read1_inv, insert_len_dis, tbl);
    }

    size_t sum = 0;
    size_t index = 0;
    BOOST_FOREACH(const int insert_num, insert_len_dis) {  
        if (insert_num > 10) {
            sum += insert_num;
            index ++;
        }
    }
    if (index > 0) {
        INSERT_SIZE = double(sum) / index;
        sum = 0;
        BOOST_FOREACH(const int insert_num, insert_len_dis) {
            if (insert_num > 10)
                sum += (insert_num - INSERT_SIZE)*(insert_num - INSERT_SIZE);
        }
        DELTA =  std::sqrt(double(sum) / index);
    }
    LOG4CXX_INFO(logger, boost::format("INSERT=%d DELTA=%d") % INSERT_SIZE % DELTA);
}

void PairReadSet::estimateOnePR(const std::string& read1, const std::string& read2, std::vector<int>& insert_len_dis, const KmerTable &tbl) {
    Kmer kmer1,kmer2;
    std::pair<size_t, long> left,right;
    for (int i=0; i<read1.size()-_k+1 && i<read2.size()-_k+1; i++) {
        kmer1 = read1.substr(i, _k);
        kmer2 = read2.substr(i, _k);

        left = tbl.findPos(kmer1);
        right = tbl.findPos(kmer2);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
        LOG4CXX_TRACE(logger, boost::format("INSERT SIZE pair kmre2=%s %d %d") % kmer2 % right.first % right.second);

        if (left.first != -1  && right.first != -1 && left.first == right.first) {
            if (right.second - left.second > INSERT_SIZE / 4 
                    && right.second - left.second < 2*INSERT_SIZE) 
                insert_len_dis.push_back(right.second - left.second);
        }
    }

}

struct PairCmp{
    bool operator()(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        return a.second < b.second;
    }
};

void PairReadSet::findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const std::vector<Component>& components) {
    std::set< std::pair<size_t, size_t>, PairCmp > find_component;
    for (int i=0; i<read1.size()-_k+1 && i<read2.size()-_k+1; i++) {
        Kmer kmer1 = read1.substr(i, _k);
        Kmer kmer2 = read2.substr(i, _k);

        std::pair<size_t, long> left = tbl.findPos(kmer1);
        std::pair<size_t, long> right = tbl.findPos(kmer2);
        
        if (left.first != right.first && left.first != -1 && right.first != -1) {
            PAIR_KMER_NUM ++;
            LOG4CXX_TRACE(logger, boost::format("pair kmre1=%s %d %d") % kmer1 % left.first % left.second);
            LOG4CXX_TRACE(logger, boost::format("pair kmre2=%s %d %d") % kmer2 % right.first % right.second);
        
            long len_tmp = INSERT_SIZE + left.second - right.second;
            long left_len = components[left.first].getLen();
            long overlap = left_len - len_tmp - _k + 1;
            if (overlap > 0) {
                if (overlap > INSERT_SIZE)
                    continue;
                len_tmp = left_len - _k + 1;
            }
            if (find_component.find(std::make_pair(left.first, right.first)) != find_component.end()) {
                graph.addEdge(left.first, right.first, len_tmp, false);
            } else {
                find_component.insert( std::make_pair(left.first, right.first));
                graph.addEdge(left.first, right.first, len_tmp, true);
            }
        }
    }
}

std::string PairReadSet::make_complement(std::string seq){
    std::reverse(seq.begin(), seq.end());
    static std::map< char, char > mapping = boost::assign::map_list_of
        ('A', 'T')
        ('C', 'G')
        ('G', 'C')
        ('T', 'A')
        ('N', 'N');

    for (size_t i = 0; i < seq.size(); ++i) {
        seq[i] = mapping[seq[i]];
    }
    return seq;
}

std::ostream& operator<<(std::ostream& os, const PairReadSet& p_r) {
    for (int i=0; i<p_r._pair_reads.size(); i++) {
        os << p_r._pair_reads[i].set1 <<"\t"<<p_r._pair_reads[i].set2 << std::endl; 
    }
    return os;
}


