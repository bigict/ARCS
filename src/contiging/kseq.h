#ifndef kseq_h_
#define kseq_h_

#include <string>
#include <iostream>


class DNASeq {
public:
    DNASeq() {}
    DNASeq(const std::string& name, const std::string& seq, const std::string& quality) : name(name), seq(seq), _quality(quality) {}
    virtual ~DNASeq() {}

    std::string name;
    std::string seq;
    int quality() const;
private:
    friend class DNASeqReader;
    friend std::ostream& operator << (std::ostream& os, const DNASeq& seq);

    std::string _quality;
};

class DNASeqReader {
public:
    DNASeqReader(std::istream& stream) : _stream(stream) {}
    
    bool read(DNASeq& sequence);
private:
    std::istream& _stream;
};

#endif // kseq_h_
