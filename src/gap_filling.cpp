#include "gap_filling.h"
#include "gap_filler.h"
#include "constant.h"

#include <string>
#include <sstream>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GapFilling"));
GapFilling GapFilling::_runner;

extern int EXTEND;
extern int STEP;
extern int K;
extern int OVERLAP;
extern int MU;
extern int var;

int EXTEND = -1;
int STEP = -1;
int K = -1;
int OVERLAP = -1;
int MU = -1;
int var = -1;

int GapFilling::run(const Properties& options, const Arguments& arguments) {
    int r = 0;
    if ((r = checkOptions(options)) != 0) {
        return r;
    }
    LOG4CXX_DEBUG(logger, "gap_filling begin");
		
    std::string initial_contig_file_name;
    std::string scaffold_file_name;
    std::string condensed_contig_file_name;
    std::string work_dir;
    K = options.get< size_t >("K", kKmerSize);
    work_dir = options.get< std::string >("d", ".");
    if (options.find("MAX_OVERLAP") != options.not_found()) {
        OVERLAP = options.get< size_t >("MAX_OVERLAP");
    }
    if (options.find("C") != options.not_found()) {
        condensed_contig_file_name = options.get< std::string >("C");
    }
    if (options.find("l") != options.not_found()) {
        scaffold_file_name = options.get< std::string >("l");
    }
    if (options.find("I") != options.not_found()) {
        initial_contig_file_name = options.get< std::string >("I");
    }
	

	if(chdir(work_dir.c_str()))
	{
		cout << "No such file or directory!! Change work directory failed!!!" << endl;
		exit(EXIT_FAILURE);
	}

    MU = options.get< size_t >("INSERT_SIZE");
    var = options.get< double >("DELTA");

	STEP = MU + 3 * var;
	EXTEND = 200;
	
	GapFiller gf;
    // load data
    {
        typedef bool(GapFiller::*LoadDataPtr)(const std::string&);
        typedef std::tr1::tuple< std::string, LoadDataPtr > File2FuncPtr;
        std::vector< File2FuncPtr > file2func_list = boost::assign::list_of
            (std::tr1::make_tuple(condensed_contig_file_name, &GapFiller::input_contigs))
            (std::tr1::make_tuple(initial_contig_file_name,   &GapFiller::input_debruijn))
            (std::tr1::make_tuple(scaffold_file_name,         &GapFiller::input_scaffold))
            ;
        BOOST_FOREACH(const File2FuncPtr& file2func, file2func_list) {
            std::string file = std::tr1::get< 0 >(file2func);
            if (!boost::bind(std::tr1::get< 1 >(file2func), &gf, _1)(file)) {
                LOG4CXX_ERROR(logger, boost::format("failed to load %s") % file);
                return 2;
            }
        }
    }

	gf.fill();

    // write data
    {
        std::string file = boost::str(boost::format("%dmer.scaf_seq_with_gaps") % K);
        std::ofstream stream(file.c_str());
        if (!stream) {
            LOG4CXX_ERROR(logger, boost::format("Create %s error!!!") % file);
            return 2;
        }
        stream << gf;
    }

    LOG4CXX_DEBUG(logger, "gap_filling end");

    return 0;
}

GapFilling::GapFilling() : Runner("c:s:d:K:O:C:I:l:h", boost::assign::map_list_of('O', "MAX_OVERLAP")) {
    RUNNER_INSTALL("gap_filling", this, "gap_filling");
}

int GapFilling::printHelps() const {
    std::cout << "arcs gap_filling -K [kmer] -O [overlpa_for_me] -C [condensed_contig_file_name] -I [initial_contig_file_name] -l [line_component_file]" << std::endl;
    return 256;
}

int GapFilling::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
