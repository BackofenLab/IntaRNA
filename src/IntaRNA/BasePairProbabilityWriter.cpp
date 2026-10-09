#include "IntaRNA/BasePairProbabilityWriter.h"
#include <iomanip>
#include <sstream>
#include <iostream>
#include <memory>
#include <boost/iostreams/filtering_stream.hpp>

namespace IntaRNA {
std::string BasePairProbabilityWriter::block(const BasePairProbabilities & result,
		const RnaSequence & target,const RnaSequence & query,const std::string & sep)
{
	const bool empty=result.status()==BasePairProbabilities::Status::empty;
	if (!empty && result.status()!=BasePairProbabilities::Status::nonempty)
		throw std::logic_error("bpProb output requires successful explicit finalization");
	if (result.rawMasses().size1()!=target.size() || result.rawMasses().size2()!=query.size())
		throw std::invalid_argument("bpProb output sequence dimensions do not match result");
	// Validate every normalized value before opening a file or writing a header.
	const auto p=empty?Matrix<Z_type>():result.probabilities();
	std::ostringstream data;
	data<<std::setprecision(std::numeric_limits<Z_type>::max_digits10)<<"bpProb";
	for(size_t j=0;j<query.size();++j) data<<sep<<query.asString().at(j)<<'_'<<query.getInOutIndex(j);
	data<<'\n';
	for(size_t i=0;i<target.size();++i) {
		data<<target.asString().at(i)<<'_'<<target.getInOutIndex(i);
		for(size_t j=0;j<query.size();++j) { data<<sep;if(empty) data<<"NA";else data<<p(i,j); }
		data<<'\n';
	}
	return data.str();
}
void BasePairProbabilityWriter::emit(std::ostream & out,const std::string & data) {
	// Exceptions must leave the OpenMP structured block on the same thread.
	std::exception_ptr failure;
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		try {
			out<<data;out.flush();
			if(!out) throw std::runtime_error("base-pair probability output write failed");
		} catch (...) { failure=std::current_exception(); }
	}
	if(failure) std::rethrow_exception(failure);
}
void BasePairProbabilityWriter::write(std::ostream & out,const BasePairProbabilities & result,
		const RnaSequence & t,const RnaSequence & q,const std::string & sep) {
	emit(out,block(result,t,q,sep));
}
void BasePairProbabilityWriter::writeFile(const std::string & name,const BasePairProbabilities & result,
		const RnaSequence & t,const RnaSequence & q,const std::string & sep) {
	emitFile(name,block(result,t,q,sep));
}
void BasePairProbabilityWriter::emitFile(const std::string & name,const std::string & data) {
	std::ostream * out=newOutputStream(name);
	if(!out) throw std::runtime_error("could not open base-pair probability output '"+name+"'");
	// Standard streams are borrowed; file streams are owned here. Explicit reset
	// closes gzip before success is returned, exposing compression/close failures.
	std::unique_ptr<std::ostream> file(out==&std::cout || out==&std::cerr?nullptr:out);
	emit(*out,data);
	if(auto filtered=dynamic_cast<boost::iostreams::filtering_ostream *>(out)) filtered->reset();
}
}
