// Link real APIs rather than the legacy configure macros' misdeclared exit().
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <sstream>
extern "C" {
#include <ViennaRNA/model.h>
}
int main() {
	std::stringstream data;
	boost::archive::binary_oarchive archive(data);
	archive << 42;
	boost::iostreams::filtering_ostream compressed;
	compressed.push(boost::iostreams::gzip_compressor());
	compressed.push(data);
	compressed << "RNA";
	compressed.reset();
	boost::program_options::options_description options("probe");
	options.add_options()("query", "RNA sequence");
	vrna_md_t model;
	vrna_md_set_default(&model);
	return boost::filesystem::exists("/") && boost::regex_match("RNA", boost::regex("RNA")) ? 0 : 1;
}
