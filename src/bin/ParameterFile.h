#ifndef INTARNA_PARAMETERFILE_H_
#define INTARNA_PARAMETERFILE_H_

#include <boost/program_options.hpp>
#include <istream>
#include <string>

/** Reads legacy or JSON CLI parameters and owns any embedded input files. */
class ParameterFile {
public:
	ParameterFile() = default;
	~ParameterFile();
	ParameterFile(const ParameterFile&) = delete;
	ParameterFile& operator=(const ParameterFile&) = delete;

	/**
	 * Parse parameters through the existing option definitions. Embedded text
	 * lives in private temporary files until this object is destroyed.
	 * @param input parameter file stream (format detected from its contents)
	 * @param options supported CLI options, including their validators
	 * @return parsed options for boost::program_options::store
	 * @throws std::exception on malformed parameters or failed file operations
	 */
	boost::program_options::parsed_options
	parse(std::istream& input, const boost::program_options::options_description& options);

	/**
	 * Read a file's personality before constructing personality-specific defaults.
	 * @param filename legacy or JSON parameter file
	 * @return personality name, or empty if not specified
	 */
	static std::string personality(const std::string& filename);

private:
	//! private directory containing generated filenames only; never user paths
	std::string directory;

	/** Materialize text with a generated basename; return its absolute path. */
	std::string writeText(const std::string& name, const std::string& text);
};

#endif
