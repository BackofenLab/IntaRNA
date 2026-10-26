#include "ParameterFile.h"
#include "nlohmann/json.hpp"

#include <boost/filesystem.hpp>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
using Json = nlohmann::json;
using namespace boost::program_options;

// Read once, remove a UTF-8 BOM, then detect JSON independently of the suffix.
std::string readText(std::istream& input)
{
	std::string text{std::istreambuf_iterator<char>(input), {}};
	if (input.bad()) throw std::runtime_error("Could not read parameter file");
	if (text.starts_with("\xef\xbb\xbf")) text.erase(0, 3);
	return text;
}

bool isJson(const std::string& text)
{
	const auto first = text.find_first_not_of(" \t\r\n");
	return first != std::string::npos && (text[first] == '{' || text[first] == '[');
}

Json readJson(const std::string& text)
{
	// Reject duplicate keys and excessive nesting rather than silently replacing
	// settings or consuming an unbounded call stack on malformed input.
	std::vector<std::set<std::string>> keys;
	auto callback = [&keys](int depth, Json::parse_event_t event, Json& value) {
		if (depth > 32) throw std::runtime_error("JSON parameter nesting exceeds 32 levels");
		if (event == Json::parse_event_t::object_start) keys.emplace_back();
		if (event == Json::parse_event_t::key && !keys.back().insert(value.get<std::string>()).second)
			throw std::runtime_error("Duplicate JSON key: " + value.get<std::string>());
		if (event == Json::parse_event_t::object_end) keys.pop_back();
		return true;
	};
	Json document = Json::parse(text, callback);
	if (!document.is_object()) throw std::runtime_error("JSON parameter file must be an object");
	for (const auto& [key, value] : document.items()) {
		if (key != "parameters" && key != "schema" && key != "_comment" && key != "gui")
			throw std::runtime_error("Unknown JSON document field: " + key);
	}
	if (document.value("schema", std::string()) != "IntaRNA-parameters-v1")
		throw std::runtime_error("Expected schema IntaRNA-parameters-v1");
	if (!document.contains("parameters") || !document["parameters"].is_object())
		throw std::runtime_error("JSON parameters must be an object");
	return document["parameters"];
}

std::string scalar(const Json& value)
{
	if (value.is_string()) {
		const auto text = value.get<std::string>();
		if (text.find('\0') != std::string::npos) throw std::runtime_error("NUL in JSON parameter text");
		return text;
	}
	if (value.is_boolean() || value.is_number()) return value.dump();
	throw std::runtime_error("JSON parameter values must be strings, numbers, booleans or embedded text");
}
}

ParameterFile::~ParameterFile()
{
	// Cleanup must also succeed during exception unwinding and after failed input.
	if (!directory.empty()) {
		boost::system::error_code ignored;
		boost::filesystem::remove_all(directory, ignored);
	}
}

std::string ParameterFile::writeText(const std::string& name, const std::string& text)
{
	// Never interpret embedded data or JSON keys as filesystem paths.
	if (directory.empty()) {
		const auto path = boost::filesystem::temp_directory_path()
			/ boost::filesystem::unique_path("intarna-json-%%%%-%%%%-%%%%-%%%%");
		if (!boost::filesystem::create_directory(path)) throw std::runtime_error("Could not create embedded-input directory");
		directory = path.string();
		boost::filesystem::permissions(path, boost::filesystem::owner_all);
	}
	const auto path = (boost::filesystem::path(directory) / name).string();
	std::ofstream file(path, std::ios::binary);
	file.exceptions(std::ios::badbit | std::ios::failbit);
	file.write(text.data(), text.size());
	file.close();
	return path;
}

boost::program_options::parsed_options
ParameterFile::parse(std::istream& input, const boost::program_options::options_description& options)
{
	const auto text = readText(input);
	if (!isJson(text)) {
		std::istringstream legacy(text);
		return boost::program_options::parse_config_file(legacy, options);
	}
	const auto parameters = readJson(text);
	boost::program_options::parsed_options parsed(&options);
	const std::set<std::string> fileInputs = {"query", "target", "qAccFile", "tAccFile", "qShape", "tShape", "energyVRNA"};
	for (const auto& [key, value] : parameters.items()) {
		// Reuse Program_options type conversion and validation without converting
		// through config-file syntax (which would reinterpret # or newlines).
		const auto* definition = options.find_nothrow(key, false);
		if (!definition || definition->long_name() != key || key == "parameterFile")
			throw std::runtime_error("Unknown or recursive JSON parameter: " + key);
		auto append = [&](const Json& item) {
			boost::program_options::option option;
			option.string_key = key;
			if (definition->semantic()->max_tokens() == 0) {
				if (!item.is_boolean()) throw std::runtime_error(key + " requires a JSON boolean");
				if (!item.get<bool>()) return;
			} else if (item.is_object()) {
				if (!fileInputs.contains(key)) throw std::runtime_error(key + " does not accept embedded files");
				if (item.size() == 1 && item.contains("text") && item["text"].is_string()) {
					option.value.push_back(writeText(key + ".txt", scalar(item["text"])));
				} else if ((key == "qAccFile" || key == "tAccFile") && item.size() == 1 && item.contains("texts") && item["texts"].is_array() && !item["texts"].empty()) {
					// Multi-sequence accessibility follows the existing -sN naming
					// convention; indices refer to original FASTA record numbers.
					const auto& texts = item["texts"];
					for (size_t i = 0; i < texts.size(); ++i) {
						if (!texts[i].is_string()) throw std::runtime_error(key + ": texts must contain strings");
						writeText(key + "-s" + std::to_string(i + 1) + ".txt", scalar(texts[i]));
					}
					option.value.push_back((boost::filesystem::path(directory) / (key + ".txt")).string());
				} else throw std::runtime_error(key + ": expected {\"text\":\"...\"} (or accessibility {\"texts\":[...]})");
			} else option.value.push_back(scalar(item));
			parsed.options.push_back(std::move(option));
		};
		if (value.is_array()) {
			if (!definition->semantic()->is_composing() || value.empty()) throw std::runtime_error(key + " does not accept an empty or repeated-value array");
			for (const auto& item : value) append(item);
		} else append(value);
	}
	return parsed;
}

std::string ParameterFile::personality(const std::string& filename)
{
	std::ifstream input(filename, std::ios::binary);
	if (!input) throw std::runtime_error("Could not open parameter file: " + filename);
	const auto text = readText(input);
	if (isJson(text)) {
		const auto parameters = readJson(text);
		if (!parameters.contains("personality")) return {};
		if (!parameters["personality"].is_string()) throw std::runtime_error("personality must be a string");
		return scalar(parameters["personality"]);
	}
	// Legacy parameter files receive the same personality/default precedence.
	boost::program_options::options_description options;
	options.add_options()("personality", boost::program_options::value<std::string>());
	std::istringstream legacy(text);
	boost::program_options::variables_map values;
	store(boost::program_options::parse_config_file(legacy, options, true), values);
	return values.count("personality") ? values["personality"].as<std::string>() : "";
}
