#include "zanalyzer_config.h"
#include "json11.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<std::string> ZAnalyzerConfig::validate() const {
    std::vector<std::string> errors;

    if (input_file.empty())
        errors.push_back("input_file is required");
    if (z_min >= z_max)
        errors.push_back("zmin must be less than zmax");

    return errors;
}

ZAnalyzerConfig ZAnalyzerConfig::from_json(const std::string& filepath, std::vector<std::string>& errors) {
    ZAnalyzerConfig config;

    std::ifstream file(filepath.c_str());
    if (!file.is_open()) {
        errors.push_back("cannot open config file: " + filepath);
        return config;
    }

    std::stringstream buf;
    buf << file.rdbuf();
    std::string contents = buf.str();

    std::string parse_err;
    json11::Json json = json11::Json::parse(contents, parse_err);
    if (!parse_err.empty()) {
        errors.push_back("JSON parse error: " + parse_err);
        return config;
    }

    if (!json.is_object()) {
        errors.push_back("config file must contain a JSON object");
        return config;
    }

    for (auto& kv : json.object_items()) {
        const std::string& key = kv.first;
        const json11::Json& val = kv.second;

        if (key == "input_file")
            config.input_file = val.string_value();
        else if (key == "target_length")
            config.target_length = val.int_value();
        else if (key == "tolerance")
            config.tolerance = val.int_value();
        else if (key == "zmin")
            config.z_min = val.number_value();
        else if (key == "zmax")
            config.z_max = val.number_value();
        else if (key == "q")
            config.q = val.int_value();
        else if (key == "steps_between_samples")
            config.steps_between_samples = val.int_value();
        else if (key == "warmup")
            config.warmup = val.int_value();
        else if (key == "seed")
            config.seed = val.int_value();
        else
            std::cerr << "Warning: unknown config key '" << key << "' (ignored)" << std::endl;
    }

    return config;
}
