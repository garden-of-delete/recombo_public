#include "bfacf_config.h"
#include "json11.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<std::string> BfacfConfig::validate() const {
    std::vector<std::string> errors;

    if (input_file.empty())
        errors.push_back("input_file is required");
    if (output_file.empty())
        errors.push_back("output_file is required");
    if (z <= 0.0)
        errors.push_back("z must be > 0");
    if (steps <= 0)
        errors.push_back("steps (-n) must be > 0");
    if (save_interval < 0)
        errors.push_back("save_interval (-k) must be >= 0");
    if (save_interval > steps)
        errors.push_back("save_interval (-k) must be <= steps (-n)");
    if (output_format != "cube" && output_format != "text" && output_format != "newsud")
        errors.push_back("output_format (-fmt) must be one of: cube, text, newsud");

    return errors;
}

BfacfConfig BfacfConfig::from_json(const std::string& filepath, std::vector<std::string>& errors) {
    BfacfConfig config;

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
        else if (key == "output_file")
            config.output_file = val.string_value();
        else if (key == "z")
            config.z = val.number_value();
        else if (key == "q")
            config.q = val.int_value();
        else if (key == "steps")
            config.steps = val.int_value();
        else if (key == "save_interval")
            config.save_interval = val.int_value();
        else if (key == "block_file_size")
            config.block_file_size = val.int_value();
        else if (key == "seed")
            config.seed = val.int_value();
        else if (key == "output_format")
            config.output_format = val.string_value();
        else if (key == "suppress_output")
            config.suppress_output = val.bool_value();
        else
            std::cerr << "Warning: unknown config key '" << key << "' (ignored)" << std::endl;
    }

    return config;
}
