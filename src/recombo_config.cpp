#include "recombo_config.h"
#include "json11.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

bool RecomboConfig::parse_recombo_paras(const std::string& paras) {
    if (paras.length() < 2 || paras.length() > 3)
        return false;

    char prefix = paras[0];
    if (prefix == 'i')
        sequence_type = 0;
    else if (prefix == 'd')
        sequence_type = 1;
    else
        return false;

    if (paras.length() == 2) {
        char code = paras[1];
        if (code == 's')      recombo_type = 0;
        else if (code == 'p') recombo_type = 1;
        else if (code == 'n') recombo_type = 2;
        else if (code == 'a') recombo_type = 3;
        else return false;
    } else {
        if (paras[1] != 's')
            return false;
        char code = paras[2];
        if (code == 'p')      recombo_type = 4;
        else if (code == 'n') recombo_type = 5;
        else if (code == 'a') recombo_type = 6;
        else return false;
    }
    return true;
}

std::vector<std::string> RecomboConfig::validate() const {
    std::vector<std::string> errors;

    if (input_file.empty())
        errors.push_back("input_file is required");
    if (output_file.empty())
        errors.push_back("output_file is required");
    if (min_arc <= 0)
        errors.push_back("-minarc must be > 0");
    if (max_arc <= 0)
        errors.push_back("-maxarc must be > 0");
    if (min_arc >= max_arc)
        errors.push_back("-minarc must be less than -maxarc");
    if (sequence_type == -1 || recombo_type == -1)
        errors.push_back("-recomboParas is required");
    if (output_format != "cube" && output_format != "text" && output_format != "newsud")
        errors.push_back("-fmt must be one of: cube, text, newsud");

    return errors;
}

RecomboConfig RecomboConfig::from_cli(int argc, const char* const* argv, std::vector<std::string>& errors) {
    RecomboConfig config;

    if (argc < 3) {
        errors.push_back("usage: recombo input_file output_file [options]");
        return config;
    }

    config.input_file = argv[1];
    config.output_file = argv[2];

    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-config") {
            errors.push_back("-config and CLI options are mutually exclusive");
            return config;
        }
        else if (arg == "-minarc" && i + 1 < argc)  { config.min_arc = atoi(argv[++i]); }
        else if (arg == "-maxarc" && i + 1 < argc)   { config.max_arc = atoi(argv[++i]); }
        else if (arg == "-targetlength" && i + 1 < argc) { config.target_length = atoi(argv[++i]); }
        else if (arg == "-seed" && i + 1 < argc)    { config.seed = atoi(argv[++i]); }
        else if (arg == "-fmt" && i + 1 < argc)     { config.output_format = argv[++i]; }
        else if (arg == "+s")                        { config.suppress_output = true; }
        else if (arg == "-recomboParas" && i + 1 < argc) {
            std::string paras = argv[++i];
            if (!config.parse_recombo_paras(paras))
                errors.push_back("unrecognized -recomboParas value '" + paras + "'");
        }
        else {
            errors.push_back("unrecognized option '" + arg + "'");
        }
    }

    return config;
}

RecomboConfig RecomboConfig::from_json(const std::string& filepath, std::vector<std::string>& errors) {
    RecomboConfig config;

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
        else if (key == "min_arc")
            config.min_arc = val.int_value();
        else if (key == "max_arc")
            config.max_arc = val.int_value();
        else if (key == "target_length")
            config.target_length = val.int_value();
        else if (key == "seed")
            config.seed = val.int_value();
        else if (key == "output_format")
            config.output_format = val.string_value();
        else if (key == "suppress_output")
            config.suppress_output = val.bool_value();
        else if (key == "recombo_paras") {
            if (!config.parse_recombo_paras(val.string_value()))
                errors.push_back("unrecognized recombo_paras: " + val.string_value());
        }
        else if (key == "sequence_type") {
            if (val.is_string()) {
                std::string s = val.string_value();
                if (s == "inverted") config.sequence_type = 0;
                else if (s == "direct") config.sequence_type = 1;
                else errors.push_back("unknown sequence_type: " + s);
            } else {
                config.sequence_type = val.int_value();
            }
        }
        else if (key == "recombo_type") {
            if (val.is_string()) {
                std::string s = val.string_value();
                if (s == "standard")           config.recombo_type = 0;
                else if (s == "positive")      config.recombo_type = 1;
                else if (s == "negative")      config.recombo_type = 2;
                else if (s == "automatic")     config.recombo_type = 3;
                else if (s == "standard_positive") config.recombo_type = 4;
                else if (s == "standard_negative") config.recombo_type = 5;
                else if (s == "standard_automatic") config.recombo_type = 6;
                else errors.push_back("unknown recombo_type: " + s);
            } else {
                config.recombo_type = val.int_value();
            }
        }
        else {
            std::cerr << "Warning: unknown config key '" << key << "' (ignored)" << std::endl;
        }
    }

    return config;
}
