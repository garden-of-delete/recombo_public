#include "mmc_config.h"
#include "json11.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

bool MmcConfig::parse_recombo_paras(const std::string& paras) {
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
        // length == 3, e.g. "dsp", "dsn", "dsa"
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

std::vector<std::string> MmcConfig::validate() const {
    std::vector<std::string> errors;

    if (input_file.empty())
        errors.push_back("input_file is required");
    if (output_file.empty())
        errors.push_back("output_file is required");

    if (num_samples <= 0 && time_limit_hours <= 0)
        errors.push_back("at least one halting criterion required: num_samples (-n) > 0 or time_limit_hours (-t) > 0");

    std::string valid_modes = "asbfmr";
    if (valid_modes.find(sample_mode) == std::string::npos)
        errors.push_back("sample_mode (-mode) must be one of: a, s, b, f, m, r");

    if (num_chains > 1) {
        if (z_min <= 0)
            errors.push_back("z_min (-zmin) must be > 0 when using multiple chains");
        if (z_max <= 0)
            errors.push_back("z_max (-zmax) must be > 0 when using multiple chains");
        if (z_min >= z_max)
            errors.push_back("z_min (-zmin) must be less than z_max (-zmax) when using multiple chains");
    }

    if (sample_mode == 'f' || sample_mode == 'r') {
        if (sequence_type == -1 || recombo_type == -1)
            errors.push_back("reconnection mode (f/r) requires -recomboParas to be specified");
    }

    if ((min_arc > 0 || max_arc > 0) && min_arc >= max_arc)
        errors.push_back("min_arc (-minarc) must be less than max_arc (-maxarc)");

    if (target_recombo_length > 0 && min_arc <= 3)
        errors.push_back("-targetlength requires -minarc > 3");

    return errors;
}

MmcConfig MmcConfig::from_json(const std::string& filepath, std::vector<std::string>& errors) {
    MmcConfig config;

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

    // Key names match stamp_log() output
    for (auto& kv : json.object_items()) {
        const std::string& key = kv.first;
        const json11::Json& val = kv.second;

        if (key == "input_file")
            config.input_file = val.string_value();
        else if (key == "output_file")
            config.output_file = val.string_value();
        else if (key == "zmin")
            config.z_min = val.number_value();
        else if (key == "zmax")
            config.z_max = val.number_value();
        else if (key == "q")
            config.q = val.int_value();
        else if (key == "target_swap_ratio")
            config.swap_ratio = val.number_value();
        else if (key == "swap_interval")
            config.swap_interval = val.int_value();
        else if (key == "n_samples")
            config.num_samples = val.int_value();
        else if (key == "steps_between_samples")
            config.steps_between_samples = val.int_value();
        else if (key == "warmup")
            config.warmup_steps = val.int_value();
        else if (key == "initial_n_chains")
            config.num_chains = val.int_value();
        else if (key == "sample_mode") {
            std::string s = val.string_value();
            if (!s.empty())
                config.sample_mode = s[0];
        }
        else if (key == "seed")
            config.seed = val.int_value();
        else if (key == "min_arc")
            config.min_arc = val.int_value();
        else if (key == "max_arc")
            config.max_arc = val.int_value();
        else if (key == "target_recombo_length")
            config.target_recombo_length = val.int_value();
        else if (key == "block_file_size")
            config.block_file_size = val.int_value();
        else if (key == "time_limit")
            config.time_limit_hours = val.int_value();
        else if (key == "suppress_output")
            config.suppress_output = val.bool_value();
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
        else if (key == "recombo_paras") {
            if (!config.parse_recombo_paras(val.string_value()))
                errors.push_back("unrecognized recombo_paras: " + val.string_value());
        }
        else {
            std::cerr << "Warning: unknown config key '" << key << "' (ignored)" << std::endl;
        }
    }

    return config;
}

std::string MmcConfig::to_json() const {
    json11::Json json = json11::Json::object {
        {"input_file", input_file},
        {"output_file", output_file},
        {"zmin", z_min},
        {"zmax", z_max},
        {"q", q},
        {"target_swap_ratio", swap_ratio},
        {"swap_interval", swap_interval},
        {"n_samples", num_samples},
        {"steps_between_samples", steps_between_samples},
        {"warmup", (int)warmup_steps},
        {"initial_n_chains", num_chains},
        {"sample_mode", std::string(1, sample_mode)},
        {"seed", seed},
        {"block_file_size", block_file_size},
        {"time_limit", time_limit_hours},
        {"suppress_output", suppress_output},
        {"min_arc", min_arc},
        {"max_arc", max_arc},
        {"target_recombo_length", target_recombo_length},
        {"sequence_type", sequence_type},
        {"recombo_type", recombo_type},
    };
    return json.dump();
}
