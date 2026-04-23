#include "mmc_config.h"

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
