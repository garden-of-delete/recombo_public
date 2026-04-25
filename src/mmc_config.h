#pragma once

#include <string>
#include <vector>

struct MmcConfig {
    std::string input_file;
    std::string output_file;
    double z_min = 0.0;
    double z_max = 0.0;
    int q = 1;
    double swap_ratio = 0.8;
    int swap_interval = 5;
    int num_samples = 0;
    int steps_between_samples = 20000;
    long int warmup_steps = 100000;
    int num_chains = 1;
    char sample_mode = 's';
    int seed = 0;               // 0 = use system time
    int min_arc = 0;
    int max_arc = 0;
    int target_recombo_length = 0;
    int block_file_size = 0;
    int time_limit_hours = 0;
    bool suppress_output = false;
    int sequence_type = -1;     // -1 = unset, 0 = inverted, 1 = direct
    int recombo_type = -1;      // -1 = unset, 0-6

    // Returns empty vector if valid, otherwise list of error strings.
    std::vector<std::string> validate() const;

    // Parse recomboParas string (e.g. "is", "ds", "dsp") into sequence_type and recombo_type.
    // Returns true on success, false if the string is unrecognized.
    bool parse_recombo_paras(const std::string& paras);

    // Load config from a JSON file. Keys match stamp_log() output.
    // Returns errors via the errors vector; config is partially populated on failure.
    static MmcConfig from_json(const std::string& filepath, std::vector<std::string>& errors);

    // Parse CLI arguments (positional + flags). Expects argv[0]=program, argv[1]=input, argv[2]=output.
    static MmcConfig from_cli(int argc, const char* const* argv, std::vector<std::string>& errors);

    // Serialize config to a JSON string.
    std::string to_json() const;
};
