#pragma once

#include <string>
#include <vector>

struct RecomboConfig {
    std::string input_file;
    std::string output_file;
    int min_arc = 0;
    int max_arc = 0;
    int target_length = 0;      // 0 = no length constraint
    int sequence_type = -1;     // -1 = unset, 0 = inverted, 1 = direct
    int recombo_type = -1;      // -1 = unset, 0-6
    int seed = 0;               // 0 = use system time
    std::string output_format = "cube";  // "cube", "text", "newsud"
    bool suppress_output = false;

    std::vector<std::string> validate() const;
    bool parse_recombo_paras(const std::string& paras);

    static RecomboConfig from_cli(int argc, const char* const* argv, std::vector<std::string>& errors);
    static RecomboConfig from_json(const std::string& filepath, std::vector<std::string>& errors);
};
