#pragma once

#include <string>
#include <vector>

struct BfacfConfig {
    std::string input_file;
    std::string output_file;
    double z = 0.0;
    int q = 1;
    int steps = 0;
    int save_interval = 0;      // save every k steps; 0 = save only final
    int block_file_size = 0;    // conformations per .b file; 0 = all in one file
    int seed = 0;               // 0 = use system time
    std::string output_format = "cube";  // "cube", "text", "newsud"
    bool suppress_output = false;

    std::vector<std::string> validate() const;

    static BfacfConfig from_cli(int argc, const char* const* argv, std::vector<std::string>& errors);
    static BfacfConfig from_json(const std::string& filepath, std::vector<std::string>& errors);
};
