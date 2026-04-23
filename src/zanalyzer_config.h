#pragma once

#include <string>
#include <vector>

struct ZAnalyzerConfig {
    std::string input_file;
    int target_length = 0;
    int tolerance = 0;
    double z_min = 0.0;
    double z_max = 0.0;
    int q = 0;
    int steps_between_samples = 0;
    int warmup = 0;
    int seed = 0;               // 0 = use system time

    std::vector<std::string> validate() const;

    static ZAnalyzerConfig from_json(const std::string& filepath, std::vector<std::string>& errors);
};
