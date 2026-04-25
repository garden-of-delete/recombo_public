#pragma once

#include "json11.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

inline json11::Json load_test_data(const std::string& filename) {
    std::string path = std::string(TEST_DATA_DIR) + "/" + filename;
    std::ifstream file(path.c_str());
    if (!file.is_open())
        throw std::runtime_error("cannot open test data file: " + path);
    std::stringstream buf;
    buf << file.rdbuf();
    std::string err;
    json11::Json json = json11::Json::parse(buf.str(), err);
    if (!err.empty())
        throw std::runtime_error("JSON parse error in " + path + ": " + err);
    return json;
}
