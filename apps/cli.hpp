// Minimal "--key value" argument parser shared by the command-line tools.
#pragma once

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

namespace cli {

class Args {
 public:
  Args(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      std::string key = argv[i];
      if (key.rfind("--", 0) != 0) fail("unexpected argument: " + key);
      key = key.substr(2);
      if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
        values_[key] = argv[++i];
      else
        values_[key] = "";
    }
  }

  bool has(const std::string& key) const { return values_.count(key) > 0; }

  std::string str(const std::string& key, const std::string& def) const {
    auto it = values_.find(key);
    return it == values_.end() ? def : it->second;
  }

  double num(const std::string& key, double def) const {
    auto it = values_.find(key);
    if (it == values_.end()) return def;
    char* end = nullptr;
    const double v = std::strtod(it->second.c_str(), &end);
    if (end == it->second.c_str() || *end != '\0') fail("invalid number for --" + key);
    return v;
  }

  [[noreturn]] static void fail(const std::string& msg) {
    std::cerr << "error: " << msg << "\n";
    std::exit(2);
  }

 private:
  std::map<std::string, std::string> values_;
};

}  // namespace cli
