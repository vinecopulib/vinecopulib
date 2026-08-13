// Copyright © 2016-2025 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include "rscript.hpp"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <random>
#include <stdexcept>
#include <string>

//! Skips the enclosing test when the build has no Rscript.
#define SKIP_WITHOUT_RSCRIPT()                                                 \
  do {                                                                         \
    if (!VINECOPULIB_HAS_RSCRIPT) {                                            \
      GTEST_SKIP() << "Rscript was not found at configure time";               \
    }                                                                          \
  } while (false)

namespace test_r_parity {

//! @brief A directory for one test's exchange files with R, removed on scope
//! exit.
//!
//! The R scripts write their output under `VINECOPULIB_TEST_OUTDIR`. Giving
//! every test its own directory is what makes concurrent test processes safe:
//! with a shared name they overwrite and delete each other's files.
class RTempDir
{
public:
  RTempDir()
  {
    // The random component separates concurrent test processes; the counter
    // separates fixtures within one process.
    static const auto salt = std::to_string(std::random_device{}());
    static std::atomic<unsigned> counter{ 0 };
    path_ = std::filesystem::temp_directory_path() /
            ("vinecopulib-test-" + salt + "-" + std::to_string(counter++));
    std::filesystem::create_directories(path_);
    set_env("VINECOPULIB_TEST_OUTDIR", path_.string());
  }

  ~RTempDir()
  {
    std::error_code ec; // a failed cleanup must not throw from a destructor
    std::filesystem::remove_all(path_, ec);
  }

  RTempDir(const RTempDir&) = delete;
  RTempDir& operator=(const RTempDir&) = delete;

  //! Path of a file the R script wrote (or that the test writes for itself).
  std::string file(const std::string& name) const
  {
    return (path_ / name).string();
  }

  //! Runs `Rscript <script> <args>`, throwing if it exits non-zero.
  void run(const std::string& script, const std::string& args = "") const
  {
    std::string cmd = std::string(RSCRIPT) + script + " " + args;
#ifdef _WIN32
    // system() goes through `cmd /c`, which keeps quote characters only when
    // there are exactly two of them. The interpreter path and the script path
    // make four, so cmd would strip the outermost pair and mangle the rest;
    // an extra enclosing pair is what it strips instead.
    cmd = "\"" + cmd + "\"";
#endif
    if (system(cmd.c_str()) != 0) {
      throw std::runtime_error("Rscript failed: " + cmd);
    }
  }

private:
  static void set_env(const std::string& name, const std::string& value)
  {
#ifdef _WIN32
    _putenv_s(name.c_str(), value.c_str());
#else
    setenv(name.c_str(), value.c_str(), 1);
#endif
  }

  std::filesystem::path path_;
};
}
