#pragma once

#include <stdio.h>
#include <time.h>
#include <iostream> 
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <netdb.h>
#define HOST_NAME_MAX 80

#include <vector> 
#include <regex>

namespace util
{
    std::string time_and_user(); ///< (Brief description)

    void print_time() ; ///< (Brief description)

    void compile_time(const char[11], const char[11], const char[11], const char[11]) ; ///< (Brief description)

    void print(std::vector<std::string> &s); ///< (Brief description)

    void print(std::vector<int> &v); ///< (Brief description)

    void print(std::vector<float> &v); ///< (Brief description)

    std::vector<std::string> strip_white_space(std::vector<std::string> &s) ; ///< (Brief description)

    void print_run_details() ; ///< (Brief description)

    void pp(const std::string &s) ; ///< (Brief description)

    bool has_only_spaces(const std::string &str) ; ///< (Brief description)

    std::vector<int> unique_vector_int(const std::vector<int> & in); ///< (Brief description)

    std::vector<std::string> unique_vector_string(const std::vector<std::string> & in); ///< (Brief description)
}
