#ifndef chaospp_io_h
#define chaospp_io_h

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <stdarg.h>  // for va_start, etc
#include <memory>    // for unique_ptr
#include <iterator>  // for istream_iterator
#include <limits>    // for numeric_limits


//! splits a list of strings by the delimiter
//! e.g. split("2,3,4", ",") returns {2,3,4}
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delimiter)) {
        elems.push_back(item);
    }
    return elems;
}

//! joins a list of strings in a string separated by s
//! e.g. join({2,3,4}, ",") returns "2,3,4"
template <typename T>
std::string join(std::vector<T> const& strings, std::string const& s) {
    std::stringstream result;
    for (unsigned int i = 0; i < strings.size(); i ++) {
        result << strings[i];
        if (i != strings.size() - 1)
            result << s;
    }
    return result.str();
}


// helper to format strings, see http://stackoverflow.com/a/8098080/931303
std::string format(const std::string fmt_str, ...) {
    int final_n, n = ((int)fmt_str.size()) * 2; /* reserve 2 times as much as the length of the fmt_str */
    std::string str;
    std::unique_ptr<char[]> formatted;
    va_list ap;
    while(1) {
        formatted.reset(new char[n]); /* wrap the plain char array into the unique_ptr */
        strcpy(&formatted[0], fmt_str.c_str());
        va_start(ap, fmt_str);
        final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
        va_end(ap);
        if (final_n < 0 || final_n >= n)
            n += abs(final_n - n + 1);
        else
            break;
    }
    return std::string(formatted.get());
}


namespace io {

    // Exports a file of arbitrary columns
    template <typename T>
    void save(std::vector<std::vector<T> > const& data, std::string file_name) {
        std::ofstream file;
        file.open(file_name.c_str());
        if(!file.is_open())
        {
            std::cout << "file \"" << file_name << "\" not found" << std::endl;
            exit(1);
        }

        file.precision(std::numeric_limits<T>::digits10 + 1);
        for(unsigned int i = 0; i < data.size(); i++)
        {
            for(unsigned int j = 0; j < data[i].size(); j++)
            {
                file << data[i][j];
                if (j != data[i].size() - 1)
                    file << " ";
            }
            file << std::endl;
        }
        file.close();
    }

    // Reads a file of arbitrary columns
    std::vector<std::vector<double> > load(std::string file_name) {
        std::vector<std::vector<double> > data;

        std::ifstream file;
        file.open(file_name.c_str());
        if(!file.is_open())
        {
            std::cout << "file \"" << file_name << "\" not found" << std::endl;
            exit(1);
        }

        file.precision(std::numeric_limits<double>::digits10 + 1);

        while(true)
        {
            if(!file)
                break;

            std::string line;
            getline(file, line);
            std::istringstream iss(line);

            std::vector<double> row;
            copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(),
                 back_inserter(row));

            if (row.size())
                data.push_back(row);
        }
        file.close();

        return data;
    }
}

#endif
