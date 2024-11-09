#ifndef ResultsSaver_hpp
#define ResultsSaver_hpp
/*
Copyright 2019, Brown University, Providence, RI.
                        All Rights Reserved
Permission to use, copy, modify, and distribute this software and
its documentation for any purpose other than its incorporation into a
commercial product or service is hereby granted without fee, provided
that the above copyright notice appear in all copies and that both
that copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.
BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
// Cole Foster
// 2022-06-06

#include <string.h>
#include <sys/stat.h>
#include <vector>
#include <chrono>
#include <fstream>

/**
 * @brief Class to write results to a .txt file. Used for optimization
 *
 */
class ResultsSaver {
  public:
    ResultsSaver(){};
    ResultsSaver(std::string filename) : _filename(filename) {
        printf("Initiated File Saver for: %s\n", _filename.c_str());
        _getDateTime();  // retrieve datetime

        // looping through filenames
        if (_fileExists(_filename)) {
            printf("Whoops! filename: '%s' exists, appending date and time to end...\n", _filename.c_str());
            _filename = std::string(filename).append(_dateTime);
            printf("New filename: %s\n", _filename.c_str());
        }

        // initialize
        _file.open(_filename);
        _file.close();
    };
    ResultsSaver(ResultsSaver const& ref) {
        _filename = ref._filename;
        _dateTime = ref._dateTime;
        _del = ref._del;
    };

    ~ResultsSaver(){};

    // member vars
    std::string _filename = "";
    std::string _dateTime = "";
    std::ofstream _file;
    char _del = ',';

    inline void _getDateTime() {
        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string s(30, '\0');
        std::strftime(&s[0], s.size(), "_%Y-%m-%d_%H:%M:%S", std::localtime(&now));
        _dateTime = s;
        return;
    }

    inline bool _fileExists(const std::string& filename) {
        struct stat buf;
        if (stat(filename.c_str(), &buf) != -1) {
            return true;
        }
        return false;
    }

    inline void printLine(std::string const line) {
        _file.open(_filename.c_str(), std::ios_base::app);  // append to file, instead of overwrite
        _file << line.c_str() << '\n';
        _file.close();
    }

    inline void printDateTime() { printLine(_dateTime); }

    // saves vector, as an entire new line, '\n' after
    template <class T>
    inline void saveVectorLine(std::vector<T> vec) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        for (int i = 0; i < (int)vec.size(); i++) _file << vec[i] << _del;
        _file << '\n';
        _file.close();
    }

    // saves vector using delimiter of ';' and wrapped by { }
    // ASSUMES FILE IS ALREADY OPEN
    template <class T>
    inline void saveVectorElement(std::vector<T> vec) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << "{;";
        for (int i = 0; i < (int)vec.size(); i++) _file << vec[i] << ';';
        _file << "}" << _del;
        _file.close();
    }

    template <class T>
    inline void saveElement(T element) {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << element << _del;
        _file.close();
    }

    inline void newLine() {
        _file.open(_filename, std::ios_base::app);  // append to file, instead of overwrite
        _file << '\n';
        _file.close();
    }
};

#endif  // ResultsSaver_hpp