/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <iostream>
#include <vector>
#include <string>       // std::string
#include <sstream>

class MSG {

    private:

        const int message_threshold_;

        const int verbosity_treshold_;

        const char* label_;

        const std::string scope_;

        std::ostringstream state;

    public:

        MSG(int verbosity_treshold, const char* label, int message_threshold, std::string scope):
            message_threshold_(message_threshold),
            verbosity_treshold_(verbosity_treshold),
            label_(label),
            scope_((scope=="")? scope: std::string(" ") + scope)
        {};

        MSG(const MSG& msg):
            message_threshold_(msg.message_threshold_),
            verbosity_treshold_(msg.verbosity_treshold_),
            label_(msg.label_),
            scope_(msg.scope_),
            state(msg.state.str())
        {};

        template<typename T>
        MSG &operator<<(T &in) {
            state << in;
            return *this;
        }

        // template<int>
        // MSG &operator<<(int &in) {
        //   state << "int";
        //   return *this;
        // }


        ~MSG() {
            if (message_threshold_ > verbosity_treshold_ ) {
                std::cout
                    << "[" << label_  << scope_ <<  "] "
                    << state.str()
                    << std::endl;
            }

            if (label_ == "CRITICAL")  exit(EXIT_FAILURE);
        };
};

class Logger {
    // NOTE Verbosity levels
    // 0 - vvverbose
    // 1 - vverbose
    // 2 - verbose
    // 3 - info
    // 4 - warning
    // 5 - critical

    private:

        int verbosity_treshold_;

        std::string scope_;

        std::vector<std::chrono::system_clock::time_point> invoked_stack_;

    public:

        Logger(
                int verbosity_treshold=3,
                std::string scope=""):
            verbosity_treshold_(verbosity_treshold),
            scope_(scope)
            {}

        void set_scope(std::string scope) {scope_ = scope;};

        MSG info_begin() {
            invoked_stack_.push_back(std::chrono::system_clock::system_clock::now());
            return MSG(verbosity_treshold_, "INFO START", 3, scope_);
        };

        MSG info_end() {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::system_clock::now();
            MSG msg(verbosity_treshold_, "INFO END", 3, scope_);
            std::chrono::system_clock::time_point invoked {invoked_stack_.back()};
            invoked_stack_.pop_back();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - invoked).count();
            msg << time << "ms ";
            return msg;
        };

        MSG info() {
            return MSG(verbosity_treshold_, "INFO", 3, scope_);
        };

        MSG verbose(int i) {
            return MSG(verbosity_treshold_, "INFO", i, scope_);
        };

        MSG warn() {
            return MSG(verbosity_treshold_, "WARN", 4, scope_);
        };

        MSG critical() {
            return MSG(verbosity_treshold_, "CRITICAL", 5, scope_);
        };

  // TODO implement an assert statement
        MSG check(bool b) {
          if (b) {return MSG(verbosity_treshold_, "INFO", 0, scope_);}
          else {return MSG(verbosity_treshold_, "CRITICAL", 5, scope_);}
        };
};

#endif
