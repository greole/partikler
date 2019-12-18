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

enum MSGTYPE {
    INFO,
    INFOSTART,
    INFOEND,
    WARN,
    CRITICAL
};

class MSG {

    private:

        const int message_threshold_;

        const int verbosity_treshold_;

        const MSGTYPE label_;

        const std::string scope_;

        std::ostringstream state_;

    public:

        MSG(int verbosity_treshold, const MSGTYPE label, int message_threshold, std::string scope):
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
            state_(msg.state_.str())
        {};

        std::ostringstream& get_state() {return state_;};

        template<typename T>
        MSG &operator<<(const T &in) {
            state_ << in;
            return *this;
        }

        ~MSG() {
            std::string label_str;

            switch (label_) {
            case INFO:
                label_str = "INFO";
                break;
            case INFOSTART:
                label_str = "INFO START";
                break;
            case INFOEND:
                label_str = "INFO END";
                break;
            case WARN:
                label_str = "WARN";
                break;
            case CRITICAL:
                label_str = "CRITICAL";
                break;
            }

            if (message_threshold_ > verbosity_treshold_) {
                std::cout
                    << "[" << label_str  << scope_ <<  "] "
                    << state_.str()
                    << std::endl;
            }

            if (label_ == CRITICAL)  exit(EXIT_FAILURE);
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
            return MSG(verbosity_treshold_, INFOSTART, 3, scope_);
        };

        MSG info_end() {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::system_clock::now();
            MSG msg(verbosity_treshold_, INFOEND, 3, scope_);
            std::chrono::system_clock::time_point invoked {invoked_stack_.back()};
            invoked_stack_.pop_back();
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - invoked).count();
            msg << time << "ms ";
            return msg;
        };

        MSG info() {
            return MSG(verbosity_treshold_, INFO, 3, scope_);
        };

        MSG verbose(int i) {
            return MSG(verbosity_treshold_, INFO, i, scope_);
        };

        MSG warn() {
            return MSG(verbosity_treshold_, WARN, 4, scope_);
        };

        MSG critical() {
            return MSG(verbosity_treshold_, CRITICAL, 5, scope_);
        };

  // TODO implement an assert statement
        MSG check(bool b) {
          if (b) {return MSG(verbosity_treshold_, INFO, 0, scope_);}
          else {return MSG(verbosity_treshold_, CRITICAL, 5, scope_);}
        };
};

// // Free functions
// template <class T>
// MSG& operator<<(MSG& obj, T& b)
// {
//     auto state = obj.get_state();
//     state << b;
//     return obj;
// }

#endif
