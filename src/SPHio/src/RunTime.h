#ifndef RUNTIME_H
#define RUNTIME_H

#include <chrono>
#include "git.h"

class RunTime {

    class MSG {

        private:

            const int message_threshold_;

            const int verbosity_treshold_;

            const char* label_;

            std::ostringstream state;

        public:

            MSG(int verbosity_treshold, char* label, int message_threshold):
                message_threshold_(message_threshold),
                verbosity_treshold_(verbosity_treshold),
                label_(label)
            {
            };


            //void register_field(std::vector<float> f);

            MSG(const RunTime::MSG& msg):
                message_threshold_(msg.message_threshold_),
                verbosity_treshold_(msg.verbosity_treshold_),
                label_(msg.label_),
                state(msg.state.str())
            {};

            template<typename T>
            MSG &operator<<(T &in) {
                state << in;
                return *this;
            };

            ~MSG() {
                if (message_threshold_ > verbosity_treshold_ ) {
                    std::cout  << "[" << label_ << "] " << state.str() << std::endl;
                }
            };
    };

    private:

        // 0 - vvverbose
        // 1 - vverbose
        // 2 - verbose
        // 3 - info
        // 4 - warning
        // 5 - critical
        int verbosity_treshold_;

        std::chrono::system_clock::time_point invoked;

    public:

        RunTime(int verbosity_treshold):
            verbosity_treshold_(verbosity_treshold)
        {
            std::cout << "[INFO] Starting SPH RunTime " << std::endl;
            print_git_state();
        };

        void print_git_state() {
            if(GIT_RETRIEVED_STATE) {
                std::cout
                    << "[INFO] Git commit "
                    << GIT_HEAD_SHA1
                    << std::endl;
                if(GIT_IS_DIRTY) {
                    std::cerr
                        << "[WARN] There were uncommitted changes."
                        << std::endl;
                }
            }
            else {
                std::cerr
                    << "[WARN] Failed to get the current git state."
                    << "Is this a git repo?"
                    << std::endl;
            }
        };

        MSG info_begin() {
            invoked = std::chrono::system_clock::system_clock::now();
            return MSG(verbosity_treshold_, "INFO START", 3);
        };

        MSG info_end() {
            std::chrono::system_clock::time_point end = std::chrono::system_clock::system_clock::now();
            MSG msg(verbosity_treshold_, "INFO END", 3);
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - invoked).count();
            msg << time << "ms ";
            return msg;
        };


        MSG info() {
            return MSG(verbosity_treshold_, "INFO", 3);
        };

        MSG verbose(int i) {
            return MSG(verbosity_treshold_, "INFO", i);
        };

        MSG warn(int i) {
            return MSG(verbosity_treshold_, "WARN", 4);
        };

        MSG critical(int i) {
            return MSG(verbosity_treshold_, "CRITICAL", 5);
        };
};

#endif
