// Reports what the OpenMP module sees, and checks that a parallel region really
// runs on the requested number of threads. The pytest beside this file parses
// the key = value lines.

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "omp/omp_environment.h"

namespace {

/// take_threads_option on a copy of `arguments`; what is left of argv, joined.
bool take(std::vector<std::string> arguments, int* threads, std::string* rest) {
    std::vector<char*> argv;
    for (std::string& a : arguments) {
        argv.push_back(a.data());
    }
    int argc = static_cast<int>(argv.size());
    const bool valid = cglbm::omp::take_threads_option(&argc, argv.data(), threads);
    rest->clear();
    for (int n = 1; n < argc; ++n) {
        *rest += (n > 1 ? " " : "") + std::string(argv[n]);
    }
    return valid;
}

/// The --threads option the programs that read their own arguments take.
void report_take_threads_option() {
    int threads = 0;
    std::string rest;
    const bool valid = take({"program", "E8", "--threads=3", "1e4"}, &threads, &rest);
    std::cout << "take_threads_valid = " << (valid ? 1 : 0) << "\n";
    std::cout << "take_threads_value = " << threads << "\n";
    std::cout << "take_threads_rest = " << rest << "\n";
    threads = 0;
    const bool absent = take({"program", "E8"}, &threads, &rest);
    std::cout << "take_threads_absent = " << ((absent && threads == 0 && rest == "E8") ? 1 : 0)
              << "\n";
    int rejected = 0;
    for (const char* bad :
         {"--threads=0", "--threads=-2", "--threads=x", "--threads=", "--threads=2.5"}) {
        threads = 7;
        if (!take({"program", bad}, &threads, &rest) && threads == 7 && rest.empty()) {
            ++rejected;
        }
    }
    std::cout << "take_threads_rejected = " << rejected << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    const int requested = (argc > 1) ? std::atoi(argv[1]) : 0;

    const int in_force = cglbm::omp::set_thread_count(requested);

    std::cout << "available = " << (cglbm::omp::available() ? 1 : 0) << "\n";
    std::cout << "requested = " << requested << "\n";
    std::cout << "max_threads = " << in_force << "\n";
    std::cout << "describe = " << cglbm::omp::describe() << "\n";

    // Count the threads that actually enter a parallel region, and check that
    // every thread id appears exactly once.
    int observed = 0;
    int id_sum = 0;
#pragma omp parallel reduction(+ : observed, id_sum)
    {
        observed += 1;
        id_sum += cglbm::omp::thread_id();
    }

    std::cout << "observed_threads = " << observed << "\n";
    // 0 + 1 + ... + (n-1)
    std::cout << "expected_id_sum = " << observed * (observed - 1) / 2 << "\n";
    std::cout << "id_sum = " << id_sum << "\n";

    const double start = cglbm::omp::wall_time();
    const double elapsed = cglbm::omp::wall_time() - start;
    std::cout << "wall_time_monotonic = " << (elapsed >= 0.0 ? 1 : 0) << "\n";

    report_take_threads_option();

    return 0;
}
