// test-batchCXX: constant-volume batch reactor with Cantera (IdealGasReactor) for the cases of
// test/batch/cases.txt, the counterpart of test-batchF. Run from test/batch:
//   test-batchCXX                        asks for the mode, then runs every case
//   test-batchCXX verification           every case, 1000 steps: <case>/batch-CXX.dat
//   test-batchCXX performance            every case, 1 step: the times in comp-batch-cantera.dat
//   test-batchCXX --reference [case...]  the references of test-batchF check, reference/<case>.dat
//                                        (every case, or the cases named), 1000 steps
// Verification and references at rtol 1e-10, atol 1e-15: with an absolute tolerance of 1e-7 on the mass
// fractions CVODE does not resolve the radicals that start at 0 (Gerlinger did not ignite in 2e-4 s).
// Cantera reads element-standard-entropies.yaml (the standard entropies of the elements) from the working
// directory, test/batch.

#include "cantera/zerodim.h"
#include "cantera/base/global.h"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace Cantera;

struct BatchCase {
    std::string name, yaml, composition;
    double tend = 0, p = 0, T = 0;
};

// The cases of cases.txt (format in its header); the species:value tokens become one Cantera
// composition string.
std::vector<BatchCase> readCases(const std::string& file)
{
    std::ifstream in(file);
    if (!in) {
        throw std::runtime_error("cannot open " + file + " (run from test/batch)");
    }
    std::vector<BatchCase> cases;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream tokens(line);
        BatchCase c;
        std::string general, y;
        if (!(tokens >> c.name) || c.name[0] == '#') {
            continue;
        }
        if (!(tokens >> c.yaml >> general >> c.tend >> c.p >> c.T)) {
            throw std::runtime_error(file + ": a case is: case yaml general tend p T species:Y ...: " + line);
        }
        while (tokens >> y) {
            c.composition += (c.composition.empty() ? "" : ", ") + y;
        }
        cases.push_back(c);
    }
    return cases;
}

// Integrates one case and returns (time, T) at nsteps equal steps up to tend, and the wall time.
std::vector<std::pair<double, double>> runCase(const BatchCase& c, int nsteps, double rtol, double atol,
                                               double& elapsed)
{
    auto sol = newSolution("../../database/" + c.name + "/" + c.yaml);
    sol->thermo()->setState_TPY(c.T, c.p, c.composition);
    auto reactor = newReactorBase("IdealGasReactor", sol);
    ReactorNet sim(reactor);
    sim.setTolerances(rtol, atol);
    sim.setMaxSteps(1000000);

    std::vector<std::pair<double, double>> timeTemp;
    timeTemp.reserve(nsteps);
    double dt = c.tend / nsteps;
    double time = 0.0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < nsteps; i++) {
        time += dt;
        sim.advance(time);
        timeTemp.emplace_back(time, reactor->temperature());
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration<double>(t1 - t0).count();
    return timeTemp;
}

int main(int argc, char** argv)
{
    const int nstepsVerification = 1000;
    const double rtolTight = 1e-10, atolTight = 1e-15;   // verification and references
    const double rtolPerf = 1e-7, atolPerf = 1e-7;       // performance: the RT = AT of test-batchF

    try {
        std::vector<std::string> args(argv + 1, argv + argc);
        std::string mode = args.empty() ? "" : args[0];
        if (mode.empty()) {
            int sim_type = 0;
            std::cout << "What kind of simulation do you want to run?\n";
            std::cout << "1) verification\n";
            std::cout << "2) performance\n";
            std::cin >> sim_type;
            if (sim_type == 1) {
                mode = "verification";
            } else if (sim_type == 2) {
                mode = "performance";
            } else {
                std::cerr << "Choose 1 or 2!\n";
                return EXIT_FAILURE;
            }
        }
        if (mode != "verification" && mode != "performance" && mode != "--reference") {
            std::cerr << "usage: test-batchCXX [verification | performance | --reference [case...]]\n";
            return 2;
        }
        bool reference = (mode == "--reference");
        bool performance = (mode == "performance");
        int nsteps = performance ? 1 : nstepsVerification;
        double rtol = performance ? rtolPerf : rtolTight;
        double atol = performance ? atolPerf : atolTight;
        std::vector<std::string> only(args.begin() + (args.empty() ? 0 : 1), args.end());

        std::cout << "Running with nstep = " << nsteps << "\n";

        std::vector<BatchCase> cases = readCases("cases.txt");
        std::vector<std::pair<std::string, double>> summaryData;
        for (const auto& name : only) {
            bool found = false;
            for (const auto& c : cases) {
                found = found || c.name == name;
            }
            if (!found) {
                throw std::runtime_error("case " + name + " is not in cases.txt");
            }
        }

        for (const auto& c : cases) {
            if (!only.empty() && std::find(only.begin(), only.end(), c.name) == only.end()) {
                continue;
            }
            double elapsed = 0;
            auto timeTemp = runCase(c, nsteps, rtol, atol, elapsed);
            std::cout << c.name << " Cantera-CXX time = " << std::scientific << elapsed << std::endl;
            summaryData.emplace_back(c.name, elapsed);

            std::string file;
            if (reference) {
                std::filesystem::create_directories("reference");
                file = "reference/" + c.name + ".dat";
            } else {
                std::filesystem::create_directories(c.name);
                file = c.name + "/batch-CXX.dat";
            }
            std::ofstream out(file);
            if (reference) {
                out << "# test-batchF check reference: " << c.name << ", test-batchCXX --reference, Cantera "
                    << CANTERA_VERSION << "\n"
                    << "# database/" << c.name << "/" << c.yaml << ", IdealGasReactor, rtol " << rtol
                    << ", atol " << atol << "\n"
                    << "# T0 = " << c.T << " K, p0 = " << c.p << " Pa, Y = " << c.composition << "\n"
                    << "# time [s]  T [K]\n";
            }
            out << std::setprecision(17);
            for (const auto& [time, temp] : timeTemp) {
                out << time << "\t" << temp << "\n";
            }
            if (reference) {
                std::cout << "  wrote " << file << std::endl;
            }
        }

        if (!reference) {
            std::ofstream summaryFile("comp-batch-cantera.dat");
            for (const auto& [name, time] : summaryData) {
                summaryFile << name << "\t" << std::scientific << time << "\n";
            }
        }

        appdelete();
        return 0;

    } catch (std::exception& err) {
        std::cout << err.what() << std::endl;
        appdelete();
        return -1;
    }
}
