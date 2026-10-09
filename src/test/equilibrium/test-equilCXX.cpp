// test-equilCXX: the Cantera references of the equilibrium sweeps of test-CEA, written to
// reference/<mechanism>.dat (run from test/equilibrium; plotted against <mechanism>/FLINT-CEA.txt by
// eq-verification.py). The same problem as CEA_solve: equilibrium at the internal energy and the density
// of the initial state (UV), from the temperature, density and mass fractions of test-CEA:
//   WD, ZK, TSR-GP-24  T = 1000 K, rho = 3.25 kg/m3, Y_O2 = of/(of+1), Y_CH4 = 1/(of+1),
//                      1000 mixture ratios of from 0.01 to 100 (log spaced)
//   Ecker              T = 3000 K, Y = H2O 0.5, CL2 0.2, H2 0.2, O2 0.1,
//                      1000 pressures from 1e-5 to 100 bar (log spaced), rho = p/(R T)
// Columns: the swept variable (of, or p [Pa]) and the equilibrium temperature [K].
//   test-equilCXX [mechanism...]   every sweep, or the ones named

#include "cantera/core.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace Cantera;

struct Sweep {
    std::string name, yaml;
    bool pressure;   // false: mixture ratio of O2/CH4 at fixed density; true: pressure at fixed composition
};

int main(int argc, char** argv)
{
    const int N = 1000;
    const std::vector<Sweep> sweeps = {
        {"WD", "WD.yaml", false},
        {"ZK", "ZK.yaml", false},
        {"TSR-GP-24", "TSR-GP-24.yaml", false},
        {"Ecker", "ecker.yaml", true},
    };
    std::vector<std::string> only(argv + 1, argv + argc);

    try {
        std::filesystem::create_directories("reference");
        for (const auto& sw : sweeps) {
            if (!only.empty() && std::find(only.begin(), only.end(), sw.name) == only.end()) {
                continue;
            }
            auto sol = newSolution("../../database/" + sw.name + "/" + sw.yaml);
            auto gas = sol->thermo();
            std::string file = "reference/" + sw.name + ".dat";
            std::ofstream out(file);
            out << "# test-CEA reference: " << sw.name << ", test-equilCXX, Cantera " << CANTERA_VERSION << "\n"
                << "# database/" << sw.name << "/" << sw.yaml << ", equilibrium at constant U and V (UV)\n"
                << (sw.pressure ? "# p [Pa]  T_eq [K]\n" : "# of = Y_O2/Y_CH4  T_eq [K]\n")
                << std::setprecision(17);
            int failed = 0;
            for (int i = 0; i < N; i++) {
                double x;
                if (sw.pressure) {
                    x = 1e5 * 1e-5 * std::pow(10.0, i * std::log10(100.0 / 1e-5) / (N - 1));
                    gas->setState_TPY(3000.0, x, "H2O:0.5, CL2:0.2, H2:0.2, O2:0.1");
                } else {
                    x = 0.01 * std::pow(100.0 / 0.01, double(i) / (N - 1));
                    std::vector<double> Y(gas->nSpecies(), 0.0);
                    Y[gas->speciesIndex("O2")] = x / (x + 1);
                    Y[gas->speciesIndex("CH4")] = 1 / (x + 1);
                    gas->setMassFractions(Y.data());
                    gas->setState_TD(1000.0, 3.25);
                }
                try {
                    gas->equilibrate("UV");
                    out << x << "\t" << gas->temperature() << "\n";
                } catch (CanteraError& err) {
                    failed++;
                }
            }
            std::cout << sw.name << ": wrote " << file << " (" << N - failed << " states, " << failed
                      << " without equilibrium)" << std::endl;
        }
        appdelete();
        return 0;
    } catch (std::exception& err) {
        std::cout << err.what() << std::endl;
        appdelete();
        return -1;
    }
}
