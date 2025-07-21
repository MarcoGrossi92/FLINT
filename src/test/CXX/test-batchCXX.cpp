#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include "example_utils.h"
#include <vector>
#include <fstream>
#include <utility>

using namespace Cantera;
using std::cout;
using std::endl;
#include <chrono>

std::vector<std::pair<double, double>> simulateReactorTemperature(ReactorNet& sim, ThermoPhase& gas, int nsteps, double dt)
{
    std::vector<std::pair<double, double>> timeTemp;

    for (int i = 0; i <= nsteps; i++) {
        double tm = i * dt;
        sim.advance(tm);
        double T = gas.temperature();
        timeTemp.emplace_back(tm, T);
    }

    return timeTemp;
}

int main()
{
    try {
        // Create the solution and thermo object
        auto sol1 = newSolution("WD/WD.yaml");
        auto& gas = *sol1->thermo();
        gas.setState_TPY(1000.0, 100000.0, "CH4:0.20, O2:0.8");

        // Reactor setup
        IdealGasReactor reactor(sol1);
        ReactorNet sim;
        sim.addReactor(reactor);
        sim.setTolerances(1e-7, 1e-7);

        // Integration settings
        int nsteps = 1000;
        double dt = 8.e-3 / nsteps;

        // Run the simulation
        auto t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp1 = simulateReactorTemperature(sim, gas, nsteps, dt);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t1 - t0;
        std::cout << "* Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;

        // Write results for first simulation
        std::ofstream outFile("WD-batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp1) {
            outFile << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Troyes
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol2 = newSolution("Troyes/troyes.yaml");
        auto& gas2 = *sol2->thermo();
        gas2.setState_TPY(1000.0, 100000.0, "H2:0.00534, CL2:0.18796, N2:0.80670");

        // Reactor setup
        IdealGasReactor reactor2(sol2);
        ReactorNet sim2;
        sim2.addReactor(reactor2);
        sim2.setTolerances(1e-7, 1e-7);

        // Integration settings
        int nsteps2 = 1000;
        double dt2 = 0.02 / nsteps2;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp2 = simulateReactorTemperature(sim2, gas2, nsteps2, dt2);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "* Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;

        // Write results for second simulation
        std::ofstream outFile2("Troyes-batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp2) {
            outFile2 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Ecker
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol3 = newSolution("Ecker/ecker.yaml");
        auto& gas3 = *sol3->thermo();
        gas3.setState_TPY(1000.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        IdealGasReactor reactor3(sol3);
        ReactorNet sim3;
        sim3.addReactor(reactor3);
        sim3.setTolerances(1e-7, 1e-7);

        // Integration settings
        int nsteps3 = 1000;
        double dt3 = 0.02 / nsteps3;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp3 = simulateReactorTemperature(sim3, gas3, nsteps3, dt3);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "* Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;

        // Write results for second simulation
        std::ofstream outFile3("Ecker-batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp3) {
            outFile3 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Cross
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol4 = newSolution("Cross/cross.yaml");
        auto& gas4 = *sol4->thermo();
        gas4.setState_TPY(1005.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        IdealGasReactor reactor4(sol4);
        ReactorNet sim4;
        sim4.addReactor(reactor4);
        sim4.setTolerances(1e-7, 1e-7);
        sim4.setMaxSteps(100000);

        // Integration settings
        int nsteps4 = 1000;
        double dt4 = 0.02 / nsteps4;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        auto timeTemp4 = simulateReactorTemperature(sim4, gas4, nsteps4, dt4);
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "* Cantera-CXX time = "
          << std::scientific << elapsed.count()
          << std::endl;

        // Write results for second simulation
        std::ofstream outFile4("Cross-batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp4) {
            outFile4 << time << "\t" << temp << "\n";
        }

        appdelete();
        return 0;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}
