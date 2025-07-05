#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include "example_utils.h"
#include <vector>
#include <fstream>
#include <utility>

using namespace Cantera;
using std::cout;
using std::endl;

std::vector<std::pair<double, double>> kinetics1(ReactorNet& sim, ThermoPhase& gas, int nsteps, double dt)
{
    std::vector<std::pair<double, double>> timeTemp;
    clock_t t0 = clock();

    for (int i = 0; i <= nsteps; i++) {
        double tm = i * dt;
        sim.advance(tm);
        double T = gas.temperature();
        timeTemp.emplace_back(tm, T);
        std::cout << "time = " << tm << " s, T = " << T << " K" << std::endl;
    }

    clock_t t1 = clock();

    double tmm = 1.0 * (t1 - t0) / CLOCKS_PER_SEC;
    cout << " Tfinal = " << gas.temperature() << endl;
    cout << " time = " << tmm << endl;
    cout << " number of residual function evaluations = "
         << sim.integrator().nEvals() << endl;
    cout << " time per evaluation = " << tmm / sim.integrator().nEvals()
         << endl << endl;

    return timeTemp;
}

int main()
{
    try {
        // Create the solution and thermo object
        auto sol = newSolution("~/Codici/ATLAS/database/Chemistry/WD.yaml");
        auto& gas = *sol->thermo();
        gas.setState_TPY(1000.0, OneAtm, "CH4:0.20, O2:0.8");

        // Reactor setup
        IdealGasReactor reactor(sol);
        ReactorNet sim;
        sim.addReactor(reactor);
        sim.setTolerances(1e-7, 1e-7);

        // Integration settings
        int nsteps = 1000;
        double dt = 8.e-3 / nsteps;

        // Run the simulation
        auto timeTemp = kinetics1(sim, gas, nsteps, dt);

        // Write results
        std::ofstream outFile("time_temperature.dat");
        outFile << "# Time (s)\tTemperature (K)\n";
        for (const auto& [time, temp] : timeTemp) {
            outFile << time << "\t" << temp << "\n";
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
