#include "cantera/zerodim.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/core.h"
#include <iostream>

extern "C" {

void run_batch_reactor(double* final_time, double* T_out, double* Y_out, int* nsp)
{
    using namespace Cantera;

    auto sol = newSolution("h2o2.yaml");
    auto gas = sol->thermo();

    gas->setState_TPX(1000.0, 101325.0, "CH4:1, O2:2, N2:7.52");

    size_t nSpecies = gas->nSpecies();

    auto reactor = new IdealGasReactor(gas);
    // ReactorNet net;
    // net.addReactor(*reactor);

    // double t = 0.0;
    // double tf = *final_time;
    // while (t < tf) {
    //     t = net.step();
    // }

    // *T_out = gas->temperature();
    // auto Y = gas->massFractions();
    // for (size_t i = 0; i < nSpecies; ++i) {
    //     Y_out[i] = Y[i];
    // }
    // *nsp = static_cast<int>(nSpecies);

    //delete reactor;
    //delete gas;
}

}
