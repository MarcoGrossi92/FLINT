#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"
#include "utils.h"
#include <vector>
#include <fstream>
#include <utility>
#include "cantera/base/AnyMap.h"
#include <string>
#include "cantera/base/global.h"


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

        // AnyMap stats = sim.solverStats();
        // std::cout << "\n=== CVODE Solver Stats ===\n";
        // std::cout << "Steps taken:                 " << stats["steps"].asInt() << "\n";
        // std::cout << "RHS evaluations:             " << stats["rhs_evals"].asInt() << "\n";
        // std::cout << "Nonlinear iterations:        " << stats["nonlinear_iters"].asInt() << "\n";
        // std::cout << "Nonlinear convergence fails: " << stats["nonlinear_conv_fails"].asInt() << "\n";
        // std::cout << "Error test failures:         " << stats["err_test_fails"].asInt() << "\n";
        // std::cout << "Last method order used:      " << stats["last_order"].asInt() << "\n";
        // std::cout << "Stability limit reductions:  " << stats["stab_order_reductions"].asInt() << "\n";

        // std::cout << "Jacobian evaluations:        " << stats["jac_evals"].asInt() << "\n";
        // std::cout << "Linear solver setups:        " << stats["lin_solve_setups"].asInt() << "\n";
        // std::cout << "Linear RHS evaluations:      " << stats["lin_rhs_evals"].asInt() << "\n";
        // std::cout << "Linear iterations:           " << stats["lin_iters"].asInt() << "\n";
        // std::cout << "Linear convergence failures: " << stats["lin_conv_fails"].asInt() << "\n";

        // std::cout << "Preconditioner evaluations:  " << stats["prec_evals"].asInt() << "\n";
        // std::cout << "Preconditioner solves:       " << stats["prec_solves"].asInt() << "\n";
        // std::cout << "Jt-vector setup evals:       " << stats["jt_vec_setup_evals"].asInt() << "\n";
        // std::cout << "Jt-vector product evals:     " << stats["jt_vec_prod_evals"].asInt() << "\n";

        // if (stats.hasKey("step_solve_fails")) {
        //     std::cout << "Step solve failures:         " << stats["step_solve_fails"].asInt() << "\n";
        // }
        
        double T = gas.temperature();
        timeTemp.emplace_back(tm, T);
    }

    return timeTemp;
}

int main()
{
    try {

        std::vector<std::pair<std::string, double>> summaryData;

        std::chrono::high_resolution_clock::time_point t0;
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::duration<double> elapsed = t1 - t0;

        int nsteps = 0;
        double dt = 0;
        int sim_type = 0;

        std::cout << "What kind of simulation do you want to run?\n";
        std::cout << "1) verification\n";
        std::cout << "2) performance\n";
        std::cin >> sim_type;

        if (sim_type == 1) {
            nsteps = 1000;
        } else if (sim_type == 2) {
            nsteps = 1;
        } else {
            std::cerr << "Choose 1 or 2!\n";
            exit(EXIT_FAILURE);
        }

        std::cout << "Running with nstep = " << nsteps << "\n";

        //-------------------------------------------------------------------------------------------------
        // WD
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol1 = newSolution("WD/INPUT/WD.yaml");
        auto gas = sol1->thermo();
        gas->setState_TPY(1000.0, 100000.0, "CH4:0.20, O2:0.8");

        // Reactor setup
        auto reactor = newReactorBase("IdealGasReactor", sol1);
        ReactorNet sim(reactor);
        //sim.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 8.e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp;
        timeTemp.reserve(nsteps);
        double time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim.advance(time);
            double T = reactor->temperature();
            timeTemp.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = t1 - t0;
        std::cout << "WD Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("WD", elapsed.count());

        // Write results for first simulation
        std::ofstream outFile("WD/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp) {
            outFile << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Troyes
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol2 = newSolution("Troyes/INPUT/troyes.yaml");
        auto gas2 = sol2->thermo();

        gas2->setState_TPY(1000.0, 100000.0,
            "H2:0.00534, CL2:0.18796, N2:0.80670");

        auto reactor2 = newReactorBase("IdealGasReactor", sol2);
        ReactorNet sim2(reactor2);
        sim2.setTolerances(1e-12, 1e-15);

        dt = 5.0e-3 / nsteps;

        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp2;
        timeTemp2.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim2.advance(time);
            double T = reactor2->temperature();
            timeTemp2.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "Troyes Cantera-CXX time = "
                << std::scientific << elapsed.count() << std::endl;

        summaryData.emplace_back("Troyes", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile2("Troyes/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp2) {
            outFile2 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Ecker
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol3 = newSolution("Ecker/INPUT/ecker.yaml");
        auto gas3 = sol3->thermo();
        gas3->setState_TPY(1000.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        auto reactor3 = newReactorBase("IdealGasReactor", sol3);
        ReactorNet sim3(reactor3);
        sim3.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp3;
        timeTemp3.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim3.advance(time);
            double T = reactor3->temperature();
            timeTemp3.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "Ecker Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;

        summaryData.emplace_back("Ecker", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile3("Ecker/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp3) {
            outFile3 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Cross
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol4 = newSolution("Cross/INPUT/cross.yaml");
        auto gas4 = sol4->thermo();
        gas4->setState_TPY(1010.0, 100000.0, "H2:0.00534534, CL2:0.18798856, N2:0.8066661");

        // Reactor setup
        auto reactor4 = newReactorBase("IdealGasReactor", sol4);
        ReactorNet sim4(reactor4);
        sim4.setTolerances(1e-7, 1e-7);
        sim4.setMaxSteps(100000);

        // Integration settings
        dt = 1.0e-2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp4;
        timeTemp4.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim4.advance(time);
            double T = reactor4->temperature();
            timeTemp4.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "Cross Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("Cross", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile4("Cross/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp4) {
            outFile4 << time << "\t" << temp << "\n";
        }


        //-------------------------------------------------------------------------------------------------
        // Smooke
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol5 = newSolution("Smooke/INPUT/smooke.yaml");
        auto gas5 = sol5->thermo();
        gas5->setState_TPY(1300.0, 100000.0, "CH4:0.0552, O2:0.2201, N2: 0.7247");

        // Reactor setup
        auto reactor5 = newReactorBase("IdealGasReactor", sol5);
        ReactorNet sim5(reactor5);
        sim5.setTolerances(1e-7, 1e-7);
        sim5.setMaxSteps(1000000);

        // Integration settings
        dt = 0.2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp5;
        timeTemp5.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim5.advance(time);
            double T = reactor5->temperature();
            timeTemp5.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "Smooke Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("Smooke", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile5("Smooke/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp5) {
            outFile5 << time << "\t" << temp << "\n";
        }


        //-------------------------------------------------------------------------------------------------
        // CORIA-CNRS
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol6 = newSolution("CORIA/INPUT/coria.yaml");
        auto gas6 = sol6->thermo();
        gas6->setState_TPY(1300.0, 100000.0, "CH4:0.2, O2:0.8");

        // Reactor setup
        auto reactor6 = newReactorBase("IdealGasReactor", sol6);
        ReactorNet sim6(reactor6);
        sim6.setTolerances(1e-7, 1e-7);
        sim6.setMaxSteps(1000000);

        // Integration settings
        dt = 0.005 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp6;
        timeTemp6.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim6.advance(time);
            double T = reactor6->temperature();
            timeTemp6.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "CORIA Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("CORIA", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile6("CORIA/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp6) {
            outFile6 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-CDF-13
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol7 = newSolution("TSR-CDF-13/INPUT/TSR-CDF-13.yaml");
        auto gas7 = sol7->thermo();
        gas7->setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        auto reactor7 = newReactorBase("IdealGasReactor", sol7);
        ReactorNet sim7(reactor7);
        sim7.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp7;
        timeTemp7.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim7.advance(time);
            double T = reactor7->temperature();
            timeTemp7.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "TSR-CDF-13 Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;

        summaryData.emplace_back("TSR-CDF-13", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile7("TSR-CDF-13/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp7) {
            outFile7 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Pelucchi
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol8 = newSolution("Pelucchi/INPUT/pelucchi.yaml");
        auto gas8 = sol8->thermo();
        gas8->setState_TPY(1250.0, 100000, "CO:0.00859, O2:0.00606, H2O:0.00365, HCL:0.00025, N2:0.98044");

        // Reactor setup
        auto reactor8 = newReactorBase("IdealGasReactor", sol8);
        ReactorNet sim8(reactor8);
        sim8.setTolerances(1e-12, 1e-15);

        // Integration settings
        dt = 5.0e-2 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp8;
        timeTemp8.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim8.advance(time);
            double T = reactor8->temperature();
            timeTemp8.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "Pelucchi Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;

        summaryData.emplace_back("Pelucchi", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile8("Pelucchi/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp8) {
            outFile8 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // ZK
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol9 = newSolution("ZK/INPUT/ZK.yaml");
        auto gas9 = sol9->thermo();
        gas9->setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        auto reactor9 = newReactorBase("IdealGasReactor", sol9);
        ReactorNet sim9(reactor9);
        sim9.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 2.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp9;
        timeTemp9.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim9.advance(time);
            double T = reactor9->temperature();
            timeTemp9.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "ZK Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("ZK", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile9("ZK/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp9) {
            outFile9 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-GP-24
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol10 = newSolution("TSR-GP-24/INPUT/TSR-GP-24.yaml");
        auto gas10 = sol10->thermo();
        gas10->setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        auto reactor10 = newReactorBase("IdealGasReactor", sol10);
        ReactorNet sim10(reactor10);
        sim10.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 2.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp10;
        timeTemp10.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim10.advance(time);
            double T = reactor10->temperature();
            timeTemp10.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "TSR-GP-24 Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("TSR-GP-24", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile10("TSR-GP-24/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp10) {
            outFile10 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // TSR-Rich-31
        //-------------------------------------------------------------------------------------------------

        // Create the solution and thermo object
        auto sol11 = newSolution("TSR-Rich-31/INPUT/TSR-Rich-31.yaml");
        auto gas11 = sol11->thermo();
        gas11->setState_TPY(1300.0, 500000, "CH4:1, O2:4");

        // Reactor setup
        auto reactor11 = newReactorBase("IdealGasReactor", sol11);
        ReactorNet sim11(reactor11);
        sim11.setTolerances(1e-7, 1e-7);

        // Integration settings
        dt = 5.0e-3 / nsteps;

        // Run the simulation
        t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<double, double>> timeTemp11;
        timeTemp11.reserve(nsteps);
        time = 0.0;
        for (int i = 0; i < nsteps; i++) {
            time += dt;
            sim11.advance(time);
            double T = reactor11->temperature();
            timeTemp11.emplace_back(time, T);
        }
        t1 = std::chrono::high_resolution_clock::now();

        elapsed = t1 - t0;
        std::cout << "TSR-Rich-31 Cantera-CXX time = "
          << std::scientific << elapsed.count() << std::endl;
        summaryData.emplace_back("TSR-Rich-31", elapsed.count());

        // Write results for second simulation
        std::ofstream outFile11("TSR-Rich-31/OUTPUT/batch-CXX.dat");
        for (const auto& [time, temp] : timeTemp11) {
            outFile11 << time << "\t" << temp << "\n";
        }

        //-------------------------------------------------------------------------------------------------
        // Summary output
        //-------------------------------------------------------------------------------------------------


        std::ofstream summaryFile("comp-batch-cantera.dat");
        for (const auto& [nspec, time] : summaryData) {
            summaryFile << nspec << "\t" << std::scientific << time << "\n";
        }

        appdelete();
        return 0;

    } 
    
    catch (CanteraError& err) {
    std::cout << err.what() << std::endl;
    appdelete();
    return -1;}
}
