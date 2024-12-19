#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip> 
#include <omp.h>

#include "LBM.hpp"



const int maxSteps = 10000; // Maximum number of time steps
const double Re = 10000;
const double u_lid = 0.5;

const int MAX_THREADS = 12;
const int ITERATIONS_PER_FRAME = 10;
const int ITERATIONS_PER_PROGRESS_UPDATE = 100;

const int NX = 100; // Dimension in the x-direction
const int NY = 100; // Dimension in the y-direction

int main() {
    // Create the output file for velocity
    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "could not opene/create 'vel_data.txt'.\n";
        return 1;
    }
    file << NX << "\n" << NY << "\n";

    for(int threadsCount = 1; threadsCount <= MAX_THREADS; threadsCount++) {
        LBM lbm(NX, NY, u_lid, Re);
        omp_set_num_threads(threadsCount);

        auto startTime = std::chrono::high_resolution_clock::now();
        lbm.evolution(maxSteps);
        auto endTime = std::chrono::high_resolution_clock::now();

        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        std::cout << "Threads: " << threadsCount << " - " << "Time: " << elapsedTime << "ms" << std::endl;
    }

    file.close();

    std::cout << "\n";

    return 0;
}
