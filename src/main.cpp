#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip> 
#include <omp.h>

#include "LBM.hpp"



const int maxSteps = 100; // Maximum number of time steps
const double Re = 10000;
const double u_lid = 0.5;

const int MAX_THREADS = 12;
const int WORK_PER_THREAD = 20000;

const int ITERATIONS_PER_FRAME = 10;
const int ITERATIONS_PER_PROGRESS_UPDATE = 100;


int main() {

    for(int threadsCount = 1; threadsCount <= MAX_THREADS; threadsCount++) {
        int NX = sqrt(WORK_PER_THREAD * threadsCount);
        int NY = NX;
        while (NX*NY < WORK_PER_THREAD*threadsCount) NY++;
        
        LBM lbm(NX, NY, u_lid, Re);
        omp_set_num_threads(threadsCount);

        auto startTime = std::chrono::high_resolution_clock::now();
        lbm.evolution(maxSteps);
        auto endTime = std::chrono::high_resolution_clock::now();

        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        std::cout << "Threads: " << threadsCount << " - " << "Cavity size: " << NY << "x" << NX << " - " << "Time: " << elapsedTime << "ms" << std::endl;
    }

    std::cout << "\n";

    return 0;
}
