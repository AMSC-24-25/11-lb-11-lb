#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <iomanip> 
#include <thread>


const int maxSteps = 20000; // Maximum number of time steps

const int ITERATIONS_PER_FRAME = 25;
const int ITERATIONS_PER_PROGRESS_UPDATE = 100;

// D2Q9
const int D = 2;
const int Q = 9;      

const int NX = 100; // Dimension in the x-direction
const int NY = 100; // Dimension in the y-direction
const double U = 0.50; // Velocity imposed on the upper boundary

// Velocity configuration for the D2Q9 model
int e[Q * D] = {0, 0, 1, 0, 0, 1,-1, 0,0, -1,1, 1,-1, 1,-1, -1,1, -1};

// Equilibrium function weights for each direction
double w[Q] = {
    4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

// Main domain variables
double rho[(NX+1)*(NY+1)], u[(NX+1)*(NY+1)*D],
    f[(NX+1)*(NY + 1)*Q], F[(NX+1)*(NY + 1)*Q];

// Global variables for iterations and model parameters
int i, j, k, ip, jp, n, P, Z;
double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;


inline int& direction(unsigned int d, unsigned int v){
    return e[d*D + v];
}

inline double& density(unsigned int x, unsigned int y){
    return rho[(NX+1)*y+x];
}

inline double& velocity(unsigned int x, unsigned int y, unsigned int d){
    return u[((NX+1)*y+x)*D + d];
}

inline double& field(unsigned int x, unsigned int y,unsigned int d){
    return f[((NX+1)*y+x)*Q + d];
}

inline double& field_2(unsigned int x, unsigned int y,unsigned int d){
    return F[((NX+1)*y+x)*Q + d];
}

void init();
double feq(unsigned int k, unsigned int x, unsigned int y);

void evolution();
void compute_collision();  
void compute_macro_quantities();
void apply_boundary_conditions();

int main() {
    // Create the output file for velocity
    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "Creazione di vel_data.txt fallita.\n";
        return 1;
    }
    file << NX << "\n" << NY << "\n";

    init(); // initialization

    auto startTime = std::chrono::high_resolution_clock::now();

    for (n = 1; n <= maxSteps; n++) {
        evolution(); // System evolution

        //Every ITERATIONS_PER_FRAME steps, save velocity data
        if (n % ITERATIONS_PER_FRAME == 0) {
            for (int i = 0; i <= NX; ++i) {
                for (int j = 0; j <= NY; ++j) {
                    double vx = velocity(i,j,0); 
                    double vy = velocity(i,j,1);
                    double v = sqrt(vx*vx + vy*vy); 
                    file << v << "\n";
                }
            }
        }

        // Update the progress bar
        if (n % ITERATIONS_PER_PROGRESS_UPDATE == 0 || n == maxSteps - 1) {
            float progress = (static_cast<float>(n) / maxSteps);
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            
            double estimatedTotalTime = elapsedTime / progress;
            int remainingTime = estimatedTotalTime - elapsedTime;

            progress *= 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "% completed "
                      << "| Elapsed Time: " << elapsedTime << "s, "
                      << "Remaining Time (estimated): " << static_cast<int>(remainingTime) << "s"
                      << std::flush;
        }
    }

    file.close();

    std::cout << "\n";

    return 0;
}

void init() {
    dx = 1.0;         // Spatial step
    dy = 1.0;         // Spatial step
    Lx = dx * double(NY); // Domain length in y
    Ly = dy * double(NX); // Domain length in x
    dt = dx;          // Time step
    c = dx / dt;      // Speed of sound (in lattice units)
    rho0 = 1.0;       // Initial density
    Re = 10000;        // Reynolds number
    niu = U * Lx / Re; // Kinematic viscosity
    tau_f = 3.0 * niu + 0.5; // Relaxation time

    // Initialize variables in the domain
    for (i = 0; i <= NX; i++) {
        for (j = 0; j <= NY; j++) {
            for(k = 0; k < D; k++) {
                velocity(i,j,k) = 0; // Initial velocity
            }

            density(i,j) = rho0; // Initial density
			velocity(i, NY, 0) = U;  // Velocity imposed on the upper boundary

            for (k = 0; k < Q; k++) {
                field(i,j,k) = feq(k, i, j); // Initial equilibrium function
            }
			velocity(i,NY,0) = U; // Velocity imposed on the upper boundary
        }
    }
}

// Compute the equilibrium function
double feq(unsigned int k, unsigned int x, unsigned int y) { 
    double eu{0.0}, uv{0.0};
    for(int a = 0; a < D; a++) {
        eu += direction(k,a) * velocity(x,y,a);  
        uv += velocity(x,y,a) * velocity(x,y,a);
    }

    return w[k] * density(x,y) * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
}

void evolution() {
    compute_collision();      // Collision
    compute_macro_quantities(); // Update macroscopic density and velocities
    apply_boundary_conditions(); // Apply boundary conditions
}

// Collision
void compute_collision() {
    for (i = 1; i < NX; i++) {
        for (j = 1; j < NY; j++) {
            for (k = 0; k < Q; k++) {
                ip = i - direction(k,0); // Node from which the contribution comes
                jp = j - direction(k,1);

                field_2(i,j,k) = field(ip,jp,k) +(feq(k, ip, jp) - field(ip,jp,k)) / tau_f; // Collision
            }
        }
    }
}

// Compute macroscopic quantities
void compute_macro_quantities() {
    for (i = 0; i < NX; i++) {
        for (j = 0; j < NY; j++) {
            density(i,j) = 0;

            for(k = 0; k < D; k++) velocity(i,j,k) = 0;
            
            for (k = 0; k < Q; k++) {
                field(i,j,k) = field_2(i,j,k);
                density(i,j) += field(i,j,k);
            
                velocity(i,j,0) += direction(k,0) * field(i,j,k);
                velocity(i,j,1) += direction(k,1) * field(i,j,k);
            }

            for(k = 0; k < D; k++) velocity(i,j,k) /= density(i,j);
        }
    }
}

// Apply boundary conditions
void apply_boundary_conditions() {
    for (j = 0; j <= NY; j++) { // Left and right boundaries
        // Left boundary (x = 0)
        for(k = 0; k < D; k++) velocity(0,j,k) = 0.0;
        
        density(0,j) = density(1,j);
        for (k = 0; k < Q; k++) {
            field(0,j,k) = feq(k, 0, j) + field(1,j,k) - feq(k, 1, j);
        }

        // Right boundary (x = NX)
        for(k = 0; k < D; k++) velocity(NX,j,k) = 0.0;

        density(NX,j) = density(NX-1,j);
        for (k = 0; k < Q; k++) {
            field(NX,j,k) = feq(k, NX, j) + field(NX-1,j,k) - feq(k, NX-1, j);
        }
    }
	
    for (i = 0; i <= NX; i++) { // Bottom and top boundaries
        velocity(i,NY,0) = U; 
        for(k = 1; k < D; k++) velocity(i,NY,k) = 0.0;

        density(i,NY) = density(i,NY - 1);
        for (k = 0; k < Q; k++) {
            field(i,NY,k) = feq(k, i, NY) + field(i,NY-1,k) - feq(k, i, NY-1);
        }

        // Top boundary (y = NY)
        for(k = 0; k < D; k++) velocity(i,0,k) = 0.0;

        density(i,0) = density(i,1);
        for (k = 0; k < Q; k++) {
            field(i,0,k) = feq(k, i, 0) + field(i,1,k) - feq(k, i, 1);
        }
    }
}