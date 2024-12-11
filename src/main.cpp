#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <iomanip> 

const size_t nx=30, ny=30;
const double Re = 100;
const double L = 1.0;
const int q = 9;

const double u_lid = 0.05;
const double rho_0 = 1.0;

const double nu = 0.1;
const double omega = 0.8;//1.0 / (3.0 * nu + 0.5); // 1/tau

const double cs2 = 1.0/3.0;

// velocity directions D2Q9
const int ex[q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

// weights D2Q9
const double w[q] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 
                     1.0/36, 1.0/36, 1.0/36, 1.0/36};

double equilibrium(int i, double rho, double ux, double uy) {
    double eu = ex[i] * ux + ey[i] * uy;
    double u2 = ux * ux + uy * uy;
    return w[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * u2);
}

void applyNoSlip(std::vector<std::vector<std::vector<double>>>& f, int nx, int ny) {
    // y = 0
    for (int x = 0; x < nx; ++x) {
        f[x][0][2] = f[x][0][4];
        f[x][0][5] = f[x][0][7];
        f[x][0][6] = f[x][0][8];
    }
    
    // y = ny-1
    for (int x = 0; x < nx; ++x) {
        f[x][ny-1][4] = f[x][ny-1][2];
        f[x][ny-1][7] = f[x][ny-1][5];
        f[x][ny-1][8] = f[x][ny-1][6];
    }

    for (int y = 0; y < ny; ++y) {
        // x = 0
        f[0][y][1] = f[0][y][3];
        f[0][y][5] = f[0][y][7];
        f[0][y][8] = f[0][y][6];

        // x = nx-1
        f[nx-1][y][3] = f[nx-1][y][1];
        f[nx-1][y][7] = f[nx-1][y][5];
        f[nx-1][y][6] = f[nx-1][y][8];
    }
}

void applyLidVelocity(std::vector<std::vector<std::vector<double>>>& f,
                      std::vector<std::vector<double>>& rho,
                      double u_lid, int nx, int ny) {
    for (int x = 0; x < nx; ++x) {
        double ux = u_lid;
        double uy = 0; 

        double rho_wall = f[x][ny-1][0] + f[x][ny-1][2] + f[x][ny-1][4] +
                          2 * (f[x][ny-1][1] + f[x][ny-1][5] + f[x][ny-1][8]);

        f[x][ny-1][3] = f[x][ny-1][1] - (2.0 / 3.0) * rho_wall * ux;
        f[x][ny-1][6] = f[x][ny-1][8] - (1.0 / 6.0) * rho_wall * ux;
        f[x][ny-1][7] = f[x][ny-1][5] - (1.0 / 6.0) * rho_wall * ux;
    }
}

void save_results(const std::vector<std::vector<double>>& rho,
                  const std::vector<std::vector<double>>& ux,
                  const std::vector<std::vector<double>>& uy,
                  const std::string& filename) {
    std::ofstream file(filename);
    file << "y,x,rho,ux,uy\n";
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            file << x << "," << y << "," << rho[x][y] << "," << ux[x][y] << "," << uy[x][y] << "\n";
        }
    }
    file.close();
}

std::vector<std::vector<std::vector<double>>> initialize() {
    std::vector<std::vector<std::vector<double>>> f(nx, 
        std::vector<std::vector<double>>(ny, 
        std::vector<double>(q, 0.0)));

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            double rho = 1.0;
            double ux = 0.0;
            double uy = 0.0;
            for (int i = 0; i < q; ++i) {
                f[x][y][i] = equilibrium(i, rho, ux, uy);
            }
        }
    }
    return f;
}

void collide(std::vector<std::vector<std::vector<double>>>& f, 
             std::vector<std::vector<double>>& rho, 
             std::vector<std::vector<double>>& ux, 
             std::vector<std::vector<double>>& uy) {
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int i = 0; i < q; ++i) {
                double feq = equilibrium(i, rho[x][y], ux[x][y], uy[x][y]);
                f[x][y][i] += -omega * (f[x][y][i] - feq);
            }
        }
    }
}

void stream(std::vector<std::vector<std::vector<double>>>& f) {
    std::vector<std::vector<std::vector<double>>> f_new = f;

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int i = 0; i < q; ++i) {
                int bx = (x + ex[i] + nx) % nx;
                int by = (y + ey[i] + ny) % ny;
                f_new[bx][by][i] = f[x][y][i];
            }
        }
    }
    f = f_new;
}

void computeRhoAndVelocity(const std::vector<std::vector<std::vector<double>>>& f,
                           std::vector<std::vector<double>>& rho,
                           std::vector<std::vector<double>>& ux,
                           std::vector<std::vector<double>>& uy) {
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            rho[x][y] = 0.0;
            ux[x][y] = 0.0;
            uy[x][y] = 0.0;
            for (int i = 0; i < q; ++i) {
                rho[x][y] += f[x][y][i];
                ux[x][y] += f[x][y][i] * ex[i];
                uy[x][y] += f[x][y][i] * ey[i];
            }
            ux[x][y] /= rho[x][y];
            uy[x][y] /= rho[x][y];
        }
    }
}

#include <ctime>

int main() {
    auto f = initialize();
    

    std::vector<std::vector<double>> rho(nx, std::vector<double>(ny, 1.0));
    std::vector<std::vector<double>> ux(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> uy(nx, std::vector<double>(ny, 0.0));

    const int maxSteps = 1000;

    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "vel_data.txt creation failed.\n";
        return 1;
    }

    file << nx << "\n" << ny << "\n";

    for (int t = 0; t < maxSteps; ++t) {
        stream(f);

        applyNoSlip(f, nx, ny);
        applyLidVelocity(f, rho, u_lid, nx, ny);
        computeRhoAndVelocity(f, rho, ux, uy);

        collide(f, rho, ux, uy);

        
        if(t%10 == 0) {
            for (int i = 0; i < ny; ++i) {
                for (int j = 0; j < nx; ++j) {
                    double vx = ux[j][i]; 
                    double vy = uy[j][i];
                    double v = sqrt(vx*vx + vy*vy); 
                    file << v << "\n";
                }
            }
        }

        if(t%50 == 0 || t == maxSteps-1) {
            float progress = (static_cast<float>(t+1) / maxSteps) * 100.0f;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2)
                    << progress << "% completed" << std::flush;   
        }
    }

    file.close();

    std::cout << "\n";

    save_results(rho, ux, uy, "output.csv");

    return 0;
}

// #include "LBM.hpp"
// int main(){
//     size_t nx=50, ny=50;
//     double Re = 100;
//     double L = 1.0;
//     double u_lid = 1.0;
//     double nu = u_lid*L/Re;
//     double tau = 3.0*nu + 0.5;

//     LBM<2> solver(nx,ny,1, tau);

//     size_t max_steps = 100;
//     for (size_t step = 0; step < max_steps; step++)
//     {
//         solver.simulate();

//         // Boundary conditions
//         for (size_t x=0; x<nx ; x++)
//         {
//             solver.setBoundaryVelocity( x , ny-1, {u_lid, 0.0});
//         }

//         // Error and printing check
//         int percent = (step * 100) / max_steps;
//         if (percent % 5 == 0) {
//             std::cout << percent << "% completato.\n";
//         }
//         solver.writeVTK("output.vtk");
//     }
//     return 0;
//     //g++ -std=c++17 -O2 -Wall *.cpp -o main
// }