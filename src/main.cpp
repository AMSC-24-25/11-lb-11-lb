#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <iomanip> 

const double Re = 100;
const double L = 1.0; //SCALE
const unsigned int nx=30*L, ny=nx;
const int q = 9;

const size_t m_size_ndir = sizeof(double) * nx*ny*q;
const size_t m_size_scalar = sizeof(double) * nx * ny ;

const double u_lid = 0.05;
const double rho_0 = 1.0;

const double nu = 1.0/6.0;

const double tau = (3.0 * nu + 0.5);
const double omega = 1/tau;

const double cs2 = 1.0/3.0;

//LATTICE WEIGHTS
const double w0=4.0/9.0;
const double ws=1.0/9.0;
const double wd= 1.0/36.0;

//direction component
const double w[q] ={w0,ws,ws,ws,ws,wd,wd,wd,wd};

const int ex[q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

inline size_t scalar_index(unsigned int x, unsigned int y){
    return nx*y+x;
}

inline size_t field_index(unsigned int x, unsigned int y,unsigned int d){
    return nx*(ny*d+y)+x;
}

void init_equilibrium(double *f, double *r, double *u, double *v){
    for(unsigned int y=0;y<ny;++y){
        for(unsigned int x=0;x<nx;++x){
            double rho=r[scalar_index(x,y)];
            double ux = u[scalar_index(x,y)];
            double uy = v[scalar_index(x,y)];

            for(unsigned int i = 0; i<q;++i){
                double cdotu=ex[i]*ux + ey[i]*uy;
                f[field_index(x,y,i)]= w[i]*rho*(1.0 + 3.0*cdotu + 4.5*cdotu*cdotu - 1.5*(ux*ux*uy*uy));
            }
        }

    }
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



void collide(double *f, double *r , double *u , double *v) {
    for (unsigned int y = 0; y < ny; ++y) {
        for (unsigned int x = 0; x < nx; ++x) {
            double rho = r[scalar_index(x,y)];
            double ux = u[scalar_index(x,y)];
            double uy = v[scalar_index(x,y)];

            for (unsigned int i = 0; i < q; ++i) {
                double cdotu=ex[i]*ux + ey[i]*uy;
                double feq = w[i]*rho*(1.0 + 3.0*cdotu + 4.5*cdotu*cdotu - 1.5*(ux*ux*uy*uy));
                f[field_index(x,y,i)] = (1-omega) * f[field_index(x,y,i)] + omega * feq ;
            }
        }
    }
}

void stream(double *f_src , double *f_dst) {
   
    for (int y= 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            for (int i = 0; i < q; ++i) {
                int xmd = (nx+x-ex[i]) % nx;
                int ymd = (ny+y-ey[i]) % ny;
                f_dst[field_index(x,y,i)]= f_src[field_index(xmd,ymd,i)];
            }
        }
    }
}

void compute_rho_u(double *f, double *r , double *u , double *v) {
    for (int y= 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;
            for (int i = 0; i < q; ++i) {
                rho+=f[field_index(x,y,i)];
                ux+= f[field_index(x,y,i)]*ex[i];
                uy+= f[field_index(x,y,i)]*ey[i];
            }
            r[scalar_index(x,y)] = rho ;
            u[scalar_index(x,y)] = ux/rho;
            v[scalar_index(x,y)] = uy/rho;

         }
    }
}

#include <ctime>

int main() {
    
    double *f1 =(double*) malloc(m_size_ndir);
    double *f2 = (double*) malloc(m_size_ndir);
    double *rho = (double*) malloc(m_size_scalar);
    double *ux = (double*) malloc(m_size_scalar);
    double *uy = (double*) malloc(m_size_scalar);

    const int maxSteps = 1000;

    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "vel_data.txt creation failed.\n";
        return 1;
    }

    file << nx << "\n" << ny << "\n";

    init_equilibrium(f1,rho,ux,uy);

    for (int t = 0; t < maxSteps; ++t) {
        stream(f1,f2);

        applyNoSlip(f, nx, ny);
        applyLidVelocity(f, rho, u_lid, nx, ny);


        compute_rho_u(f2,rho,ux,uy);
        collide(f2, rho, ux, uy);

        double *temp = f1;
          f1=f2;
          f2=temp;

        
        if(t%10 == 0) {
            for (int i = 0; i < ny; ++i) {
                for (int j = 0; j < nx; ++j) {
                    double vx = ux[scalar_index(j,i)]; 
                    double vy = uy[scalar_index(j,i)];
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

    free(f1);
    free(f2);
    free(rho);
    free(ux);
    free(uy);


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