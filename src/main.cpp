#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <array>
#include <iomanip> 
#include <thread>

#include <matplot/matplot.h>

const double Re = 100;
const double L = 1.0; //SCALE
const unsigned int nx=100*L, ny=nx;
const int q = 9;

const size_t m_size_ndir = sizeof(double) * nx*ny*q;
const size_t m_size_scalar = sizeof(double) * nx * ny ;

const double u_lid = 0.5;
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
    return (nx*y+x)*q + d;
}

void init_equilibrium(double *f, double *r, double *u, double *v){
    for(unsigned int y=0;y<ny;++y){
        for(unsigned int x=0;x<nx;++x){
            double rho= rho_0;//r[scalar_index(x,y)];
            double ux = 0;// u[scalar_index(x,y)];
            double uy = 0;//v[scalar_index(x,y)];

            for(unsigned int i = 0; i<q;++i){
                double cdotu=ex[i]*ux + ey[i]*uy;
                f[field_index(x,y,i)]= w[i]*rho*(1.0 + 3.0*cdotu + 4.5*cdotu*cdotu - 1.5*(ux*ux + uy*uy));
            }
        }

    }
}

void collide(double *f, double *r , double *u , double *v) {
    const double tauinv = 1/tau;//2.0/(6.0*nu+1.0);
    const double omtauinv = 1.0-tauinv;
    for (unsigned int y = 0; y < ny; ++y) {
        for (unsigned int x = 0; x < nx; ++x) {
            double rho = r[scalar_index(x,y)];
            double ux = u[scalar_index(x,y)];
            double uy = v[scalar_index(x,y)];

            for (unsigned int i = 0; i < q; ++i) {
                double cdotu=ex[i]*ux + ey[i]*uy;
                double feq = w[i]*rho*(1.0 + 3.0*cdotu + 4.5*cdotu*cdotu - 1.5*(ux*ux + uy*uy));
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

void boundary_conditions(double *f, double *ux, double*uy, double*r)
{
    const double tauinv = 1/tau;//2.0/(6.0*nu+1.0);
    const double omtauinv = 1.0-tauinv;
    for (unsigned int x=0; x<nx;x++)
    {
        for(unsigned int y=0; y<ny; y++)
        {
            ux[scalar_index(x,0)] = u_lid;
            ux[scalar_index(0,y)] = 0.0;
            ux[scalar_index(x,ny-1)] = 0.0;
            ux[scalar_index(nx-1, y)] = 0.0;            
            double rho = r[scalar_index(x,y)];
            double u = ux[scalar_index(x,y)];
            double v = uy[scalar_index(x,y)];

            for (unsigned int i = 0; i < q; ++i) {
                double cdotu=ex[i]*u + ey[i]*v;
                double feq = w[i]*rho*(1.0 + 3.0*cdotu + 4.5*cdotu*cdotu - 1.5*(u*u + v*v));
                f[field_index(x,y,i)] = (1-omega) * f[field_index(x,y,i)] + omega * feq ;
            }
        } 
    }
}


int main() {
    
    double *f1 =(double*) malloc(m_size_ndir);
    double *f2 = (double*) malloc(m_size_ndir);
    double *rho = (double*) malloc(m_size_scalar);
    double *ux = (double*) malloc(m_size_scalar);
    double *uy = (double*) malloc(m_size_scalar);

    const int maxSteps = 10000;

    for(int i = 0; i < nx*ny; ++i){
        rho[i] = rho_0;
    }


    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "vel_data.txt creation failed.\n";
        return 1;
    }

    file << nx << "\n" << ny << "\n";

    init_equilibrium(f1,rho,ux,uy);

    
    using namespace matplot;
    auto f = figure(true);

    
    std::vector<std::vector<double>> velocities(ny, std::vector<double>(nx, 0.0));
    //auto heatmap = imshow(velocities);
    image(velocities, true);
    colorbar();

    //show();

    for (int t = 0; t < maxSteps; ++t) {

        boundary_conditions(f1,ux,uy, rho);
        stream(f1,f2);
        compute_rho_u(f2,rho,ux,uy);
        collide(f2, rho, ux, uy);

        double *temp = f1;
        f1=f2;
        f2=temp;

        if(t%50 == 0) {
            
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    double vx = ux[scalar_index(i,j)]; 
                    double vy = uy[scalar_index(i,j)];
                    double v = sqrt(vx*vx + vy*vy); 
                    //file << v << "\n";
                    velocities[j][i] = v;
                }
            }
            // Ricrea la heatmap
            image(velocities, true);
            colorbar();
            f->draw();

            // Introduci un ritardo
            std::this_thread::sleep_for(std::chrono::milliseconds(300));
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