#ifndef LBM_HPP
#define LBM_HPP

#include "vector"
#include "array"
#include "cmath"
#include "iostream"

template <size_t dim> class LBM {
    private:
        size_t nx, ny, nz; // Domain dimentions
        double tau;        // Relaxing time 
        std::vector<std::array<double, dim>> velocities, u; // Microscopic and macroscopic velocities
        std::vector<std::vector<double>> f, f_temp; // Population distribution
        std::vector<double> rho;                    // Macroscopic density

    public : 
        // Costructor
        LBM(size_t nx_, size_t ny_, size_t nz_, double tau_);
        
        void initialize();

        std::vector<std::array<double, dim>> generateVelocities();

        void simulate();

        double equilibrium(double rho, const std::array<double, dim>& u, const std::array<double, dim>& v);

        double dotProduct(const std::array<double, dim>& a, const std::array<double, dim>& b);
    
        void setBoundaryVelocity(size_t x, size_t y, std::array<double, dim> velocity);


}; 
#endif