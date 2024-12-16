#ifndef LBM_HPP
#define LBM_HPP


class LBM {
private:
    // D2Q9
    static constexpr int D = 2;
    static constexpr int Q = 9;      
    
    // Velocity configuration for the D2Q9 model
    static constexpr int e[Q * D] = {0, 0, 1, 0, 0, 1,-1, 0,0, -1,1, 1,-1, 1,-1, -1,1, -1};

    // Equilibrium function weights for each direction
    static constexpr double w[Q] = {
        4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36
        };

    const double U = 0.50; // Velocity imposed on the upper boundary

    // Main domain variables
    double  *rho, 
            *u,
            *f, 
            *F;

    // Global variables for iterations and model parameters
    int i, j, k, ip, jp, n, P, Z;
    double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;

    static inline const int& direction(unsigned int d, unsigned int v){
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

    double feq(unsigned int k, unsigned int x, unsigned int y);

    
    void compute_collision();  
    void compute_macro_quantities();
    void apply_boundary_conditions();

public : 

    const int NX; // Dimension in the x-direction
    const int NY; // Dimension in the y-direction

    LBM(unsigned int nx, unsigned int ny, double u_lid, double Re);
    ~LBM();

    void evolution(unsigned int iterations=1);

    inline const double& get_vel(unsigned int x, unsigned int y, unsigned int d) const {
        return u[((NX+1)*y+x)*D + d];
    }

}; 
#endif