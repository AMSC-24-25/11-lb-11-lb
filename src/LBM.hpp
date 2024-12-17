#ifndef LBM_HPP
#define LBM_HPP


class LBM {
private:
    // D2Q9
    static constexpr int D = 2;
    static constexpr int Q = 9;      
    
    // Velocity configuration for the D2Q9 model
    static constexpr int e[Q * D] = {1,1,  1,0,  1,-1,  0,1,  0,0,  0,-1,  -1,1,  -1,0,  -1,-1};

    // Equilibrium function weights for each direction
    static constexpr double w[Q] = {
        1.0/36.0,
        1.0/9.0,
        1.0/36.0,
        1.0/9.0,
        4.0/9.0,
        1.0/9.0,
        1.0/36.0,
        1.0/9.0,
        1.0/36.0
    };

    const double u_lid; // Velocity imposed on the upper boundary

    // Main domain variables
    double  *rho, 
            *rho2,
            *u,
            *u2,
            *f, 
            *F;

    double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;

    static inline const int& direction(unsigned int d, unsigned int v){
        return e[d*D + v];
    }

    inline double& density(unsigned int x, unsigned int y){
        return rho[(NY)*x+y];
    }

    inline double& density_2(unsigned int x, unsigned int y){
        return rho2[(NY)*x+y];
    }

    inline double& velocity(unsigned int x, unsigned int y, unsigned int d){
        return u[((NY)*x+y)*D + d];
    }

    inline double& velocity_2(unsigned int x, unsigned int y, unsigned int d){
        return u2[((NY)*x+y)*D + d];
    }

    inline double& field(unsigned int x, unsigned int y,unsigned int d){
        return f[((NY)*x+y)*Q + d];
    }

    inline double& field_2(unsigned int x, unsigned int y,unsigned int d){
        return F[((NY)*x+y)*Q + d];
    }

    double feq(unsigned int k, unsigned int x, unsigned int y);

    
    void compute();  
    void apply_boundary_conditions();

public : 

    const int NX; // Dimension in the x-direction
    const int NY; // Dimension in the y-direction

    LBM(unsigned int nx, unsigned int ny, double u_lid, double Re);
    ~LBM();

    void evolution();
    void evolution(unsigned int iterations);

    inline const double& get_vel(unsigned int x, unsigned int y, unsigned int d) const {
        return u[((NY)*x+y)*D + d];
    }

}; 
#endif