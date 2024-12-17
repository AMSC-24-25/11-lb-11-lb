#include "LBM.hpp"

#include <cstdlib>

LBM::LBM(unsigned int nx, unsigned int ny, double u_lid, double Re) : NX(nx), NY(ny), u_lid(u_lid), Re(Re) {
    dx = 1.0;         // Spatial step
    dy = 1.0;         // Spatial step
    Lx = dx * double(NY); // Domain length in y
    Ly = dy * double(NX); // Domain length in x
    dt = dx;          // Time step
    c = dx / dt;      // Speed of sound (in lattice units)
    rho0 = 1.0;       // Initial density
    //Re = 1000;        // Reynolds number
    niu = u_lid * Lx / Re; // Kinematic viscosity
    tau_f = 3.0 * niu + 0.5; // Relaxation time

    rho = (double*)malloc(sizeof(double)*(NX)*(NY));
    rho2 = (double*)malloc(sizeof(double)*(NX)*(NY));
    u = (double*)malloc(sizeof(double)*(NX)*(NY)*D);
    u2 = (double*)malloc(sizeof(double)*(NX)*(NY)*D);
    f = (double*)malloc(sizeof(double)*(NX)*(NY)*Q);
    F = (double*)malloc(sizeof(double)*(NX)*(NY)*Q);



    // Initialize variables in the domain
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for(int k = 0; k < D; k++) {
                velocity(i,j,k) = 0; // Initial velocity
            }

            density(i,j) = rho0; // Initial density
			velocity(i, NY, 0) = u_lid;  // Velocity imposed on the upper boundary

            for (int k = 0; k < Q; k++) {
                field(i,j,k) = feq(k, i, j); // Initial equilibrium function
            }
        }
    }
    
}

LBM::~LBM() {
    free(rho);
    free(rho2);
    free(u);
    free(u2);
    free(f);
    free(F);
}

// Compute the equilibrium function
double LBM::feq(unsigned int k, unsigned int x, unsigned int y) {  
    double eu = direction(k,0) * velocity(x,y,0) + direction(k,1) * velocity(x,y,1);  
    double uv = velocity(x,y,0) * velocity(x,y,0) + velocity(x,y,1) * velocity(x,y,1);

    return w[k] * density(x,y) * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
}

void LBM::evolution() {
    compute();      // Collision and Update macroscopic density and velocities
    apply_boundary_conditions(); // Apply boundary conditions
}


void LBM::evolution(unsigned int iterations) {
   
    #pragma omp parallel for 
    for(int iter = 0; iter < iterations; iter++) {
        compute();      // Collision and Update macroscopic density and velocities
        apply_boundary_conditions(); // Apply boundary conditions
    }
}

// Collision and Update macroscopic density and velocities
void LBM::compute() {
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
    
            density_2(i,j) = 0;

            velocity_2(i,j,0) = 0;
            velocity_2(i,j,1) = 0;

            for (int k = 0; k < Q; k++) {
                int ip = i - direction(k,0); // Node from which the contribution comes
                int jp = j - direction(k,1);

                field_2(i,j,k) = field(ip,jp,k) +(feq(k, ip, jp) - field(ip,jp,k)) / tau_f; // Collision

                density_2(i,j) += field_2(i,j,k);
            
                velocity_2(i,j,0) += direction(k,0) * field_2(i,j,k);
                velocity_2(i,j,1) += direction(k,1) * field_2(i,j,k);
            }

            velocity_2(i,j,0) /= density_2(i,j);
            velocity_2(i,j,1) /= density_2(i,j);
        }
    }

    double *t;

    t = u;
    u = u2;
    u2 = t;

    t = rho;
    rho = rho2;
    rho2 = t;

    t = f;
    f = F;
    F = t;
}


// Apply boundary conditions
void LBM::apply_boundary_conditions() {
    for (int j = 0; j < NY; j++) { // Left and right boundaries
        // Left boundary (x = 0)
        velocity(0,j,0) = 0.0;
        velocity(0,j,1) = 0.0;
        
        density(0,j) = density(1,j);
        for (int k = 0; k < Q; k++) {
            field(0,j,k) = feq(k, 0, j) + field(1,j,k) - feq(k, 1, j);
        }

        // Right boundary (x = NX-1)
        velocity(NX-1,j,0) = 0.0;
        velocity(NX-1,j,1) = 0.0;

        density(NX-1,j) = density(NX-2,j);
        for (int k = 0; k < Q; k++) {
            field(NX-1,j,k) = feq(k, NX-1, j) + field(NX-2,j,k) - feq(k, NX-2, j);
        }
    }
	
    for (int i = 0; i < NX; i++) { // Bottom and top boundaries
        velocity(i,NY-1,0) = u_lid;
        velocity(i,NY-1,1) = 0.0;

        density(i,NY-1) = density(i,NY - 2);
        for (int k = 0; k < Q; k++) {
            field(i,NY-1,k) = feq(k, i, NY-1) + field(i,NY-2,k) - feq(k, i, NY-2);
        }

        // Bottom boundary (y = 0)
        velocity(i,0,1) = 0.0;
        velocity(i,0,2) = 0.0;

        density(i,0) = density(i,1);
        for (int k = 0; k < Q; k++) {
            field(i,0,k) = feq(k, i, 0) + field(i,1,k) - feq(k, i, 1);
        }
    }
}