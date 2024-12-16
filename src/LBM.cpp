#include "LBM.hpp"

#include <cstdlib>

LBM::LBM(unsigned int nx, unsigned int ny, double u_lid, double Re) : NX(nx), NY(ny), U(u_lid), Re(Re) {
    dx = 1.0;         // Spatial step
    dy = 1.0;         // Spatial step
    Lx = dx * double(NY); // Domain length in y
    Ly = dy * double(NX); // Domain length in x
    dt = dx;          // Time step
    c = dx / dt;      // Speed of sound (in lattice units)
    rho0 = 1.0;       // Initial density
    //Re = 1000;        // Reynolds number
    niu = U * Lx / Re; // Kinematic viscosity
    tau_f = 3.0 * niu + 0.5; // Relaxation time

    rho = (double*)malloc(sizeof(double)*(NX+1)*(NY+1));
    u = (double*)malloc(sizeof(double)*(NX+1)*(NY+1)*D);
    f = (double*)malloc(sizeof(double)*(NX+1)*(NY + 1)*Q);
    F = (double*)malloc(sizeof(double)*(NX+1)*(NY + 1)*Q);



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

LBM::~LBM() {
    free(rho);
    free(u);
    free(f);
    free(F);
}

// Compute the equilibrium function
double LBM::feq(unsigned int k, unsigned int x, unsigned int y) { 
    double eu{0.0}, uv{0.0};
    for(int a = 0; a < D; a++) {
        eu += direction(k,a) * velocity(x,y,a);  
        uv += velocity(x,y,a) * velocity(x,y,a);
    }

    return w[k] * density(x,y) * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
}

void LBM::evolution(unsigned int iterations) {
    for(int iter = 0; iter < iterations; iter++) {
        compute_collision();      // Collision
        compute_macro_quantities(); // Update macroscopic density and velocities
        apply_boundary_conditions(); // Apply boundary conditions
    }
}

// Collision
void LBM::compute_collision() {
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
void LBM::compute_macro_quantities() {
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
void LBM::apply_boundary_conditions() {
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

        // Bottom boundary (y = 0)
        for(k = 0; k < D; k++) velocity(i,0,k) = 0.0;

        velocity(i,0,0) = 0;

        density(i,0) = density(i,1);
        for (k = 0; k < Q; k++) {
            field(i,0,k) = feq(k, i, 0) + field(i,1,k) - feq(k, i, 1);
        }
    }
}