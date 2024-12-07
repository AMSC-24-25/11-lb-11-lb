#include "LBM.hpp"
int main(){
    size_t nx=50, ny=50;
    double Re = 100;
    double L = 1.0;
    double u_lid = 1.0;
    double nu = u_lid*L/Re;
    double tau = 3.0*nu + 0.5;

    LBM<2> solver(nx,ny,1, tau);

    size_t max_steps = 100;
    for (size_t step = 0; step < max_steps; step++)
    {
        solver.simulate();

        // Boundary conditions
        for (size_t x=0; x<nx ; x++)
        {
            solver.setBoundaryVelocity( x , ny-1, {u_lid, 0.0});
        }

        // Error and printing check
        int percent = (step * 100) / max_steps;
        if (percent % 5 == 0) {
            std::cout << percent << "% completato.\n";
        }
        solver.writeVTK("output.vtk");
    }
    return 0;
    //g++ -std=c++17 -O2 -Wall *.cpp -o main
}