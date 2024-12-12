#include "LBM.hpp"

// consrtuctor
template <size_t dim>
LBM<dim>::LBM(size_t nx_, size_t ny_, size_t nz_, double tau_)
    : nx(nx_), ny(ny_), nz(nz_), tau(tau_) {
    initialize();
}

template <size_t dim>
void LBM<dim>::initialize()
{
    size_t total_nodes = nx*ny*(dim == 3 ? nz : 1);
    velocities = generateVelocities();
    f.resize(total_nodes, std::vector<double> (velocities.size(), 1.0));
    f_temp = f;
    rho.resize(total_nodes, 1.0);
    u.resize(total_nodes, {0.0});
}

template <size_t dim>
std::vector<std::array<double, dim>> LBM<dim>::generateVelocities()
{
    if constexpr (dim==2)
        return {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    else if constexpr (dim == 3)
        return {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {-1,0,0}, {0,-1,0}, {0,0,-1}};
    else 
    {
        std::cout << "\n!DIMENTIONS ARE NOT COMPATIBLE!" << std::endl;
        return 0;
    }
}

template <size_t dim>
void LBM<dim>::simulate() {
    // new f temp to update f_temp during the collision
    std::vector<std::vector<double>> new_f_temp(f.size(), 
                                                std::vector<double> (velocities.size(), 0.0));
    // Colllision has to be considered first
    for (size_t node; node < f.size(); node++)
    {
        // Macroscopic variables have to be evaluated at each time step starting from 0
        rho[node] = 0.0;
        u[node] ={1.0};

        // Evaluation
        for (size_t i=0; i<velocities.size(); i++)
        {
            rho[node] += f[node][i];
            for (size_t d=0; d<dim ; d++)
            {
                u[node][d] = f[node][i]*velocities[i][d];
            }
        }
        // Velocity has to be normalized
        for (size_t d = 0; d < dim; ++d) {
            u[node][d] /= rho[node];
        }
        for (size_t i=0; i<velocities.size() ; i++)
        {
            double feq = equilibrium( rho[node], u[node], velocities[i]);
            new_f_temp[node][i] = f[node][i] + (feq - f[node][i]) / tau;
        }

        // Streaming part
        for (size_t x = 0; x<nx ; x++){
            for (size_t y=0; y<ny; y++){
                for (size_t z=0; z< (dim==3 ? nz : 1); z++){
                    size_t node = x+nx*(y+ny*(z));

                     // Streaming for each direction
                    for (size_t i = 0; i < velocities.size(); ++i) {
                        int dx = velocities[i][0];
                        int dy = velocities[i][1];
                        int dz = (dim == 3) ? velocities[i][2] : 0;

                        // source node
                        size_t source_x = (x + nx - dx) % nx;
                        size_t source_y = (y + ny - dy) % ny;
                        size_t source_z = (dim == 3) ? (z + nz - dz) % nz : 0;

                        size_t source_node = source_x + nx * (source_y + ny * source_z);

                        // copy distribution after collision
                        f[node][i] = new_f_temp[source_node][i];
                    }
                }
            }
        }

    }
}

template <size_t dim>
double LBM<dim>::equilibrium(double rho, const std::array<double, dim> &u, const std::array<double, dim> &v)
{
    double uv = dotProduct(u, v);
    double uu = dotProduct(u, u);
    return rho * (1.0 + 3.0 * uv + 4.5 * uv * uv - 1.5 * uu);
}

template <size_t dim>
double LBM<dim>::dotProduct(const std::array<double, dim>& a, const std::array<double, dim>& b) {
    double result = 0.0;
    for (size_t i = 0; i < dim; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template <size_t dim>
void LBM<dim>::setBoundaryVelocity(size_t x, size_t y, std::array<double, dim> velocity) {
    size_t node = x + nx * y;
    u[node] = velocity;

    // update f to respect boundary conditions
    for (size_t i = 0; i < velocities.size(); ++i) {
        f[node][i] = equilibrium(rho[node], velocity, velocities[i]);
    }
}

// Funzione per esportare i dati in formato VTK
template <size_t dim>
void LBM<dim>::writeVTK(const std::string& filename) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Errore nell'aprire il file VTK!" << std::endl;
        return;
    }

    // Header del file VTK
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "LBM Data\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << nx << " " << ny << " " << (dim == 3 ? nz : 1) << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING 1 1 1\n";
    vtkFile << "POINT_DATA " << nx * ny * (dim == 3 ? nz : 1) << "\n";

    // Esportazione della densità
    vtkFile << "SCALARS density double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < nx * ny * (dim == 3 ? nz : 1); ++i) {
        vtkFile << rho[i] << "\n";
    }

    // Esportazione della velocità (composto di componenti X, Y e Z)
    vtkFile << "VECTORS velocity double\n";
    for (size_t i = 0; i < nx * ny * (dim == 3 ? nz : 1); ++i) {
        vtkFile << u[i][0] << " " << u[i][1];
        if (dim == 3) {
            vtkFile << " " << u[i][2];
        }
        vtkFile << "\n";
    }

    vtkFile.close();
}

template class LBM<2>; // for 2D
template class LBM<3>; // for 3D