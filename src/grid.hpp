#include <vector>
#include <array>

template <typename T>
class Grid {
public:
    Grid(int Nx, int Ny) : Nx(Nx), Ny(Ny) {
        f.resize(Ny, std::vector<std::array<T, 9>>(Nx));
        rho.resize(Ny, std::vector<T>(Nx, 1.0));
        ux.resize(Ny, std::vector<T>(Nx, 0.0));
        uy.resize(Ny, std::vector<T>(Nx, 0.0));
    }

    // Accesso alle funzioni di distribuzione
    std::array<T, 9>& getDistribution(int x, int y) {
        return f[y][x];
    }

    // Accesso a densità e velocità
    T& density(int x, int y) { return rho[y][x]; }
    T& velocityX(int x, int y) { return ux[y][x]; }
    T& velocityY(int x, int y) { return uy[y][x]; }

    int getWidth() const { return Nx; }
    int getHeight() const { return Ny; }

private:
    int Nx, Ny;
    std::vector<std::vector<std::array<T, 9>>> f; // Funzioni di distribuzione
    std::vector<std::vector<T>> rho, ux, uy;      
};
