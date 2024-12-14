#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

const int Q = 9;      // Modello D2Q9
const int NX = 100;   // Dimensione in direzione x
const int NY = 100;   // Dimensione in direzione y
const double U = 0.50; // Velocità imposta al bordo superiore

// Configurazione delle velocità per il modello D2Q9
int e[Q][2] = {
    {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};

// Pesi della funzione di equilibrio per ogni direzione
double w[Q] = {
    4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

// Variabili principali del dominio
double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], u0[NX + 1][NY + 1][2],
    f[NX + 1][NY + 1][Q], F[NX + 1][NY + 1][Q];

// Variabili globali per iterazioni e parametri del modello
int i, j, k, ip, jp, n, P, Z;
double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;

// Dichiarazione delle funzioni principali
void init();
double feq(int k, double rho, double u[2]);
void evolution();

int main() {
    int maxSteps = 100000; // Numero massimo di passi temporali

    // Creazione del file di output per velocità
    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "Creazione di vel_data.txt fallita.\n";
        return 1;
    }
    file << NX << "\n" << NY << "\n";

    init(); // Inizializzazione del dominio

    for (n = 0; n < maxSteps; n++) {
        evolution(); // Evoluzione del sistema
		for (i=0; i<=NX; i++)
			u[i][NY][0] = U;

        // Ogni 100 passi, salva i dati di velocità
        if (n % 1000 == 0) {
            for (int i = 0; i <= NX; ++i) {
                for (int j = 0; j <= NY; ++j) {
                    double vx = u[i][j][0];
                    double vy = u[i][j][1];
                    double v = sqrt(vx * vx + vy * vy);
                    file << v << "\n";
                }
            }
        }

        // Aggiorna la barra di progresso
        if (n % 50 == 0 || n == maxSteps - 1) {
            float progress = (static_cast<float>(n + 1) / maxSteps) * 100.0f;
            std::cout << "\rProgresso: " << std::fixed << std::setprecision(2)
                      << progress << "% completato" << std::flush;
        }
    }

    return 0;
}

void init() {
    dx = 1.0;         // Passo spaziale
    dy = 1.0;         // Passo spaziale
    Lx = dx * double(NY); // Lunghezza del dominio in y
    Ly = dy * double(NX); // Lunghezza del dominio in x
    dt = dx;          // Passo temporale
    c = dx / dt;      // Velocità del suono (in unità di reticolo)
    rho0 = 1.0;       // Densità iniziale
    Re = 1000;        // Numero di Reynolds
    niu = U * Lx / Re; // Viscosita cinematica
    tau_f = 3.0 * niu + 0.5; // Tempo di rilassamento

    // Inizializzazione delle variabili nel dominio
    for (i = 0; i <= NX; i++)
        for (j = 0; j <= NY; j++) {
            u[i][j][0] = 0; // Velocità iniziale in x
            u[i][j][1] = 0; // Velocità iniziale in y
            rho[i][j] = rho0; // Densità iniziale
			u[i][NY][0] = U;  // Velocità imposta al bordo superiore

            for (k = 0; k < Q; k++) {
                f[i][j][k] = feq(k, rho[i][j], u[i][j]); // Funzione di equilibrio iniziale
            }
			u[i][NY][0] = U;  // Velocità imposta al bordo superiore
        }
}

double feq(int k, double rho, double u[2]) { // Calcolo della funzione di equilibrio
    double eu, uv, feq;
    eu = (e[k][0] * u[0] + e[k][1] * u[1]); // Prodotto scalare tra velocità microscopica e macroscopica
    uv = (u[0] * u[0] + u[1] * u[1]);       // Modulo quadrato della velocità macroscopica
    feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
    return feq;
}

void evolution() {
    // Evoluzione della funzione di distribuzione
    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++)
            for (k = 0; k < Q; k++) {
                ip = i - e[k][0]; // Nodo da cui proviene il contributo
                jp = j - e[k][1];
                F[i][j][k] = f[ip][jp][k] +
                             (feq(k, rho[ip][jp], u[ip][jp]) - f[ip][jp][k]) / tau_f; // Collisione
            }

    // Calcolo delle grandezze macroscopiche
    for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++) {
            u0[i][j][0] = u[i][j][0];
            u0[i][j][1] = u[i][j][1];
            rho[i][j] = 0;
            u[i][j][0] = 0;
            u[i][j][1] = 0;

            for (k = 0; k < Q; k++) {
                f[i][j][k] = F[i][j][k];
                rho[i][j] += f[i][j][k];
                u[i][j][0] += e[k][0] * f[i][j][k];
                u[i][j][1] += e[k][1] * f[i][j][k];
            }

            u[i][j][0] /= rho[i][j];
            u[i][j][1] /= rho[i][j];
        }

	// **Condizioni al contorno**
    for (j = 0; j <= NY; j++) { // Bordo sinistro e destro
        // Bordo sinistro (x = 0)
        u[0][j][0] = 0.0;
        u[0][j][1] = 0.0;
        rho[0][j] = rho[1][j]; // Continuazione della densità
        for (k = 0; k < Q; k++) {
            f[0][j][k] = feq(k, rho[0][j], u[0][j]) + f[1][j][k] - feq(k, rho[1][j], u[1][j]);
        }

        // Bordo destro (x = NX)
        u[NX][j][0] = 0.0;
        u[NX][j][1] = 0.0;
        rho[NX][j] = rho[NX - 1][j]; // Continuazione della densità
        for (k = 0; k < Q; k++) {
            f[NX][j][k] = feq(k, rho[NX][j], u[NX][j]) + f[NX - 1][j][k] - feq(k, rho[NX - 1][j], u[NX - 1][j]);
        }
    }
	
    for (i = 0; i <= NX; i++) { // Bordo superiore e inferiore
        // Bordo superiore (y = NY)
        u[i][NY][0] = U;  // Velocità costante
        u[i][NY][1] = 0.0;
        rho[i][NY] = rho[i][NY - 1]; // Continuazione della densità
        for (k = 0; k < Q; k++) {
            f[i][NY][k] = feq(k, rho[i][NY], u[i][NY]) + f[i][NY - 1][k] - feq(k, rho[i][NY - 1], u[i][NY - 1]);
        }

        // Bordo inferiore (y = 0)
        u[i][0][0] = 0.0;
        u[i][0][1] = 0.0;
        rho[i][0] = rho[i][1]; // Continuazione della densità
        for (k = 0; k < Q; k++) {
            f[i][0][k] = feq(k, rho[i][0], u[i][0]) + f[i][1][k] - feq(k, rho[i][1], u[i][1]);
        }
    }
}

