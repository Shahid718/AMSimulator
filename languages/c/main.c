#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main () {
    // declare data
    const int Nx = 100, Ny = 100, nsteps = 3000;
    int istep, i, j, ip, im, jp, jm;
    float dx, dy, delta, sigma, Tm, alpha, aniso, theta0, K, Cp, H, mu, T0,
        lambda, b, a0, w, pmob, term1, term2, f, theta, dt, t1, t2;
    float p[Nx][Ny], pnew[Nx][Ny], T[Nx][Ny], Tnew[Nx][Ny], pdx[Nx][Ny], pdy[Nx][Ny];
    float a[Nx][Ny], da[Nx][Ny], T_lap[Nx][Ny], p_lap[Nx][Ny];

    // initialize data
    dx = 20e-9;
    dy = 20e-9;
    delta = 4.0 * dx;
    sigma = 0.37;
    alpha = 0.05;
    aniso = 4.0;
    theta0 = 0.0;
    Tm = 1728.0;
    K = 84.01;
    Cp = 5.42e+06;
    H = 2.35e+09;
    mu = 2.0;
    T0 = Tm - H / Cp * 0.3;
    lambda = 0.1;
    b = 2.0 * atanh(1.0 - 2.0 * lambda);
    a0 = sqrt(3.0 * delta * sigma / b);
    w = 6.0 * sigma * b / delta;
    pmob = b * Tm * mu / (3.0 * delta * H);
    dt = dx * dx / (5.0 * K / Cp);
		
    // get the start time
    clock_t start_time = clock();
	
    // initial microstructure
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            T[i][j] = T0;
            p[i][j] = 0;
            if ((i-Nx/2.0)*(i-Nx/2.0)+(j-Ny/2.0)*(j-Ny/2.0)-3<= 0) {
                p[i][j] = 1.0;
            }
        }
    }
	
    // evolution 
    for ( istep = 0; istep < nsteps ; istep++ ) {
        for (i = 0; i < Nx ; i++) {
            for (j = 0; j < Ny ; j++) {

                // boundary conditions
                jp = (j + 1)      % Ny;
                jm = (j - 1 + Ny) % Ny;
                ip = (i + 1)      % Nx;
                im = (i - 1 + Nx) % Nx; 
                //laplacian
                p_lap[i][j] = (p[ip][j] + p[im][j] + p[i][jm] + p[i][jp] - 4.0 * p[i][j]) / (dx * dy);
                T_lap[i][j] = (T[ip][j] + T[im][j] + T[i][jm] + T[i][jp] - 4.0 * T[i][j]) / (dx * dy);
                // gradients and angle          
                pdx[i][j] = (p[ip][j] - p[im][j]) / (2.0 * dx);
                pdy[i][j] = (p[i][jp] - p[i][jm]) / (2.0 * dy);
                theta = atan2(pdy[i][j], pdx[i][j]);
                // anisotropy and its derivative
                a[i][j] = a0 * (1.0 + alpha * cos(aniso * (theta - theta0)));
                da[i][j] = -a0 * aniso * alpha * sin(aniso * (theta - theta0));
                term1 = (a[i][jp] * da[i][jp] * pdx[i][jp] - a[i][jm] * da[i][jm] * pdx[i][jm]) / dy;
                term2 = -(a[ip][j] * da[ip][j] * pdy[ip][j] - a[im][j] * da[im][j] * pdy[im][j]) / dx;
                // Driving force
                f = -H * (T[i][j] - Tm) / Tm;
                // Phi and temperature
                pnew[i][j] = p[i][j] + (a[i][j] * a[i][j] * p_lap[i][j] + term1 + term2 +
                    4.0 * w * p[i][j] * (1.0 - p[i][j]) * (p[i][j] - 0.5 + 15.0 / (2.0 * w) * f *
                        p[i][j] * (1.0 - p[i][j]))) * dt * pmob;
                Tnew[i][j] = T[i][j] + (K / Cp) * T_lap[i][j] * dt + 30.0 * p[i][j] * p[i][j] *
                    (1.0 - p[i][j]) * (1.0 - p[i][j]) * (H / Cp) * (pnew[i][j] - p[i][j]);
            }
        }
		// update fields
  		for (i = 0; i < Nx ; i++) {
			for (j = 0; j < Ny ; j++) {
				p[i][j] = pnew[i][j];
				T[i][j] = Tnew[i][j];
			} 
		}            
   }      
    // get the end time
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf(" code time  %f  seconds\n", elapsed_time);
    return 0;
}