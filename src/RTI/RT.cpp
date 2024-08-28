#include <iostream>
#include <cmath>

// In this version, all the intermediate variables are calculated for clarity. The code is not optimized for performance.
// Periodic boundary conditions p170 of the book (Graduate Texts in Physics) Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, Erlend Magnus Viggen (auth.) - The Lattice Boltzmann Method_ Principles and Practice-Springer
// Initial conditions depend on the problem

// periodic boundary conditions on x directon
// half-way bounce-back boundary conditions on y direction with resting wall 

// Define your lattice dimensions and parameters
const int Lx = 128; // Number of lattice nodes in the x-direction
const int Ly = 1028; // Number of lattice nodes in the y-direction
const int Q = 9;  // Number of discrete velocities


// Define your simulation parameters without dimensions
const double dx = 1.; // Lattice spacing
const double dt = 1.; // Time step

const double c_dx = 1.e-5 ; // m : conversion factor from lattice units to physical units
const double c_dt = c_dx/347./sqrt(3.); // s : conversion factor from lattice units to physical units

const int numSteps = 5000000; // Number of simulation steps
const int interval = 10000; // Output interval
const double epsilon = 1.0e-10; // Small number to avoid division by zero in recoloration step 

// Speed of sound and related constants
const double cs = dx / sqrt(3.0) / dt; // Speed of sound in the lattice
const double cs2 = cs * cs;
const double cs4 = cs2 * cs2;
const double cs6 = cs4 * cs2;

const double a_g = 9.81e2/(c_dx/c_dt/c_dt); // m/s^2
//parameters
const double rho1 = 4.; // kg * m-3 Density for component 1 
const double rho2 = 1.;  // kg * m-3 Density for component 2
const double c1 = 347./(c_dx/c_dt);//dx / sqrt(3.0) / dt; // Speed of sound for component 1 in lattice units
const double c2 = 347./(c_dx/c_dt);//dx / sqrt(3.0) / dt; // Speed of sound for component 2 in lattice units
const double radius = 10.; // Radius of the droplet
const double sigma = 0/(c_dx * c_dx * c_dx / c_dt / c_dt); //1.e-3 / (347*347*3*1.e-4);// 1e-3; //s1.e-3; // surface tension
//const double sigma = (rho1 - rho2) * cs2 * radius;
const double ch_width_init = 1.1 * dx; // Characteristic width of the interface
const double ch_width_ope = 1.6 * dx;

// viscosities to calculate Relaxation times in collide step
const double nu   = 1.e-4/(c_dx * c_dx / c_dt) ; //1.0e-3*sqrt(3)/(347*1.0e-3); // kinematic viscosity in lattice units p284 of the book (Graduate Texts in Physics) Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, Erlend Magnus Viggen (auth.) - The Lattice Boltzmann Method_ Principles and Practice-Springer
const double nu_b = 1.e-4/(c_dx * c_dx / c_dt) ; //1.0e-3*sqrt(3)/(347*1.0e-3); // bulk viscosity in lattice units

//The pressure at infinity is used at the equation of state to calculate the pressure and equilibrium distribution function
const double p1_inf = rho1 * c1 * c1 - rho2 * c2 * c2 - sigma / radius; // Pressure at infinity for component 1
const double p2_inf = 0.; // Pressure at infinity for component 2

double S[Lx][Ly][Q]; // Force term
double F[Lx][Ly][2]; // External volumic force : kg * m^-2 * s^-2
//for Laplace equation test, the gravity is on the perpendicular direction to the interface: 0
//for instability test, the gravity is on the parallel direction to the interface

double rho[Lx][Ly]; // Density
double rho_mdt[Lx][Ly]; // Density at the previous time step
double u[Lx][Ly][2]; // Velocity components
double p[Lx][Ly]; // Pressure
double p_mdt[Lx][Ly]; // Pressure at the previous time step
double phi[Lx][Ly]; // Phase field function

double omega_1[Lx][Ly][Q]; // Collision term
double omega_2[Lx][Ly][Q]; // Collision term related to surface tension
double omega_3[Lx][Ly][Q]; // Recoloring term
double f[Lx][Ly][Q]; // Sum of distribution functions
double g[Lx][Ly][Q]; // Difference of distribution functions
double f_eq[Lx][Ly][Q]; // Sum of distribution functions at equilibrium

const double xi[Q][2] = {
    {0, 0},
    {1, 0}, {0, 1}, {-1, 0}, {0, -1},
    {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
};

//weights
const double w1 = 4./9.;
const double w2 = 1./9.;
const double w3 = 1./36.;
const double w[Q] = {w1, w2, w2, w2, w2, w3, w3, w3, w3};



void calEquilibrium() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double rho_local = rho[i][j];  // Local density
            double u_x = u[i][j][0];
            double u_y = u[i][j][1];
            double p_local = p[i][j];
            for (int k = 0; k < Q; k++) {
                double H0 = 1.0;
                double Hx = xi[k][0];
                double Hy = xi[k][1];

                double Hxx = Hx * Hx - cs2 ; // xi[k][0]*xi[k][0] - cs2;
                double Hyy = Hy * Hy - cs2 ; // xi[k][1]*xi[k][1] - cs2;
                double Hxy = xi[k][0]*xi[k][1]; // Hyx = Hxy

                double Hxxy = Hx * Hx * Hy - cs2 * Hy;
                double Hyyx = Hy * Hy * Hx - cs2 * Hx;
                double Hxxx = pow(Hx, 3) - cs2 * 3. * Hx;
                double Hyyy = pow(Hy, 3) - cs2 * 3. * Hy;

                double Hxxyy = Hx * Hx * Hy * Hy - cs2 * (Hx * Hx + Hy * Hy) + cs4;

                double E = w[k]*((Hxx + Hyy) / (2.*cs4) - Hxxyy / (4.*cs6));
                double term1 = rho_local * w[k] * (H0 + u_x * Hx / cs2 + u_y * Hy / cs2 + 0.5 * (u_x * u_x * Hxx + u_x * u_y * Hxy * 2. + u_y * u_y * Hyy) / cs4);
                f_eq[i][j][k] = term1 + (p_local - rho_local * cs2) * (E +  w[k] * (u_x * (Hyyx + Hxxx) + u_y * (Hyyy + Hxxy))/(2.*cs6));
                // f_eq[i][j][k] = term1 + (p_local - rho_local * cs2) * (E +  w[k] * (u_x * (Hyyx) + u_y * (Hxxy))/(2.*cs6));
                // std::cout << "k = " << k << " f_eq = " << f_eq[i][j][k] << std::endl;
            }
        }
    }
}

void calMacroscopic() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double sum_f = 0.0;
            double sum_xi_x = 0.0;
            double sum_xi_y = 0.0;
            for (int k = 0; k < Q; k++) {
                double f_local = f[i][j][k];
                sum_f += f_local;
                sum_xi_x += f_local * xi[k][0];
                sum_xi_y += f_local * xi[k][1];
            }
            rho[i][j] = sum_f;
            F[i][j][0] = 0.0; // External force in lattice units on x direction
            F[i][j][1] = - sum_f * a_g; // External force in lattice units on y direction
            u[i][j][0] = (sum_xi_x + F[i][j][0]*dt*0.5) / sum_f;//add force term Guo al. 
            //u[i][j][1] = sum_xi_y / sum_f;
            u[i][j][1] = (sum_xi_y + F[i][j][1]*dt*0.5) / sum_f;

        }
    }
}

// Function to calculate the phase field function
void calPhaseField() {
    double c1_squared = c1 * c1;
    double c2_squared = c2 * c2;
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double sum_f = 0.0;
            double sum_g = 0.0;
            for (int k = 0; k < Q; k++) {
                sum_f += f[i][j][k];
                sum_g += g[i][j][k];
            }
            phi[i][j] = sum_g / sum_f;
            // if (abs(phi[i][j]) > 1.0) {
            //     std::cout << "Error : phi = " << phi[i][j] << std::endl;
            // }
            double rho_local = rho[i][j];
            double phi_local = sum_g / sum_f;
            double c_hat_squared = (c1_squared + c2_squared) * 0.5 + phi_local * (c1_squared - c2_squared) * 0.5;
            double c_bar_squared = (c1_squared - c2_squared) * 0.5 + phi_local * (c1_squared + c2_squared) * 0.5;
            double sqrt_term = sqrt(pow((p2_inf - p1_inf + rho_local * c_bar_squared), 2) + rho_local * rho_local * (1.0 - phi_local * phi_local) * c1_squared * c2_squared);
            // if (abs(phi[i][j]) > 1.0) {
            //     std::cout << "sqrt_term = " << sqrt_term << std::endl;
            // }
            double p_local = 0.5 * (rho_local * c_hat_squared - p1_inf - p2_inf + sqrt_term);
            p_local = rho_local * cs2 - (1+phi_local)/2.*p1_inf - (1-phi_local)/2.*p2_inf;
            p[i][j] = p_local;
            }
    }
}

// Function to perform the collision step
void collide() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double rho_local = rho[i][j];  // Local density
            double p_local = p[i][j];  // Local pressure
            // Calculate relaxation times
            double tau_nu = rho_local * nu   / (p_local * dt) + 0.5;  //shear relaxation time
            double tau_b  = rho_local * nu_b / (p_local * dt) + 0.5; //bulk relaxation time
            double sum_nu_neq = 0.0, sum_b_neq = 0.0, sum_xy_neq = 0.0;

            // Calculate $f_{k, i}^{r, neq}$ with $k\inn\{\nu, b, x y\}$
            for (int k = 0; k < Q; k++) {
                double xi_x = xi[k][0], xi_y = xi[k][1];
                double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2;
                double H_xy = xi_x * xi_y;

                double f_neq = f[i][j][k] - f_eq[i][j][k] + 0.5*S[i][j][k]; // Calculate non-equilibrium part for each direction
                sum_nu_neq += f_neq * H_nu;
                sum_b_neq  += f_neq * H_b;
                sum_xy_neq += f_neq * H_xy;
            }
            // Calculate $\Omega_i^{(1)}$
            for (int k = 0; k < Q; k++) {
                double xi_x = xi[k][0], xi_y = xi[k][1];
                double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2;
                double H_xy = xi_x * xi_y;

                double f_nu_neq = H_nu / cs4 * sum_nu_neq; //optimal way to calculate f_nu_neq without using vector, only using scalar
                double f_b_neq  = H_b  / cs4 * sum_b_neq;
                double f_xy_neq = H_xy / cs4 * sum_xy_neq;
                omega_1[i][j][k] = w[k]*(1.0 - 1.0 / tau_nu) * (f_nu_neq + f_xy_neq) + w[k]*(1.0 - 1.0 / tau_b) * f_b_neq;            
            }
        }
    }
}

double S_F[Lx][Ly][Q]; // Force term for the potential volume forces eg.gravity
double S_Sp[Lx][Ly][Q]; // corrective term corrects the error stemming from the third order moment, improperly resolved on the D2Q9 lattice:
double S_t[Lx][Ly][Q]; // temporal correction term
// Function to calculate the force term
void force() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 1; j < Ly; j++) {
            double u_x = u[i][j][0];
            double u_y = u[i][j][1];
            double F_x = F[i][j][0];
            double F_y = F[i][j][1];           
            for (int k = 0; k < Q; k++) {
                double xik0 = xi[k][0];
                double xik1 = xi[k][1];
                double Hxxk = xik0*xik0 - cs2;
                double Hyyk = xik1*xik1 - cs2;
                double Hxyk = xik0*xik1;

                double term1 = (F_x * xik0 + F_y * xik1) / cs2;
                double term23 = (u_x * F_x * Hxxk + u_y * F_y * Hyyk + (u_x * F_y + u_y * F_x) * Hxyk) / cs4;

                S_F[i][j][k] = w[k] * (term1 + term23);
            }           
            double derive_x = 0., derive_y = 0.;
            for (int k = 0; k < Q; k++) {
                int xi_x = (int)xi[k][0];
                int xi_y = (int)xi[k][1];
                int ind_x = (i + xi_x + Lx) % Lx; // periodic boundary conditions
                int ind_y ;
                if (j==0 && (k==4 || k==7 || k==8)) {
                    derive_x += 0.;
                }
                else if (j==Ly-1 && (k==2 || k==5 || k==6)) {
                    derive_x += 0.;
                }
                else {
                    ind_y = j + xi_y;
                    derive_x += w[k] * xi_x * (p[ind_x][ind_y] - rho[ind_x][ind_y] * cs2) * u[ind_x][ind_y][0];
                    derive_y += w[k] * xi_y * (p[ind_x][ind_y] - rho[ind_x][ind_y] * cs2) * u[ind_x][ind_y][1];
                }
            }
            derive_x /= (dt*cs2);
            derive_y /= (dt*cs2);
            for (int k = 0; k < Q; k++) {
                double H_nu = (xi[k][0] * xi[k][0] - xi[k][1] * xi[k][1]) / 2.;
                double H_b  = (xi[k][0] * xi[k][0] + xi[k][1] * xi[k][1]) / 2. -  cs2 ;
                S_Sp[i][j][k] = w[k] * (derive_y*(3.*H_nu-H_b) + derive_x*(-3.*H_nu-H_b)) / (2.*cs4);
            }
            for (int k = 0; k < Q; k++) {
                double H_xx = xi[k][0] * xi[k][0] - cs2;
                double H_yy = xi[k][1] * xi[k][1] - cs2;
                double H_xxyy = xi[k][0] * xi[k][0] * xi[k][1] * xi[k][1] - cs2 * (xi[k][0] * xi[k][0] + xi[k][1] * xi[k][1]) + cs4;
                double E = w[k] * ((H_xx+H_yy) / (2.*cs4) - H_xxyy / (4.*cs6)) ; //w[k] * ((H_xx+H_yy) / (2.*cs4) - H_xxyy / (4.*cs6));
                S_t[i][j][k] = (p[i][j] - p_mdt[i][j] - (rho[i][j] - rho_mdt[i][j]) * cs2) * E;
                S[i][j][k] = S_F[i][j][k] + S_Sp[i][j][k] + S_t[i][j][k];
            }
        }
    }
}

void collide_surface(){
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            // Calculate color gradients Cx and Cy
            double Cx = 0.0, Cy = 0.0;
            for (int k = 0; k < Q; k++) {
                int ip = (i + (int)xi[k][0] + Lx) % Lx;  // Periodic boundary conditions on x
                int jp = j + (int)xi[k][1];
                if (j==0 && (k==4 || k==7 || k==8)) {
                    Cx += 0.;
                }
                else if (j==Ly-1 && (k==2 || k==5 || k==6)) {
                    Cx += 0.;
                }
                else {
                    Cx += w[k] * xi[k][0] * phi[ip][jp];
                    Cy += w[k] * xi[k][1] * phi[ip][jp];
                }
            }
            Cx *= 1./cs2; // Approximation of the spatial gradient
            Cy *= 1./cs2; // Approximation of the spatial gradient

            double norm_C = sqrt(Cx * Cx + Cy * Cy); // Magnitude of the color gradient vector
            for (int k = 0; k < Q; k++) {
                double xi_x = xi[k][0], xi_y = xi[k][1];
                double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2;
                double H_xy = xi_x * xi_y;
                double tau_nu = rho[i][j] * nu / (p[i][j] * dt) + 0.5;  //shear relaxation time
                double tau_b  = rho[i][j] * nu_b / (p[i][j] * dt) + 0.5; //bulk relaxation time
                if (norm_C > epsilon) { // Prevent division by zero
                    omega_2[i][j][k] = sigma * w[k] / (4 * norm_C * cs4) * ((2 * Cx * Cy * H_xy + (Cx * Cx - Cy * Cy) * H_nu) / tau_nu -((Cx * Cx + Cy * Cy) * H_b) / tau_b);
                }
                else {
                    omega_2[i][j][k] = 0.0; 
                }
            }

        }
    }
}

void recolor(){
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            // Calculate the gradient of the phase field by calculating the color gradient firstly
            // Calculate color gradients Cx and Cy
            double Cx = 0.0, Cy = 0.0;
            for (int k = 0; k < Q; k++) {
                int ip = (i + (int)xi[k][0] + Lx) % Lx;  // Periodic boundary conditions on x
                int jp = j + (int)xi[k][1];
                if (j==0 && (k==4 || k==7 || k==8)) {
                    Cx += 0.;
                }
                else if (j==Ly-1 && (k==2 || k==5 || k==6)) {
                    Cx += 0.;
                }
                else {
                    Cx += w[k] * xi[k][0] * phi[ip][jp];
                    Cy += w[k] * xi[k][1] * phi[ip][jp];
                }
            }
            Cx *= 1./cs2; // Approximation of the spatial gradient
            Cy *= 1./cs2; // Approximation of the spatial gradient
            double grad_phi_x, grad_phi_y;
            grad_phi_x = Cx / dt; 
            grad_phi_y = Cy / dt; 
            double norm_grad_phi = sqrt(grad_phi_x * grad_phi_x + grad_phi_y * grad_phi_y); // Magnitude of the gradient vector
            if (norm_grad_phi > epsilon) { // Prevent division by zero
                for (int k = 0; k < Q; k++) {
                        double xi_x = xi[k][0],  xi_y = xi[k][1];
                        omega_3[i][j][k] = w[k] * p[i][j] * (1 - phi[i][j] * phi[i][j]) / (2. * ch_width_ope) * (xi_x * grad_phi_x + xi_y * grad_phi_y) / (cs2 * norm_grad_phi); ; //w[k] * p[i][j] * (1 - phi[i][j] * phi[i][j]) / (2. * ch_width) * (xi_x * grad_phi_x + xi_y * grad_phi_y) / (epsilon);
                    }
            }
            else {
                for (int k = 0; k < Q; k++) {
                    omega_3[i][j][k] = 0.0;
                }
            }
        }
    }
}

// int oppositeDirection(int k) {
//     static int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
//     return opp[k];
// }


// Function to perform streaming step
void stream() {
    for (int i=0 ; i<Lx ; i++){
        // std::cout << "i = " << i << std::endl;
        for (int j=0 ; j<Ly ; j++){
            for (int k=0; k<Q ; k++){
                    int ip = (i + (int)xi[k][0] + Lx) % Lx;
                    int jp, kp;
                    if (j==0 && (k==4 || k==7 || k==8)) {
                        jp = j;
                        kp = k -2;
                    }
                    else if (j==Ly-1 && (k==2 || k==5 || k==6)) {
                        jp = j;
                        kp = k + 2;
                    }
                    else {
                        jp = j + (int)xi[k][1];
                        kp = k;
                    } 
                    f[ip][jp][kp] = f_eq[i][j][k] + omega_1[i][j][k] + omega_2[i][j][k] + 0.5*S[i][j][k];
                    g[ip][jp][kp] = f[ip][jp][k] * phi[i][j] + omega_3[i][j][k]; 
                }
            // Save the macroscopic variables at the previous time step for the force step
            rho_mdt[i][j] = rho[i][j];
            p_mdt[i][j] = p[i][j];
        }
    }
}

#include <fstream>
//#include <iostream> // Already included in the head
#include <string>
// Function to output data in VTK format
void outputVTK(const std::string& filename, int timestep) {
    std::ofstream vtkfile;
    std::string fullFilename = filename + std::to_string(timestep) + ".vtk";
    vtkfile.open(fullFilename);

    // VTK file header and data structure
    vtkfile << "# vtk DataFile Version 4.0" << std::endl;
    vtkfile << "Lattice Boltzmann Method Data" << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_POINTS" << std::endl;
    vtkfile << "DIMENSIONS " << Lx << " " << Ly << " 1" << std::endl;
    vtkfile << "ORIGIN 0 0 0" << std::endl;
    vtkfile << "SPACING 1 1 1" << std::endl;

    // Output density
    vtkfile << "POINT_DATA " << Lx * Ly << std::endl;
    vtkfile << "SCALARS density float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << rho[i][j] << std::endl;
        }
    }

    // Output velocity
    vtkfile << "VECTORS velocity float" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << u[i][j][0] << " " << u[i][j][1] << " 0.0" << std::endl;
        }
    }

    // Output phase field
    vtkfile << "SCALARS phase_field float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << phi[i][j] << std::endl;
        }
    }

    //output pressure
    vtkfile << "SCALARS pressure float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << p[i][j] << std::endl;
        }
    }

    vtkfile.close();
}

// #include <iomanip>  // For controlling float print precision
// void outputDataCSV(const std::string& baseName, int timestep) {
void outputDataCSV(int timestep) {
    std::ofstream fileDensity, fileVelocity, filePhase, filePressure;
    std::string densityFilename = "density_" + std::to_string(timestep) + ".csv";
    std::string velocityFilename = "velocity_" + std::to_string(timestep) + ".csv";
    std::string phaseFilename = "phase_" + std::to_string(timestep) + ".csv";
    std::string pressureFilename = "pressure_" + std::to_string(timestep) + ".csv";

    fileDensity.open(densityFilename);
    fileVelocity.open(velocityFilename);
    filePhase.open(phaseFilename);
    filePressure.open(pressureFilename);

    // Set precision for float values
    // fileDensity << std::fixed << std::setprecision(8);
    // fileVelocity << std::fixed << std::setprecision(8);
    // filePhase << std::fixed << std::setprecision(8);
    // filePressure << std::fixed << std::setprecision(8);

    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            fileDensity << rho[i][j];
            fileVelocity << u[i][j][0] << "," << u[i][j][1];
            filePhase << phi[i][j];
            filePressure << p[i][j];
            if (i < Lx - 1) {
                fileDensity << ",";
                fileVelocity << ",";
                filePhase << ",";
                filePressure << ",";
            }
        }
        fileDensity << "\n";
        fileVelocity << "\n";
        filePhase << "\n";
        filePressure << "\n";
    }

    fileDensity.close();
    fileVelocity.close();
    filePhase.close();
    filePressure.close();
}

// Function to initialize the simulation for ellipsoidal droplet
const double a2 = 1.5;
const double b2 = 1./a2;

void initialize() {
    // Initialize distribution functions, density, velocity, phase field function, pressure at infinity, etc.
    // Define the center of the bubble/droplet
    int x0 = Lx / 2;
    int y0 = Ly / 2;
    // Define the radius of the droplet/bubble and the interface width
    double r =  radius * dx ; // Lx / 8. * dx; // 16 lattice units
    
    // In order to get f_eq, we need to calculate the macroscopic variables rho(with rho1, rho2 at different nodes), u, p
    for (int i=0; i<Lx ; i++){
        for (int j=0; j<Ly ; j++){
            // \phi(x, y, 0)=\tanh \frac{y-y_0-0.2 L \cos \left(\frac{-2 \pi x}{L}\right)}{W_0}
            double phi_local = tanh((j-y0-0.2*Lx*cos(-2*M_PI*i/Lx))/ch_width_ope);
            phi[i][j] = phi_local;  // Local phase field
            u[i][j][0] = 0.0; // static flow field
            u[i][j][1] = 0.0; // static flow field
            double rho_local = rho1 * (0.5 + 0.5*phi_local) + rho2 * (0.5 - 0.5*phi_local);
            rho[i][j] = rho_local;  // Loc-al density
            rho_mdt[i][j] = rho_local; // At zero time step, the value of previous step is the one of current step
            
            // // Calculate pressure from equation of state
            double p_local = rho_local*((1+phi_local)*0.5*c1*c1 + (1-phi_local)*0.5*c2*c2)- (1+phi_local)*0.5 * p1_inf - (1-phi_local)/2.*p2_inf;          
            p[i][j] = p_local;
            p_mdt[i][j] = p_local;
            F[i][j][0] = 0.0; // External force in lattice units on x direction
            F[i][j][1] = - rho_local * a_g; // External force in lattice units on y direction
            //F[i][j][1] = rho[i][j] * 1.0e-3/(3.*347.*347.); // External force in lattice units on y direction
        }
    }
    calEquilibrium();
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            for (int k=0 ; k<Q ; k++){
                f[i][j][k] = f_eq[i][j][k];
                g[i][j][k] = f[i][j][k] * phi[i][j];
            }
        }
    }
}

// Function to run the simulation
void runSimulation() {
    double r3 = pow((sqrt(a2) + sqrt(b2)) / 2.* radius, 3);
    std::cout << "sigma = " << sigma << std::endl;
    std::cout << "radius = " << radius << std::endl;
    std::cout << "p1_inf = " << p1_inf << std::endl;
    std::cout << "p2_inf = " << p2_inf << std::endl;
    std::cout << "T_theo = " << 2*3.1415 * sqrt((rho1+rho2)*r3/6./sigma) << std::endl; //T_{\text {theo }}=2 \pi \sqrt{\frac{\left(\rho_1+\rho_2\right) r^3}{6 \sigma}}
    initialize();
    //outputVTK("lbm_output_", 0);
    outputDataCSV(0);
    int i;
    for (i = 1; i < numSteps+1; i++) {
        // std::cout << "Step " << i << std::endl;
        force();
        collide();
        collide_surface();
        recolor();
        stream();
        calMacroscopic();
        calPhaseField();
        calEquilibrium();
        
        
        //Output or visualization code here
        if (i % interval == 0) {  // Output every 100 steps
            //outputVTK("lbm_output_", i);
            std::cout << "Step " << i << std::endl;
            outputDataCSV(i);
        }
        // outputVTK("lbm_output_", i);
        // outputDataCSV(i);
        }
}


int main() {
    runSimulation();
    return 0;
}

