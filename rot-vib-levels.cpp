#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
// Instalise constants as global variables

double h = 6.62607004E-34;
double k = 1.838064852E-23; 
double c = 2.99792458E8;

float RelInt(int J, int B, int T){
    // This function returns the relative intensity of a give transition
    float RelInt = (2J+1) * exp((- h * c * B * J * (J + 1))/(k * T));
    return RelInt;
}

void RoVibLevels(float B, float alpha, float omega, float omegachi) {
    // This function returns as a list the position and relative intensity of
    // the P and R branches for the rotovibrational spectrum
    
    // deltqvib gives the wavenumber of the v=0 --> v=1 transition
    float deltavib = (1.5 * omega - pow(1.5, 2) * omegachi) - (0.5 * omega - pow(0.5, 2) * omegachi);
    
    // The use of seperate B values for different vibrational modes means the
    // centrifugal distortion correction is not required
    float B0 = B - alpha * 0.5;
    float B1 = B - alpha * 1.5;

    // The following loop prints the positions and intensities of the P branch
    cout << endl << "P Branch" << endl << endl;
    for (int J = 0; J < 10; J++) {
        cout << "J " << J + 1 << "--->" << J << endl;
        float wavenumber = deltavib + (B1 * (J - 1) * J) - (B0 * J * (J + 1));
        cout << wavenumber << " cm-1" << endl;
        cout << "Relative Intensity: " << RelInt(J + 1, B1, 298) << endl;
    }
    
    // the following loop prints the positions and intensities of the R branch
    cout << endl <<"R branch" << endl << endl;
    for (int J = 0; J < 10; J++) {
        cout << "J " << J << "--->" << J + 1 << endl;
        float wavenumber = deltavib + (B1 * (J + 1) * J) - (B0 * J * (J - 1));
        cout << wavenumber << " cm-1" << endl;
        cout << "Relative Intensity: " << RelInt(J, B1, 298) << endl;
    }
}
        
int main(){
    // Calls on RoVibLevels for HCl
    cout << "Rotational-Vibrational Spectrum for H 35Cl" << endl;
    RoVibLevels(10.59341, 0.307, 2990.946, 52.8186);
    cout << endl << endl;
    // Calls on RovibLevels for DCl
    cout << "Rotational- Vibrational Spectrum for D 35Cl" << endl;
    RoVibLevels(5.44879, 0.113291, 2145.16, 27.1825);

    return 0;
}

// Two functions are defined at the start of the program. One uses the
// Boltzmann distribution to find the relative occupancy of the lower
// vibrational-rotational level and thus calculate the relative intensity of
// the transition. The second function uses the rotational constant and quantum
// numbers to calculate the wavenumber of the transitions of the R and P branch
// and prints them, along with their relative intensities. In the main()
// function, the RovibLevels function is called on for HCl and DCl.
