
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

// Initalise universal constants as global variables outside of main function

double ele_charge = 1.6021766E-19;
double eps_0 = 8.85418782E-12;
double h = 6.62607004E-34;
double c = 2.99792458E8;
// mass electron/neutron/proton
double m_e = 9.109390E-31;
double m_n = 1.674929E-27;
double m_p = 1.672623E-27;

double ReducedMass(int protons, int neutrons, int electrons) {
    // This function returns the reduced mass for an atom.
    double reducedmass;
    // if and else statments to to handle positronium.
    if(neutrons == 0 and protons == 0){
        reducedmass = (m_e * m_e) / (2 * m_e);
    }
    else{
        reducedmass = ((protons * m_p + neutrons * m_n) * electrons* m_e)/((protons*m_p)+(neutrons*m_n)+(electrons*m_e));
    }
    
    return reducedmass;
    
}

double RydbergCalc(double reducedmass) {
    // This function returns the Rydberg Constant from fundamental constants
    
    double R_calc = (pow(ele_charge, 4) * reducedmass)/(8 * (pow(eps_0, 2) * pow(h, 3) * c));
    
    return R_calc;
}

double RydbergEquation(double rydberg, int n_1, int n_2) {
    // This function returns the wavenumber of the transition from n 1 to n 2.
    double sigma = -rydberg * ((1/pow(n_2, 2)) - (1/pow(n_1, 2)));
    
    return sigma;
}

double LandeGFactor(int j, int s, int l){
    // Calculates the Lande G Factor for given values of j, s and l.
    double LGF = (3 * j * (j + 1) + s * (s + 1) - l * (l + 1)) / (2 * j * (j + 1));
    
    return LGF;
}

int main() {
    int protons;
    int neutrons;
    int i;
    int j;
    int s;
    int l;
    double R_calc;
    
    // Ask the user to input the number of protons and electron in their atom
    cout << "For a positronium state the number of protons and neutrons as 0" << endl;
    cout << "How many protons does your hydrogenic atom have?" << endl;
    cin >> protons;
    cout << endl << "How many neutrons does your hydrogenic atom have?" << endl;
    cin >> neutrons;
    
    // Calculate R_calc for the atom of interest.
    R_calc = RydbergCalc(ReducedMass(protons, neutrons, 1));
    
    // The first 5 lines of the Balmer series are printed. For other series the value of n_1 can be altered.
    
    for(i = 0; i < 5; i++){
        cout << "The " << (i + 1) << " line of the Balmer Series is: " << (RydbergEquation(R_calc, 2, i + 3) / 100) << "cm-1" << endl;
    }
    
    // Ask the user for variables for Lande-G Factor Calculation
    cout << "Lande-G Factor Calculation:" << endl;
    cout << "Enter your value of j:" << endl;
    cin >> j;
    cout << "Enter your value of s:" << endl;
    cin >> s;
    cout << "Enter your value of l:" << endl;
    cin >> l;
    // print value of g_j
    cout << "Lande-G Factor: " << LandeGFactor(j, s, l) << endl;;
    
    return 0;

    
}
