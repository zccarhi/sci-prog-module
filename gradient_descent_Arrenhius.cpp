//
//  main.cpp
//  CW2
//
//  Created by Rafiq Hilali on 15/11/2016.
//  Copyright © 2016 Rafiq Hilali. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

float R = 8.3144598; //Real Gas Constant named as global variable.


double LeastSquaresRegression(double y_values[], double x_values[], double EA_tmp, double A_tmp) {
    //this function calculates the chisquared value of a given Arrenhius plot
    int n = 0;
    double chisquared = 0;
    
    for (n = 0; n < 8; n++){
        //iterate through the x and y values and caluclate chiqsuared.
        float ybar = log(A_tmp) - (EA_tmp/R)*x_values[n];
        chisquared += pow(y_values[n] - ybar, 2) / 8;
    }
    
    return chisquared;
}

double GradientDescentA(double y_values[], double x_values[], double EA_tmp, double A_tmp, double A_step){
    //this function will use the gradient descent algorithm to determine the next iteration value of A.
    //To do this the chi_squared equation is differentiated with respect to A.
    // The gradient is then multiplied by a constant step to determine a new value of A
    int n = 0;
    double d_by_dA = 0;
    
    for (n = 0; n < 8; n++){
        d_by_dA += (-2/(8*A_tmp))*(y_values[n]-(log(A_tmp) - (EA_tmp/R)*x_values[n]));
    }
    
    double A_new = A_tmp - (d_by_dA * A_step);
    
    return A_new;
    
    
}

double GradientDescentEA(double y_values[], double x_values[], double EA_tmp, double A_tmp, double EA_step){
    //this function will use the gradient descent algorithm to determine the next iteration value of EA.
    //To do this the chi_squared equation is differentiated with respect to EA.
    // The gradient is then multiplied by a constant step to determine a new value of EA
    int n = 0;
    double d_by_dEA = 0;
    
    for (n = 0; n < 8; n++){
        d_by_dEA += ((2*x_values[n])/(8*R))*(y_values[n]-(log(A_tmp) - (EA_tmp/R)*x_values[n]));
    }
    
    double EA_new = EA_tmp - (d_by_dEA * EA_step);
    
    return EA_new;

}



int main() {
    //First, instialise the data for T and k
    
    float T[8] = {700.0, 730.0, 760.0, 790.0, 810.0, 840.0, 910.0, 1000.0};
    float k[8] = {0.011, 0.035, 0.105, 0.343, 0.789, 2.17, 20.0, 145.0};
    int i = 0;
    double y[8] = {};
    double x[8] = {};
    double EA = 15000; //Starting EA
    double A = 0.9e+12; //Starting A
    // Gradient times by step is the amount A/EA changes with each iteration.
    double A_step = 5e18;
    double EA_step = 500000;
    // Two values of chi^2 (measure of least squares regression) to check for convergence
    double chi = 0;
    double chi_new = 0;
    
    for (i=0; i < 8; i++){
        // x and y values are created in accordance with Arrenhius Plot
        y[i] = log(k[i]);
        x[i] = 1/T[i];
    };
    
    for(i=0; i<10000; i++){
        // iterates using the Gradient Descent Algorithm to the minimum chi_squared value.
        // Maximum loops is set to 10000 so infinite loop does not occur if break criteria is not met
        // When the difference between chi and chi_new is less than a certain threshold, the calculation has converged and the loop breaks
        chi_new = LeastSquaresRegression(y, x, EA, A);
        A = GradientDescentA(y, x, EA, A, A_step);
        EA = GradientDescentEA(y, x, EA, A, EA_step);
        
        if (abs(chi - chi_new) < 0.0000000001){
            float EA_kJ = EA / 1000;
            cout << "The calculation completed in " << i << " iterations." << endl;
            cout << "The chi squared value = " << chi_new << endl;
            cout <<  "The Pre exponential factor = " << A << " mols-1" << endl;
            cout <<  "The Activation Energy = " << EA_kJ << " kJ" << endl;
            break;
        }
        
        chi = chi_new;
    }
    
    return 0;

}
    
// The Gradient Descent Algorithm was used to find the least squares regression line of the sample.
// This was defined as two functions that found the derivative of the least squares regressions equation with respect to EA and A respectively.
// This was used to direct the 'step' in each iteration to a minimum value of chi squared.
// A constant step was multiplied by the derivative and subtracted from the EA/A value to gain a new EA/A value.
// Another function was defined to find the chi squared value (least squares regression).
// A loop is then set up where the Gradient Descent Algorithm is used to determine new values of EA and A.
// The chi squared value is then calculated and compared to the previous chi_squared value.
// When the difference between them is below a certain threshold the loop breaks and the calculation is complete.
// The final EA and A value are then printed


    
    
    
    

