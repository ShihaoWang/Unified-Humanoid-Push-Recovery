#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <dlib/matrix.h>

#include "Unified_Header.h"
#include "snfilewrapper.hh"

using namespace std;

double CubicSpline_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s)
{
	double a, b, c, d, y;
	a = CubicSpline_Coeff[0];			b = CubicSpline_Coeff[1];
	c = CubicSpline_Coeff[2];			d = CubicSpline_Coeff[3];
	y = a * s * s * s + b * s * s + c * s + d;
	return y;
}
double CubicSpline_1stOrder_Evaluation_fn(const std::vector<double> &CubicSpline_Coeff, double s, double T)
{ 	// This function is used to output the first order derivative value of a given spline at s position
	double a, b, c, d, y;
	a = CubicSpline_Coeff[0];			b = CubicSpline_Coeff[1];
	c = CubicSpline_Coeff[2];			d = CubicSpline_Coeff[3];
	y = (3 * a * s * s + 2 * b * s + c)/T;
	return y;
}

std::vector<double> CubicSpline_PosVelAcc4(double T, double a, double b, double c, double d, double s)
{
	//# This function is used to calcualte the position, velocity and acceleration given the spline coefficients
    //# Here T is the duration constant, s is the path variable
	std::vector<double> PVA(3);
	double Pos = a*s*s*s + b*s*s + c*s + d;
    double Vel = (3*a*s*s + 2*b*s + c)/T;
    double Acc = (6*a*s+2*b)/(T*T);
	PVA[0] = Pos; 			PVA[1] = Vel; 			PVA[2] = Acc;
	return PVA;
}
std::vector<double> CubicSpline_PosVelAcc8(double T, double x_a, double x_b, double x_c, double x_d, double xdot_a, double xdot_b, double xdot_c, double xdot_d, double s)
{
	std::vector<double> PVAVP(4);
	double Pos = x_a*s*s*s + x_b*s*s + x_c*s + x_d;
	double Vel =  xdot_a*s*s*s + xdot_b*s*s + xdot_c*s + xdot_d;
	double Acc = (3*xdot_a*s*s + 2*xdot_b*s + xdot_c)/T;
	double VelfromPos = (3*x_a*s*s + 2*x_b*s + x_c)/T;
	PVAVP[0] = Pos;				PVAVP[1] = Vel;				PVAVP[2] = Acc;				PVAVP[3] = VelfromPos;
	return PVAVP;
}
std::vector<double> CubicSpline_Coeff_fn(double T, double x_init, double x_end, double xdot_init, double xdot_end)
{	// This function is used to calcualte the coefficients for the cubic spline
	// The cubic spline is expressed to be : y(s) = a*s^3 + b*s^2 + c*s + d
	std::vector<double> CubicSpline_Coeff_vec(4);
	double a, b, c, d;
	a = 2*x_init - 2*x_end + T*xdot_end + T*xdot_init;
  	b = 3*x_end - 3*x_init - T*xdot_end - 2*T*xdot_init;
    c = T*xdot_init;
    d = x_init;
	CubicSpline_Coeff_vec[0] = a;
	CubicSpline_Coeff_vec[1] = b;
	CubicSpline_Coeff_vec[2] = c;
	CubicSpline_Coeff_vec[3] = d;
    return CubicSpline_Coeff_vec;
}
