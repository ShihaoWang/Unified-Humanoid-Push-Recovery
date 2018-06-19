#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <dlib/matrix.h>

#include "Unified_Header.h"
#include "snfilewrapper.hh"
extern double Inf;
extern double PI;
extern double pi;
using namespace std;

dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;
    double rIy = Robot_StateNDot_i.rIy;
    double theta = Robot_StateNDot_i.theta;
    double q1 = Robot_StateNDot_i.q1;
    double q2 = Robot_StateNDot_i.q2;
    double q3 = Robot_StateNDot_i.q3;
    double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	dlib::matrix<double>  T;
	T = dlib::zeros_matrix<double>(13,13);
	T(0,0) = 2.71E2/5.0;
	T(0,2) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(6.3E1/8.0E1)-cos(q9+q10+theta)*(6.3E1/8.0E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(9.1E1/8.0E1)-cos(q9+theta)*(9.1E1/8.0E1)+cos(theta)*1.3585E1-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,3) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,4) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,5) = sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(-1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,6) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,7) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,8) = sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(-1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(0,9) = cos(q7+q8+theta)*(-6.3E1/8.0E1)-cos(q7+theta)*(9.1E1/8.0E1);
	T(0,10) = cos(q7+q8+theta)*(-6.3E1/8.0E1);
	T(0,11) = cos(q9+q10+theta)*(-6.3E1/8.0E1)-cos(q9+theta)*(9.1E1/8.0E1);
	T(0,12) = cos(q9+q10+theta)*(-6.3E1/8.0E1);
	T(1,1) = 2.71E2/5.0;
	T(1,2) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(6.3E1/8.0E1)+sin(q9+q10+theta)*(6.3E1/8.0E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(9.1E1/8.0E1)+sin(q9+theta)*(9.1E1/8.0E1)-sin(theta)*1.3585E1+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,3) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,4) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,6) = sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,7) = sin(q4+q5+theta)*(9.1E1/1.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,8) = sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,5) = sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(1,9) = sin(q7+q8+theta)*(6.3E1/8.0E1)+sin(q7+theta)*(9.1E1/8.0E1);
	T(1,10) = sin(q7+q8+theta)*(6.3E1/8.0E1);
	T(1,11) = sin(q9+q10+theta)*(6.3E1/8.0E1)+sin(q9+theta)*(9.1E1/8.0E1);
	T(1,12) = sin(q9+q10+theta)*(6.3E1/8.0E1);
	T(2,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q4+q5+theta)*(9.1E1/1.0E2)-cos(q7+q8+theta)*(6.3E1/8.0E1)-cos(q9+q10+theta)*(6.3E1/8.0E1)-cos(q1+theta)*(3.9E1/2.0E1)-cos(q4+theta)*(3.9E1/2.0E1)-cos(q7+theta)*(9.1E1/8.0E1)-cos(q9+theta)*(9.1E1/8.0E1)+cos(theta)*1.3585E1-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(2,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q7+q8+theta)*(6.3E1/8.0E1)+sin(q9+q10+theta)*(6.3E1/8.0E1)+sin(q1+theta)*(3.9E1/2.0E1)+sin(q4+theta)*(3.9E1/2.0E1)+sin(q7+theta)*(9.1E1/8.0E1)+sin(q9+theta)*(9.1E1/8.0E1)-sin(theta)*1.3585E1+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(2,2) = cos(q7+q8)*(-6.93E2/8.0E2)-cos(q9+q10)*(6.93E2/8.0E2)+cos(q2)*5.915E-1+cos(q5)*5.915E-1-cos(q7)*1.25125+cos(q8)*(6.3E1/1.6E2)-cos(q9)*1.25125+cos(q10)*(6.3E1/1.6E2)+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/2.5E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/2.5E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.759;
	T(2,3) = cos(q2)*5.915E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(2,4) = cos(q2)*2.9575E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(2,5) = sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(2,6) = cos(q5)*5.915E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(2,7) = cos(q5)*2.9575E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(2,8) = sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(2,9) = cos(q7+q8)*(-4.33125E-1)-cos(q7)*6.25625E-1+cos(q8)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(2,10) = cos(q7+q8)*(-4.33125E-1)+cos(q8)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(2,11) = cos(q9+q10)*(-4.33125E-1)-cos(q9)*6.25625E-1+cos(q10)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(2,12) = cos(q9+q10)*(-4.33125E-1)+cos(q10)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(3,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-cos(q1+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(3,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sin(q1+theta)*(3.9E1/2.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(3,2) = cos(q2)*5.915E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(3,3) = cos(q2)*5.915E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(3,4) = cos(q2)*2.9575E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(3,5) = sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(4,0) = cos(q1+q2+theta)*(-9.1E1/1.0E2)-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(4,1) = sin(q1+q2+theta)*(9.1E1/1.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(4,2) = cos(q2)*2.9575E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(4,3) = cos(q2)*2.9575E-1+sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(4,4) = sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(4,5) = sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(5,0) = sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(-1.0/5.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(5,1) = sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(5,2) = sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(5,3) = sqrt(2.0)*cos(q2+q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(5,4) = sqrt(4.1E1)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q3+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(5,5) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+1.0/4.0E1;
	T(6,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-cos(q4+theta)*(3.9E1/2.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(6,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+sin(q4+theta)*(3.9E1/2.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(6,2) = cos(q5)*5.915E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(6,6) = cos(q5)*5.915E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+8.418333333333333E-1;
	T(6,7) = cos(q5)*2.9575E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(6,8) = sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(7,0) = cos(q4+q5+theta)*(-9.1E1/1.0E2)-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(7,1) = sin(q4+q5+theta)*(9.1E1/1.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(7,2) = cos(q5)*2.9575E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(7,6) = cos(q5)*2.9575E-1+sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(7,7) = sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/2.5E2)+2.785E-1;
	T(7,8) = sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(8,0) = sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(-1.0/5.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1);
	T(8,1) = sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*(2.0/2.5E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1);
	T(8,2) = sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(8,6) = sqrt(2.0)*cos(q5+q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(8,7) = sqrt(4.1E1)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(2.0)*cos(q6+3.141592653589793*(1.0/4.0))*(1.3E1/5.0E2)+1.0/4.0E1;
	T(8,8) = sqrt(4.1E1)*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*sin(8.960553845713439E-1)*(1.0/5.0E2)+1.0/4.0E1;
	T(9,0) = cos(q7+q8+theta)*(-6.3E1/8.0E1)-cos(q7+theta)*(9.1E1/8.0E1);
	T(9,1) = sin(q7+q8+theta)*(6.3E1/8.0E1)+sin(q7+theta)*(9.1E1/8.0E1);
	T(9,2) = cos(q7+q8)*(-4.33125E-1)-cos(q7)*6.25625E-1+cos(q8)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(9,9) = cos(q8)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(9,10) = cos(q8)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(10,0) = cos(q7+q8+theta)*(-6.3E1/8.0E1);
	T(10,1) = sin(q7+q8+theta)*(6.3E1/8.0E1);
	T(10,2) = cos(q7+q8)*(-4.33125E-1)+cos(q8)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(10,9) = cos(q8)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(10,10) = 1.954166666666667E-1;
	T(11,0) = cos(q9+q10+theta)*(-6.3E1/8.0E1)-cos(q9+theta)*(9.1E1/8.0E1);
	T(11,1) = sin(q9+q10+theta)*(6.3E1/8.0E1)+sin(q9+theta)*(9.1E1/8.0E1);
	T(11,2) = cos(q9+q10)*(-4.33125E-1)-cos(q9)*6.25625E-1+cos(q10)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(11,11) = cos(q10)*(6.3E1/1.6E2)+4.824166666666667E-1;
	T(11,12) = cos(q10)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(12,0) = cos(q9+q10+theta)*(-6.3E1/8.0E1);
	T(12,1) = sin(q9+q10+theta)*(6.3E1/8.0E1);
	T(12,2) = cos(q9+q10)*(-4.33125E-1)+cos(q10)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(12,11) = cos(q10)*(6.3E1/3.2E2)+1.954166666666667E-1;
	T(12,12) = 1.954166666666667E-1;
	return T;
}
dlib::matrix<double> B_q_fn()
{
	dlib::matrix<double> B_q;
	B_q = dlib::zeros_matrix<double>(13,10);
	B_q(3,0) = 1.0;
	B_q(4,1) = 1.0;
	B_q(5,2) = 1.0;
	B_q(6,3) = 1.0;
	B_q(7,4) = 1.0;
	B_q(8,5) = 1.0;
	B_q(9,6) = 1.0;
	B_q(10,7) = 1.0;
	B_q(11,8) = 1.0;
	B_q(12,9) = 1.0;
	return B_q;
}
dlib::matrix<double> C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;
	double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	double rIxdot = Robot_StateNDot_i.rIxdot;
	double rIydot = Robot_StateNDot_i.rIydot;
	double thetadot = Robot_StateNDot_i.thetadot;
	double q1dot = Robot_StateNDot_i.q1dot;
	double q2dot = Robot_StateNDot_i.q2dot;
	double q3dot = Robot_StateNDot_i.q3dot;
	double q4dot = Robot_StateNDot_i.q4dot;
	double q5dot = Robot_StateNDot_i.q5dot;
	double q6dot = Robot_StateNDot_i.q6dot;
	double q7dot = Robot_StateNDot_i.q7dot;
	double q8dot = Robot_StateNDot_i.q8dot;
	double q9dot = Robot_StateNDot_i.q9dot;
	double q10dot = Robot_StateNDot_i.q10dot;

	dlib::matrix<double> T;
	T = dlib::zeros_matrix<double>(13,1);
	T(0,0) = (thetadot*thetadot)*sin(theta)*(-1.3585E1)+(q10dot*q10dot)*sin(q9+q10+theta)*(6.3E1/8.0E1)+(q1dot*q1dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*sin(q7+q8+theta)*(6.3E1/8.0E1)+(q8dot*q8dot)*sin(q7+q8+theta)*(6.3E1/8.0E1)+(q9dot*q9dot)*sin(q9+q10+theta)*(6.3E1/8.0E1)+(thetadot*thetadot)*sin(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*sin(q7+q8+theta)*(6.3E1/8.0E1)+(thetadot*thetadot)*sin(q9+q10+theta)*(6.3E1/8.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*sin(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*sin(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*sin(q7+theta)*(9.1E1/8.0E1)+(q9dot*q9dot)*sin(q9+theta)*(9.1E1/8.0E1)+(thetadot*thetadot)*sin(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*sin(q7+theta)*(9.1E1/8.0E1)+(thetadot*thetadot)*sin(q9+theta)*(9.1E1/8.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*sin(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*sin(q7+q8+theta)*(6.3E1/4.0E1)+q10dot*thetadot*sin(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*sin(q7+q8+theta)*(6.3E1/4.0E1)+q8dot*thetadot*sin(q7+q8+theta)*(6.3E1/4.0E1)+q9dot*thetadot*sin(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*sin(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*sin(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*sin(q7+theta)*(9.1E1/4.0E1)+q9dot*thetadot*sin(q9+theta)*(9.1E1/4.0E1)+sqrt(4.1E1)*q1dot*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1);
	T(1,0) = (thetadot*thetadot)*cos(theta)*(-1.3585E1)+(q10dot*q10dot)*cos(q9+q10+theta)*(6.3E1/8.0E1)+(q1dot*q1dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q2dot*q2dot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(q4dot*q4dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q5dot*q5dot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(q7dot*q7dot)*cos(q7+q8+theta)*(6.3E1/8.0E1)+(q8dot*q8dot)*cos(q7+q8+theta)*(6.3E1/8.0E1)+(q9dot*q9dot)*cos(q9+q10+theta)*(6.3E1/8.0E1)+(thetadot*thetadot)*cos(q1+q2+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q4+q5+theta)*(9.1E1/1.0E2)+(thetadot*thetadot)*cos(q7+q8+theta)*(6.3E1/8.0E1)+(thetadot*thetadot)*cos(q9+q10+theta)*(6.3E1/8.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(2.0/2.5E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(2.0/2.5E1)-(q1dot*q1dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q2dot*q2dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q3dot*q3dot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(q4dot*q4dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q5dot*q5dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(q6dot*q6dot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q1+q2+q3+theta)*(2.0/2.5E1)-(thetadot*thetadot)*sin(q4+q5+q6+theta)*(2.0/2.5E1)+(q1dot*q1dot)*cos(q1+theta)*(3.9E1/2.0E1)+(q4dot*q4dot)*cos(q4+theta)*(3.9E1/2.0E1)+(q7dot*q7dot)*cos(q7+theta)*(9.1E1/8.0E1)+(q9dot*q9dot)*cos(q9+theta)*(9.1E1/8.0E1)+(thetadot*thetadot)*cos(q1+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q4+theta)*(3.9E1/2.0E1)+(thetadot*thetadot)*cos(q7+theta)*(9.1E1/8.0E1)+(thetadot*thetadot)*cos(q9+theta)*(9.1E1/8.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+q10dot*q9dot*cos(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*q8dot*cos(q7+q8+theta)*(6.3E1/4.0E1)+q10dot*thetadot*cos(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(9.1E1/5.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(9.1E1/5.0E1)+q7dot*thetadot*cos(q7+q8+theta)*(6.3E1/4.0E1)+q8dot*thetadot*cos(q7+q8+theta)*(6.3E1/4.0E1)+q9dot*thetadot*cos(q9+q10+theta)*(6.3E1/4.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q1dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*q3dot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q4dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*q6dot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q2dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q3dot*thetadot*cos(q1+q2+q3+theta)*(4.0/2.5E1)+q4dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q5dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)+q6dot*thetadot*cos(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*q2dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q1dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*q3dot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*q5dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q4dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*q6dot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q1dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q2dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q3dot*thetadot*sin(q1+q2+q3+theta)*(4.0/2.5E1)-q4dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q5dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)-q6dot*thetadot*sin(q4+q5+q6+theta)*(4.0/2.5E1)+q1dot*thetadot*cos(q1+theta)*(3.9E1/1.0E1)+q4dot*thetadot*cos(q4+theta)*(3.9E1/1.0E1)+q7dot*thetadot*cos(q7+theta)*(9.1E1/4.0E1)+q9dot*thetadot*cos(q9+theta)*(9.1E1/4.0E1)+sqrt(4.1E1)*q1dot*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q1dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q2dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q3dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q4dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q5dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+sqrt(4.1E1)*q6dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.5E1)+5.31702E2;
	T(2,0) = sin(q1+q2+theta)*8.9271+sin(q4+q5+theta)*8.9271+sin(q7+q8+theta)*7.725375+sin(q9+q10+theta)*7.725375+sin(q1+theta)*1.91295E1+sin(q4+theta)*1.91295E1+sin(q7+theta)*1.1158875E1+sin(q9+theta)*1.1158875E1-sin(theta)*1.3326885E2-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q10dot*q10dot)*sin(q10)*(6.3E1/3.2E2)-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7)*6.25625E-1-(q8dot*q8dot)*sin(q8)*(6.3E1/3.2E2)+(q9dot*q9dot)*sin(q9)*6.25625E-1+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)+(q10dot*q10dot)*sin(q9+q10)*4.33125E-1-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*sin(q7+q8)*4.33125E-1+(q8dot*q8dot)*sin(q7+q8)*4.33125E-1+(q9dot*q9dot)*sin(q9+q10)*4.33125E-1-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q10dot*q9dot*sin(q10)*(6.3E1/1.6E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q7dot*q8dot*sin(q8)*(6.3E1/1.6E2)-q10dot*thetadot*sin(q10)*(6.3E1/1.6E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7)*1.25125-q8dot*thetadot*sin(q8)*(6.3E1/1.6E2)+q9dot*thetadot*sin(q9)*1.25125-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q10dot*q9dot*sin(q9+q10)*(6.93E2/8.0E2)-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*q8dot*sin(q7+q8)*(6.93E2/8.0E2)+q10dot*thetadot*sin(q9+q10)*(6.93E2/8.0E2)-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+q7dot*thetadot*sin(q7+q8)*(6.93E2/8.0E2)+q8dot*thetadot*sin(q7+q8)*(6.93E2/8.0E2)+q9dot*thetadot*sin(q9+q10)*(6.93E2/8.0E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(3,0) = sin(q1+q2+theta)*8.9271+sin(q1+theta)*1.91295E1-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q2dot*q2dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1-(q2dot*q2dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*cos(q2+q3)*(1.3E1/5.0E2)-(q2dot*q2dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q3dot*q3dot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2)*5.915E-1-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q2dot*q2dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q3dot*q3dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)-q1dot*q2dot*sin(q2+q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q2dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q1dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(4,0) = sin(q1+q2+theta)*8.9271-(q3dot*q3dot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q2)*2.9575E-1-(q3dot*q3dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2)*2.9575E-1+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q2dot*q3dot*cos(q3)*(1.3E1/2.5E2)-q3dot*thetadot*cos(q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q3)*(1.3E1/2.5E2)-q2dot*q3dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2)*5.915E-1-q3dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q3dot*q3dot)*sin(q3-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q1dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q2dot*q3dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q3dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(5,0) = (q1dot*q1dot)*cos(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*cos(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q3)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q1dot*q1dot)*sin(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*sin(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*1.962E-1+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q3-8.960553845713439E-1)*6.5E-3+q1dot*q2dot*cos(q3)*(1.3E1/2.5E2)+q1dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q2dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q1dot*q2dot*sin(q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q3)*(1.3E1/2.5E2)+q2dot*thetadot*sin(q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q2+q3-8.960553845713439E-1)*6.5E-3+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q2dot*q2dot)*sin(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*q2dot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q2dot*thetadot*sin(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*thetadot*sin(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(6,0) = sin(q4+q5+theta)*8.9271+sin(q4+theta)*1.91295E1-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1-(q5dot*q5dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1-(q5dot*q5dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*cos(q5+q6)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q5+q6)*(1.3E1/5.0E2)-(q6dot*q6dot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5)*5.915E-1-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q5dot*q5dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*(q6dot*q6dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3-q4dot*q5dot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q5dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q4dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(7,0) = sin(q4+q5+theta)*8.9271-(q6dot*q6dot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q5)*2.9575E-1-(q6dot*q6dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5)*2.9575E-1+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q5dot*q6dot*cos(q6)*(1.3E1/2.5E2)-q6dot*thetadot*cos(q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q6)*(1.3E1/2.5E2)-q5dot*q6dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5)*5.915E-1-q6dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-sqrt(4.1E1)*(q6dot*q6dot)*sin(q6-8.960553845713439E-1)*6.5E-3-sqrt(4.1E1)*q4dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q5dot*q6dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)-sqrt(4.1E1)*q6dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(8,0) = (q4dot*q4dot)*cos(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*cos(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q6)*(1.3E1/5.0E2)+sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/4.0))*7.848E-1+(q4dot*q4dot)*sin(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*sin(q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*1.962E-1+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)+(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q6-8.960553845713439E-1)*6.5E-3+q4dot*q5dot*cos(q6)*(1.3E1/2.5E2)+q4dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q5dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q4dot*q5dot*sin(q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q6)*(1.3E1/2.5E2)+q5dot*thetadot*sin(q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*sin(q5+q6-8.960553845713439E-1)*6.5E-3+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q5dot*q5dot)*sin(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q5dot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q5dot*thetadot*sin(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*thetadot*sin(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3);
	T(9,0) = sin(q7+q8+theta)*7.725375+sin(q7+theta)*1.1158875E1-(q8dot*q8dot)*sin(q8)*(6.3E1/3.2E2)-(thetadot*thetadot)*sin(q7)*6.25625E-1-(thetadot*thetadot)*sin(q7+q8)*4.33125E-1-q7dot*q8dot*sin(q8)*(6.3E1/1.6E2)-q8dot*thetadot*sin(q8)*(6.3E1/1.6E2);
	T(10,0) = sin(q7+q8+theta)*7.725375+(q7dot*q7dot)*sin(q8)*(6.3E1/3.2E2)+(thetadot*thetadot)*sin(q8)*(6.3E1/3.2E2)-(thetadot*thetadot)*sin(q7+q8)*4.33125E-1+q7dot*thetadot*sin(q8)*(6.3E1/1.6E2);
	T(11,0) = sin(q9+q10+theta)*7.725375+sin(q9+theta)*1.1158875E1-(q10dot*q10dot)*sin(q10)*(6.3E1/3.2E2)-(thetadot*thetadot)*sin(q9)*6.25625E-1-(thetadot*thetadot)*sin(q9+q10)*4.33125E-1-q10dot*q9dot*sin(q10)*(6.3E1/1.6E2)-q10dot*thetadot*sin(q10)*(6.3E1/1.6E2);
	T(12,0) = sin(q9+q10+theta)*7.725375+(q9dot*q9dot)*sin(q10)*(6.3E1/3.2E2)+(thetadot*thetadot)*sin(q10)*(6.3E1/3.2E2)-(thetadot*thetadot)*sin(q9+q10)*4.33125E-1+q9dot*thetadot*sin(q10)*(6.3E1/1.6E2);
	return T;
}
dlib::matrix<double> Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;
	double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	dlib::matrix<double> T;
	T = dlib::zeros_matrix<double>(12,13);
	T(0,0) = 1.0;
    T(0,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(0,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(0,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(0,5) = sqrt(4.1E1)*sin(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
    T(1,1) = 1.0;
    T(1,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(1,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(1,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(1,5) = sqrt(4.1E1)*cos(q1+q2+q3+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
    T(2,0) = 1.0;
    T(2,2) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(2,3) = sin(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(2,4) = sin(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(2,5) = sqrt(2.0)*sin(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
    T(3,1) = 1.0;
    T(3,2) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(3,3) = cos(q1+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(3,4) = cos(q1+q2+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(1.0/1.0E1);
    T(3,5) = sqrt(2.0)*cos(q1+q2+q3+theta-3.141592653589793*(5.0/4.0))*(-1.0/1.0E1);
    T(4,0) = 1.0;
    T(4,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(4,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(4,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(4,8) = sqrt(4.1E1)*sin(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
    T(5,1) = 1.0;
    T(5,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(5,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(5,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
    T(5,8) = sqrt(4.1E1)*cos(q4+q5+q6+theta+3.141592653589793*(1.0/2.0)-8.960553845713439E-1)*(-1.0/4.0E1);
    T(6,0) = 1.0;
    T(6,2) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(6,6) = sin(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(6,7) = sin(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(6,8) = sqrt(2.0)*sin(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
    T(7,1) = 1.0;
    T(7,2) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(7,6) = cos(q4+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(7,7) = cos(q4+q5+theta+3.141592653589793*(1.0/2.0))*(-1.3E1/4.0E1)-sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(1.0/1.0E1);
    T(7,8) = sqrt(2.0)*cos(q4+q5+q6+theta+3.141592653589793*(3.0/4.0))*(-1.0/1.0E1);
    T(8,0) = 1.0;
    T(8,2) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
    T(8,9) = sin(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
    T(8,10) = sin(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
    T(9,1) = 1.0;
    T(9,2) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
    T(9,9) = cos(q7+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
    T(9,10) = cos(q7+q8+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
    T(10,0) = 1.0;
    T(10,2) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-sin(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
    T(10,11) = sin(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
    T(10,12) = sin(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
    T(11,1) = 1.0;
    T(11,2) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1)-cos(theta-3.141592653589793*(1.0/2.0))*(1.1E1/2.0E1);
    T(11,11) = cos(q9+theta+3.141592653589793*(1.0/2.0))*(-1.0/4.0)-cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(9.0/2.0E1);
    T(11,12) = cos(q9+q10+theta+3.141592653589793*(1.0/2.0))*(-9.0/2.0E1);
	return T;
}
dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i)
{
	double rIx = Robot_StateNDot_i.rIx;
	double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	double rIxdot = Robot_StateNDot_i.rIx;
	double rIydot = Robot_StateNDot_i.rIy;
	double thetadot = Robot_StateNDot_i.thetadot;
	double q1dot = Robot_StateNDot_i.q1dot;
	double q2dot = Robot_StateNDot_i.q2dot;
	double q3dot = Robot_StateNDot_i.q3dot;
	double q4dot = Robot_StateNDot_i.q4dot;
	double q5dot = Robot_StateNDot_i.q5dot;
	double q6dot = Robot_StateNDot_i.q6dot;
	double q7dot = Robot_StateNDot_i.q7dot;
	double q8dot = Robot_StateNDot_i.q8dot;
	double q9dot = Robot_StateNDot_i.q9dot;
	double q10dot = Robot_StateNDot_i.q10dot;

	dlib::matrix<double> T;
	T = dlib::zeros_matrix<double>(12,1);

	T(0,0) = (q1dot*q1dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*sin(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q1dot*q2dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*q3dot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q3dot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1);
	T(1,0) = (q1dot*q1dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q2dot*q2dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q3dot*q3dot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/4.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q1dot*q2dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*q3dot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q1dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q2dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q3dot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/2.0E1);
	T(2,0) = (q1dot*q1dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*sin(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q1+theta)*(1.3E1/4.0E1)+q1dot*q2dot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*sin(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q2dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*sin(q1+theta)*(1.3E1/2.0E1);
	T(3,0) = (q1dot*q1dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q2dot*q2dot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+q2+theta)*(1.3E1/4.0E1)+(q1dot*q1dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q2dot*q2dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(q3dot*q3dot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q1+q2+q3+theta)*(1.0/1.0E1)-(q1dot*q1dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(q2dot*q2dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(q3dot*q3dot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)-(thetadot*thetadot)*sin(q1+q2+q3+theta)*(1.0/1.0E1)+(q1dot*q1dot)*cos(q1+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q1+theta)*(1.3E1/4.0E1)+q1dot*q2dot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q2dot*thetadot*cos(q1+q2+theta)*(1.3E1/2.0E1)+q1dot*q2dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*q3dot*cos(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q2dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)+q3dot*thetadot*cos(q1+q2+q3+theta)*(1.0/5.0)-q1dot*q2dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q1dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q2dot*q3dot*sin(q1+q2+q3+theta)*(1.0/5.0)-q1dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)-q2dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)-q3dot*thetadot*sin(q1+q2+q3+theta)*(1.0/5.0)+q1dot*thetadot*cos(q1+theta)*(1.3E1/2.0E1);
	T(4,0) = (q4dot*q4dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*sin(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q4dot*q5dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*q6dot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q6dot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1);
	T(5,0) = (q4dot*q4dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+theta)*(1.3E1/4.0E1)+sqrt(4.1E1)*(q4dot*q4dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q5dot*q5dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(q6dot*q6dot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+sqrt(4.1E1)*(thetadot*thetadot)*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/4.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+theta)*(1.3E1/2.0E1)+sqrt(4.1E1)*q4dot*q5dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*q6dot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q4dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q5dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1)+sqrt(4.1E1)*q6dot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/2.0E1);
	T(6,0) = (q4dot*q4dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*sin(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*sin(q4+theta)*(1.3E1/4.0E1)+q4dot*q5dot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*sin(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q5dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*sin(q4+theta)*(1.3E1/2.0E1);
	T(7,0) = (q4dot*q4dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q5dot*q5dot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+q5+theta)*(1.3E1/4.0E1)+(q4dot*q4dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q5dot*q5dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(q6dot*q6dot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)+(thetadot*thetadot)*cos(q4+q5+q6+theta)*(1.0/1.0E1)-(q4dot*q4dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(q5dot*q5dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(q6dot*q6dot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)-(thetadot*thetadot)*sin(q4+q5+q6+theta)*(1.0/1.0E1)+(q4dot*q4dot)*cos(q4+theta)*(1.3E1/4.0E1)+(thetadot*thetadot)*cos(q4+theta)*(1.3E1/4.0E1)+q4dot*q5dot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q5dot*thetadot*cos(q4+q5+theta)*(1.3E1/2.0E1)+q4dot*q5dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*q6dot*cos(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q5dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)+q6dot*thetadot*cos(q4+q5+q6+theta)*(1.0/5.0)-q4dot*q5dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q4dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q5dot*q6dot*sin(q4+q5+q6+theta)*(1.0/5.0)-q4dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)-q5dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)-q6dot*thetadot*sin(q4+q5+q6+theta)*(1.0/5.0)+q4dot*thetadot*cos(q4+theta)*(1.3E1/2.0E1);
	T(8,0) = (thetadot*thetadot)*sin(theta)*(-1.1E1/2.0E1)+(q7dot*q7dot)*sin(q7+q8+theta)*(9.0/2.0E1)+(q8dot*q8dot)*sin(q7+q8+theta)*(9.0/2.0E1)+(thetadot*thetadot)*sin(q7+q8+theta)*(9.0/2.0E1)+(q7dot*q7dot)*sin(q7+theta)*(1.0/4.0)+(thetadot*thetadot)*sin(q7+theta)*(1.0/4.0)+q7dot*q8dot*sin(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*sin(q7+q8+theta)*(9.0/1.0E1)+q8dot*thetadot*sin(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*sin(q7+theta)*(1.0/2.0);
	T(9,0) = (thetadot*thetadot)*cos(theta)*(-1.1E1/2.0E1)+(q7dot*q7dot)*cos(q7+q8+theta)*(9.0/2.0E1)+(q8dot*q8dot)*cos(q7+q8+theta)*(9.0/2.0E1)+(thetadot*thetadot)*cos(q7+q8+theta)*(9.0/2.0E1)+(q7dot*q7dot)*cos(q7+theta)*(1.0/4.0)+(thetadot*thetadot)*cos(q7+theta)*(1.0/4.0)+q7dot*q8dot*cos(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*cos(q7+q8+theta)*(9.0/1.0E1)+q8dot*thetadot*cos(q7+q8+theta)*(9.0/1.0E1)+q7dot*thetadot*cos(q7+theta)*(1.0/2.0);
	T(10,0) = (thetadot*thetadot)*sin(theta)*(-1.1E1/2.0E1)+(q10dot*q10dot)*sin(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*sin(q9+q10+theta)*(9.0/2.0E1)+(thetadot*thetadot)*sin(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*sin(q9+theta)*(1.0/4.0)+(thetadot*thetadot)*sin(q9+theta)*(1.0/4.0)+q10dot*q9dot*sin(q9+q10+theta)*(9.0/1.0E1)+q10dot*thetadot*sin(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*sin(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*sin(q9+theta)*(1.0/2.0);
	T(11,0) = (thetadot*thetadot)*cos(theta)*(-1.1E1/2.0E1)+(q10dot*q10dot)*cos(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*cos(q9+q10+theta)*(9.0/2.0E1)+(thetadot*thetadot)*cos(q9+q10+theta)*(9.0/2.0E1)+(q9dot*q9dot)*cos(q9+theta)*(1.0/4.0)+(thetadot*thetadot)*cos(q9+theta)*(1.0/4.0)+q10dot*q9dot*cos(q9+q10+theta)*(9.0/1.0E1)+q10dot*thetadot*cos(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*cos(q9+q10+theta)*(9.0/1.0E1)+q9dot*thetadot*cos(q9+theta)*(1.0/2.0);
	return T;
}
std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{

	std::vector<double> r_pos(2);
	r_pos[0] = 0.0;
	r_pos[1] = 0.0;

	double rIx = Robot_StateNDot_i.rIx;
	double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	if (strcmp(s,"rA")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"rB")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
  		r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"rC")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"rD")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
	  	r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"rE")==0)
	{
		r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	    r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rF")==0)
	{
		r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	    r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rG")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
     	r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rH")==0)
	{
		r_pos[0] = rIx+cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	    r_pos[1] = rIy-sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rI")==0)
	{
		r_pos[0] = rIx;
	    r_pos[1] = rIy;
	}
	if (strcmp(s,"rJ")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
		r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)-sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rK")==0)
	{
		r_pos[0] = rIx+cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	    r_pos[1] = rIy-sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1);
	}
	if (strcmp(s,"rL")==0)
	{
		r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	    r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rM")==0)
	{
		r_pos[0] = rIx+cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
		r_pos[1] = rIy-sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rN")==0)
	{
		r_pos[0] = rIx+cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	    r_pos[1] = rIy-sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)-sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1);
	}
	if (strcmp(s,"rT")==0)
	{
		r_pos[0] = rIx+cos(theta-PI*(1.0/2.0))*(7.0/1.0E1);
	    r_pos[1] = rIy-sin(theta-PI*(1.0/2.0))*(7.0/1.0E1);
	}
	if (strcmp(s,"rCOM")==0)
	{
		r_pos[0] = rIx-sin(q1+q2+theta)*1.678966789667897E-2-sin(q4+q5+theta)*1.678966789667897E-2-sin(q7+q8+theta)*8.717712177121771E-3-sin(q9+q10+theta)*8.717712177121771E-3-sin(q1+theta)*3.597785977859779E-2-sin(q4+theta)*3.597785977859779E-2-sin(q7+theta)*1.775830258302583E-2-sin(q9+theta)*1.775830258302583E-2+sin(theta)*2.506457564575646E-1-sqrt(2.0)*sin(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*sin(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(4.1E1)*sin(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*sin(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4;
	    r_pos[1] = rIy-cos(q1+q2+theta)*1.678966789667897E-2-cos(q4+q5+theta)*1.678966789667897E-2-cos(q7+q8+theta)*8.717712177121771E-3-cos(q9+q10+theta)*8.717712177121771E-3-cos(q1+theta)*3.597785977859779E-2-cos(q4+theta)*3.597785977859779E-2-cos(q7+theta)*1.775830258302583E-2-cos(q9+theta)*1.775830258302583E-2+cos(theta)*2.506457564575646E-1-sqrt(4.1E1)*cos(q1+q2+q3+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(4.1E1)*cos(q4+q5+q6+theta-8.960553845713439E-1)*3.690036900369004E-4-sqrt(2.0)*cos(q1+q2+q3+theta+PI*(1.0/4.0))*1.476014760147601E-3-sqrt(2.0)*cos(q4+q5+q6+theta+PI*(1.0/4.0))*1.476014760147601E-3;
	}
	return r_pos;
}
std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s)
{
	double PI = 3.1415926535897932384626;
	std::vector<double> T(2);

	T[0] = 0.0;
	T[1] = 0.0;

	// double rIx = Robot_StateNDot_i.rIx;
	// double rIy = Robot_StateNDot_i.rIy;
	double theta = Robot_StateNDot_i.theta;
	double q1 = Robot_StateNDot_i.q1;
	double q2 = Robot_StateNDot_i.q2;
	double q3 = Robot_StateNDot_i.q3;
	double q4 = Robot_StateNDot_i.q4;
	double q5 = Robot_StateNDot_i.q5;
	double q6 = Robot_StateNDot_i.q6;
	double q7 = Robot_StateNDot_i.q7;
	double q8 = Robot_StateNDot_i.q8;
	double q9 = Robot_StateNDot_i.q9;
	double q10 = Robot_StateNDot_i.q10;

	double rIxdot = Robot_StateNDot_i.rIxdot;
	double rIydot = Robot_StateNDot_i.rIydot;
	double thetadot = Robot_StateNDot_i.thetadot;
	double q1dot = Robot_StateNDot_i.q1dot;
	double q2dot = Robot_StateNDot_i.q2dot;
	double q3dot = Robot_StateNDot_i.q3dot;
	double q4dot = Robot_StateNDot_i.q4dot;
	double q5dot = Robot_StateNDot_i.q5dot;
	double q6dot = Robot_StateNDot_i.q6dot;
	double q7dot = Robot_StateNDot_i.q7dot;
	double q8dot = Robot_StateNDot_i.q8dot;
	double q9dot = Robot_StateNDot_i.q9dot;
	double q10dot = Robot_StateNDot_i.q10dot;

	if (strcmp(s,"vA")==0)
	{
		T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*sin(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q3dot*cos(q1+q2+q3+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"vB")==0)
	{
		T[0] = rIxdot-q1dot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(sin(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*sin(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
		T[1] = rIydot-q1dot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q1+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-q2dot*(cos(q1+q2+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q3dot*cos(q1+q2+q3+theta-PI*(5.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"vC")==0)
	{
		T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*sin(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
		T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(4.1E1)*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1))-sqrt(4.1E1)*q6dot*cos(q4+q5+q6+theta+PI*(1.0/2.0)-8.960553845713439E-1)*(1.0/4.0E1);
	}
	if (strcmp(s,"vD")==0)
	{
		T[0] = rIxdot-q4dot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(sin(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(sin(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*sin(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
		T[1] = rIydot-q4dot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-thetadot*(cos(q4+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-q5dot*(cos(q4+q5+theta+PI*(1.0/2.0))*(1.3E1/4.0E1)+sqrt(2.0)*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1))-sqrt(2.0)*q6dot*cos(q4+q5+q6+theta+PI*(3.0/4.0))*(1.0/1.0E1);
	}
	if (strcmp(s,"vE")==0)
	{
		T[0] = rIxdot-q7dot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*sin(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
		T[1] = rIydot-q7dot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q7+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q8dot*cos(q7+q8+theta+PI*(1.0/2.0))*(9.0/2.0E1);
	}
	if(strcmp(s,"vF")==0)
	{
		T[0] = rIxdot-q9dot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(sin(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+sin(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*sin(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
		T[1] = rIydot-q9dot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1))-thetadot*(cos(q9+theta+PI*(1.0/2.0))*(1.0/4.0)+cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1)+cos(theta-PI*(1.0/2.0))*(1.1E1/2.0E1))-q10dot*cos(q9+q10+theta+PI*(1.0/2.0))*(9.0/2.0E1);
	}
	if(strcmp(s,"vI")==0)
	{
		T[0] = rIxdot;
		T[1] = rIydot;
	}
	return T;
}
double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i)
{
	dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full, Jac_Full_Trans;
	Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
	std::vector<double> Robot_StateVec_i;
	Robot_StateVec_i = StateNDot2StateVec(Robot_StateNDot_i);
	dlib::matrix<double> Robot_StatedotDlib, Robot_StatedotDlib_Trans;	Robot_StatedotDlib = dlib::zeros_matrix<double>(13,1);
	for (int i = 0; i < 13; i++) {Robot_StatedotDlib(i) = Robot_StateVec_i[i+13];}
	Robot_StatedotDlib_Trans = dlib::trans(Robot_StatedotDlib);
	double T;
	T = 0.5 * Robot_StatedotDlib_Trans * D_q * Robot_StatedotDlib;
	return T;
}

void Dynamics_Matrices(const Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full)
{
	D_q = D_q_fn(Node_StateNDot);
	B_q = B_q_fn();
	C_q_qdot = C_q_qdot_fn(Node_StateNDot);
	Jac_Full = Jac_Full_fn(Node_StateNDot);
}

dlib::matrix<double> Dynamics_RHS_Matrix_fn(dlib::matrix<double> &Jac_Full, dlib::matrix<double> &B_q)
{
	dlib::matrix<double> Jac_Full_Trans, Dynamics_RHS_Matrix;
	Jac_Full_Trans = dlib::trans(Jac_Full);

	const int Dynamics_RHS_Matrix_Row = Jac_Full_Trans.nr();
	const int Dynamics_RHS_Matrix_Col = Jac_Full_Trans.nc() + B_q.nc();
	Dynamics_RHS_Matrix = dlib::zeros_matrix<double>(Dynamics_RHS_Matrix_Row, Dynamics_RHS_Matrix_Col);
	for (int i = 0; i < Dynamics_RHS_Matrix_Col; i++) {
		if (i<Jac_Full_Trans.nc())
		{

			dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(Jac_Full_Trans,i);
		}
		else{
			dlib::set_colm(Dynamics_RHS_Matrix, i) = dlib::colm(B_q,i-Jac_Full_Trans.nc());
		}
	}
	return Dynamics_RHS_Matrix;
}
