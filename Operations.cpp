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

extern dlib::matrix<double> xlow_vec;
extern dlib::matrix<double> xupp_vec;
extern dlib::matrix<double> ctrl_low_vec;
extern dlib::matrix<double> ctrl_upp_vec;
extern int Grids;
extern double Time_Seed;
extern double mini;

Robot_StateNDot::Robot_StateNDot(){// A default constructor
	rIx = 0;			rIy = 0.7230;			theta = -0.0900;
	q1 = 0.3768;		q2 = 0.0045;			q3 = -0.2913;			q4 = -1.0015;			q5 = 0.1500;
	q6 = 0.2698;		q7 = -0.6600;			q8 = -0.6251;			q9 = 0.6900;			q10 = -0.2951;
	rIxdot = 0.2000;	rIydot = -0.0605;		thetadot = -0.2100;
	q1dot = -0.1239;	q2dot = 1.3108;			q3dot = -0.9768;		q4dot = -1.4999;		q5dot = 2.0000;
	q6dot = -1.2999;	q7dot = 1.0000;			q8dot = -2.0000;		q9dot = -1.5708;		q10dot = -1.5000;
}
Robot_StateNDot::Robot_StateNDot(std::vector<double> &Robot_AngleNRate){	// An evaluated constructor
    rIx = Robot_AngleNRate[0];		rIy = Robot_AngleNRate[1];		theta = Robot_AngleNRate[2];
    q1 = Robot_AngleNRate[3];		q2 = Robot_AngleNRate[4];		q3 = Robot_AngleNRate[5];		q4 = Robot_AngleNRate[6];		q5 = Robot_AngleNRate[7];
    q6 = Robot_AngleNRate[8];		q7 = Robot_AngleNRate[9];		q8 = Robot_AngleNRate[10];		q9 = Robot_AngleNRate[11];		q10 = Robot_AngleNRate[12];

    rIxdot = Robot_AngleNRate[13];	rIydot = Robot_AngleNRate[14];	thetadot = Robot_AngleNRate[15];
    q1dot = Robot_AngleNRate[16];	q2dot = Robot_AngleNRate[17];	q3dot = Robot_AngleNRate[18];	q4dot = Robot_AngleNRate[19];	q5dot = Robot_AngleNRate[20];
    q6dot = Robot_AngleNRate[21];	q7dot = Robot_AngleNRate[22];	q8dot = Robot_AngleNRate[23];	q9dot = Robot_AngleNRate[24];	q10dot = Robot_AngleNRate[25];
}
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i)
{
	// This function is used to plot the robot configuration
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

	std::vector<double> x(2), y(2);
	// AB
	x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
	// CD
	x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// CJ
	x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// JD
	x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// KJ
	x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// IK
	x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
	// BG
	x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// AG
	x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// GH
	x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
	// HI
	x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
	// IT
	x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
	// LM
	x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
	// ME
	x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
	// LN
	x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
	// NF
	x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);
	const char* filename = "./init.png";
	std::cout<<"Saving result to "<<filename<<std::endl;
	plt::save(filename);
}
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name)
{
	// This function is used to plot the robot configuration
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rJ = Ang_Pos_fn(StateNDot_Init_i, "rJ");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rL = Ang_Pos_fn(StateNDot_Init_i, "rL");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

	std::vector<double> x(2), y(2);
	// AB
	x.at(0) = rA[0];	x.at(1) = rB[0];	y.at(0) = rA[1];	y.at(1) = rB[1];	plt::plot(x,y);
	// CD
	x.at(0) = rC[0];	x.at(1) = rD[0];	y.at(0) = rC[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// CJ
	x.at(0) = rC[0];	x.at(1) = rJ[0];	y.at(0) = rC[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// JD
	x.at(0) = rJ[0];	x.at(1) = rD[0];	y.at(0) = rJ[1];	y.at(1) = rD[1];	plt::plot(x,y);
	// KJ
	x.at(0) = rK[0];	x.at(1) = rJ[0];	y.at(0) = rK[1];	y.at(1) = rJ[1];	plt::plot(x,y);
	// IK
	x.at(0) = rI[0];	x.at(1) = rK[0];	y.at(0) = rI[1];	y.at(1) = rK[1];	plt::plot(x,y);
	// BG
	x.at(0) = rB[0];	x.at(1) = rG[0];	y.at(0) = rB[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// AG
	x.at(0) = rA[0];	x.at(1) = rG[0];	y.at(0) = rA[1];	y.at(1) = rG[1];	plt::plot(x,y);
	// GH
	x.at(0) = rG[0];	x.at(1) = rH[0];	y.at(0) = rG[1];	y.at(1) = rH[1];	plt::plot(x,y);
	// HI
	x.at(0) = rH[0];	x.at(1) = rI[0];	y.at(0) = rH[1];	y.at(1) = rI[1];	plt::plot(x,y);
	// IT
	x.at(0) = rI[0];	x.at(1) = rT[0];	y.at(0) = rI[1];	y.at(1) = rT[1];	plt::plot(x,y);
	// LM
	x.at(0) = rL[0];	x.at(1) = rM[0];	y.at(0) = rL[1];	y.at(1) = rM[1];	plt::plot(x,y);
	// ME
	x.at(0) = rM[0];	x.at(1) = rE[0];	y.at(0) = rM[1];	y.at(1) = rE[1];	plt::plot(x,y);
	// LN
	x.at(0) = rL[0];	x.at(1) = rN[0];	y.at(0) = rL[1];	y.at(1) = rN[1];	plt::plot(x,y);
	// NF
	x.at(0) = rN[0];	x.at(1) = rF[0];	y.at(0) = rN[1];	y.at(1) = rF[1];	plt::plot(x,y);

	std::string pre_filename = "./";
    std::string post_filename = ".png";
	std::string filename = pre_filename + name + post_filename;

	std::cout<<"Saving result to "<<filename<<std::endl;
	plt::save(filename);
}

std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2)
{
	int Vec_Length = vec1.size();

	std::vector<double> vec_minus;

	for (int i = 0; i < Vec_Length; i++)
	{
		vec_minus.push_back(vec1[i] - vec2[i]);
	}
	return vec_minus;
}

std::vector<double> Sigma2Pos(std::vector<double> &sigma, int EqOrIneq)
{
	std::vector<double> sigma_pos;
	if(EqOrIneq==0)
	{
		sigma_pos.push_back(sigma[0]);
		sigma_pos.push_back(sigma[0]);
		sigma_pos.push_back(sigma[1]);
		sigma_pos.push_back(sigma[1]);
		sigma_pos.push_back(sigma[2]);
		sigma_pos.push_back(sigma[3]);
	}
	else
	{
		sigma_pos.push_back(!sigma[0]);
		sigma_pos.push_back(!sigma[0]);
		sigma_pos.push_back(!sigma[1]);
		sigma_pos.push_back(!sigma[1]);
		sigma_pos.push_back(!sigma[2]);
		sigma_pos.push_back(!sigma[3]);
	}
	return sigma_pos;
}
std::vector<double> Sigma2Vel(std::vector<double> &sigma)
{	std::vector<double> sigma_pos;
	sigma_pos.push_back(sigma[0]);
	sigma_pos.push_back(sigma[0]);
	sigma_pos.push_back(sigma[0]);
	sigma_pos.push_back(sigma[0]);
	sigma_pos.push_back(sigma[1]);
	sigma_pos.push_back(sigma[1]);
	sigma_pos.push_back(sigma[1]);
	sigma_pos.push_back(sigma[1]);
	sigma_pos.push_back(sigma[2]);
	sigma_pos.push_back(sigma[2]);
	sigma_pos.push_back(sigma[3]);
	sigma_pos.push_back(sigma[3]);
	return sigma_pos;
}
dlib::matrix<double> Diag_Matrix_fn(std::vector<double> &diag_vec)
{
	// This function is used to generate a diagonal matrix within diag_vec to its diagonal elements
	int dim = diag_vec.size();
	dlib::matrix<double> Diag_Matrix;
	Diag_Matrix = dlib::zeros_matrix<double>(dim,dim);
	for (int i = 0; i < dim; i++)
	{
		Diag_Matrix(i,i) = diag_vec[i];
	}
	return Diag_Matrix;
}

dlib::matrix<double> Quadratic_Minus(dlib::matrix<double> &Mat_A, dlib::matrix<double> &Mat_B)
{
	// This function is used to take care of the situation where the SNOPT cannot be easily used to optimize the simple minus constraint
	dlib::matrix<double> Matrix_Minus, Matrix_result; Matrix_result = Mat_A;
	Matrix_Minus = Mat_A - Mat_B;
	for (int i = 0; i < Matrix_Minus.nr(); i++) {
		Matrix_result(i) = Matrix_Minus(i) * Matrix_Minus(i);
	}
	return Matrix_result;
}

std::vector<double> StateNDot2StateVec(const Robot_StateNDot &Robot_StateNDot_i)
{
	std::vector<double> StateVec(26);
	StateVec[0] = Robot_StateNDot_i.rIx;
	StateVec[1] = Robot_StateNDot_i.rIy;
	StateVec[2] = Robot_StateNDot_i.theta;
	StateVec[3] = Robot_StateNDot_i.q1;
	StateVec[4] = Robot_StateNDot_i.q2;
	StateVec[5] = Robot_StateNDot_i.q3;
	StateVec[6] = Robot_StateNDot_i.q4;
	StateVec[7] = Robot_StateNDot_i.q5;
	StateVec[8] = Robot_StateNDot_i.q6;
	StateVec[9] = Robot_StateNDot_i.q7;
	StateVec[10] = Robot_StateNDot_i.q8;
	StateVec[11] = Robot_StateNDot_i.q9;
	StateVec[12] = Robot_StateNDot_i.q10;

	StateVec[0+13] = Robot_StateNDot_i.rIxdot;
	StateVec[1+13] = Robot_StateNDot_i.rIydot;
	StateVec[2+13] = Robot_StateNDot_i.thetadot;
	StateVec[3+13] = Robot_StateNDot_i.q1dot;
	StateVec[4+13] = Robot_StateNDot_i.q2dot;
	StateVec[5+13] = Robot_StateNDot_i.q3dot;
	StateVec[6+13] = Robot_StateNDot_i.q4dot;
	StateVec[7+13] = Robot_StateNDot_i.q5dot;
	StateVec[8+13] = Robot_StateNDot_i.q6dot;
	StateVec[9+13] = Robot_StateNDot_i.q7dot;
	StateVec[10+13] = Robot_StateNDot_i.q8dot;
	StateVec[11+13] = Robot_StateNDot_i.q9dot;
	StateVec[12+13] = Robot_StateNDot_i.q10dot;
	return StateVec;
}
Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec)
{
	Robot_StateNDot Robot_StateNDot_i;
	Robot_StateNDot_i.rIx = StateVec[0];
	Robot_StateNDot_i.rIy = StateVec[1];
	Robot_StateNDot_i.theta = StateVec[2];
	Robot_StateNDot_i.q1 = StateVec[3];
	Robot_StateNDot_i.q2 = StateVec[4];
	Robot_StateNDot_i.q3 = StateVec[5];
	Robot_StateNDot_i.q4 = StateVec[6];
	Robot_StateNDot_i.q5 = StateVec[7];
	Robot_StateNDot_i.q6 = StateVec[8];
	Robot_StateNDot_i.q7 = StateVec[9];
	Robot_StateNDot_i.q8 = StateVec[10];
	Robot_StateNDot_i.q9 = StateVec[11];
	Robot_StateNDot_i.q10 = StateVec[12];

	Robot_StateNDot_i.rIxdot = StateVec[13];
	Robot_StateNDot_i.rIydot = StateVec[14];
	Robot_StateNDot_i.thetadot = StateVec[15];
	Robot_StateNDot_i.q1dot = StateVec[16];
	Robot_StateNDot_i.q2dot = StateVec[17];
	Robot_StateNDot_i.q3dot = StateVec[18];
	Robot_StateNDot_i.q4dot = StateVec[19];
	Robot_StateNDot_i.q5dot = StateVec[20];
	Robot_StateNDot_i.q6dot = StateVec[21];
	Robot_StateNDot_i.q7dot = StateVec[22];
	Robot_StateNDot_i.q8dot = StateVec[23];
	Robot_StateNDot_i.q9dot = StateVec[24];
	Robot_StateNDot_i.q10dot = StateVec[25];

	return Robot_StateNDot_i;

}

dlib::matrix<double> ONES_VECTOR_fn(int Dim)
{
	dlib::matrix<double> Ones_vector;
	Ones_vector = dlib::ones_matrix<double>(Dim,1);
	return Ones_vector;
}

std::vector<double> Opt_Soln_Load()
{
	// This function is used to load the computed optimal solution for data analysis
	std::vector<double> Opt_Seed;
	ifstream Opt_Soln_File;              // This is to read the initial angle and angular velocities
	Opt_Soln_File.open("Opt_Soln.txt");
	if(Opt_Soln_File.is_open())
	{
		double data_each_line = 0.0;

		while(Opt_Soln_File>>data_each_line)
		{
			Opt_Seed.push_back(data_each_line);
		}
		Opt_Soln_File.close();
	}
	return Opt_Seed;
}

void Sigma_TransNGoal(std::vector<double> & sigma_i, std::vector<double> & sigma_i_child,std::vector<double> &sigma_trans, std::vector<double> & sigma_goal, int &Self_Opt_Flag)
{	double sigma_result = 0.0;
	Self_Opt_Flag = 0;
	for (int i = 0; i < sigma_i.size(); i++) {
		sigma_result = sigma_result + sigma_i_child[i] - sigma_i[i];
	}
	if(sigma_result>0)
	{
		sigma_trans = sigma_i;
		sigma_goal = sigma_i_child;
	}
	else
	{
		sigma_trans = sigma_i_child;
		sigma_goal = sigma_i_child;
	}
	if(sigma_result==0)
	{
		Self_Opt_Flag = 1;
	}
}
