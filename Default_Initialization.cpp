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
double rIxlow = -Inf;                  double rIxupp = Inf;					double rIylow = -Inf;                  double rIyupp = Inf;
double thetalow = -pi;                 double thetaupp = pi;				double q1low = -2.1817;                double q1upp = 0.733;
double q2low = 0.0;                    double q2upp = 2.618;				double q3low = -1.3;                   double q3upp = 0.733;
double q4low = -2.1817;                double q4upp = 0.733;				double q5low = 0.0;                    double q5upp = 2.618;
double q6low = -1.3;                   double q6upp = 0.733;				double q7low = -3.14;                  double q7upp = 1.047;
double q8low = -2.391;                 double q8upp = 0.0;					double q9low = -3.14;                  double q9upp = 1.047;
double q10low = -2.391;                double q10upp = 0.0;					double AngRateLow = -3.0;              double AngRateHgh = 3.0;

double rIxdotlow = -Inf;               double rIxdotupp = Inf;				double rIydotlow = -Inf;               double rIydotupp = Inf;
double thetadotlow = -Inf;             double thetadotupp = Inf;			double q1dotlow = AngRateLow;          double q1dotupp = AngRateHgh;
double q2dotlow = AngRateLow;          double q2dotupp = AngRateHgh;		double q3dotlow = AngRateLow;          double q3dotupp = AngRateHgh;
double q4dotlow = AngRateLow;          double q4dotupp = AngRateHgh;		double q5dotlow = AngRateLow;          double q5dotupp = AngRateHgh;
double q6dotlow = AngRateLow;          double q6dotupp = AngRateHgh;		double q7dotlow = AngRateLow;          double q7dotupp = AngRateHgh;
double q8dotlow = AngRateLow;          double q8dotupp = AngRateHgh;		double q9dotlow = AngRateLow;          double q9dotupp = AngRateHgh;
double q10dotlow = AngRateLow;         double q10dotupp = AngRateHgh;

double tau1_max = 100;             		double tau2_max = 100;				double tau3_max = 100;					double tau4_max = 100;
double tau5_max = 100;             		double tau6_max = 100;				double tau7_max = 60;              		double tau8_max = 50;
double tau9_max = 60;             		double tau10_max = 50;

extern dlib::matrix<double> xlow_vec;                       extern dlib::matrix<double> xupp_vec;
extern dlib::matrix<double> ctrl_low_vec;                   extern dlib::matrix<double> ctrl_upp_vec;
extern int Grids;                                           extern double Time_Seed;
extern double mini;                                         extern dlib::matrix<double>  Envi_Map;
extern dlib::matrix<double> Envi_Map_Normal;                extern dlib::matrix<double> Envi_Map_Tange; // The Normal and tangential vector of the plane

void Envi_Map_Defi()
{
	// 	This function is used to define the environment map for the simulation
	// Whenever this function gets called, it will return the array with the
	// environment obstacle information
	//
	// This is the default flat ground
	// This map is defined in a polyline manner with the first value denoting
	// the line length and the second value denoting the relative angle
	Envi_Map = dlib::ones_matrix<double>(2,4);
	Envi_Map(0,0) = -5.0;
	Envi_Map(0,1) = 0.0;
	Envi_Map(0,2) = 5.0;
	Envi_Map(0,3) = 0.0;

	Envi_Map(1,0) = 5.0;
	Envi_Map(1,1) = 0.0;
	Envi_Map(1,2) = 5.0;
	Envi_Map(1,3) = 3.0;
}

void Envi_Map_Normal_Cal(dlib::matrix<double> &Envi_Map)
{
	// This function is used to calculate the surface normal vector and tangential vector
	int NumOfObs = Envi_Map.nr();
	int Dimension = Envi_Map.nc()/2;
	Envi_Map_Normal = dlib::zeros_matrix<double>(NumOfObs,Dimension);
	Envi_Map_Tange = dlib::zeros_matrix<double>(NumOfObs,Dimension);
	double Envi_Map_Edge_A_x, Envi_Map_Edge_A_y, Envi_Map_Edge_B_x, Envi_Map_Edge_B_y;
	double Slope_Angle;
	for (int i = 0; i < NumOfObs; i++)
	{
		Envi_Map_Edge_A_x = Envi_Map(i, 0);
		Envi_Map_Edge_A_y = Envi_Map(i, 1);
		Envi_Map_Edge_B_x = Envi_Map(i, 2);
		Envi_Map_Edge_B_y = Envi_Map(i, 3);
		Slope_Angle = atan2(Envi_Map_Edge_B_y - Envi_Map_Edge_A_y, Envi_Map_Edge_B_x - Envi_Map_Edge_A_x);
		Envi_Map_Tange(i,0) = cos(Slope_Angle);
		Envi_Map_Tange(i,1) = sin(Slope_Angle);

		Envi_Map_Normal(i,0) = cos(Slope_Angle + pi/2.0);
		Envi_Map_Normal(i,1) = sin(Slope_Angle + pi/2.0);
	}
	return;
}

std::vector<double> Default_Init(const std::vector<double> &sigma_i)
{
	// This function is used to initialize the whole optimization process
	// First, is to substitute the map info into the Envi_Map matrix
	// Second, is to give the proper bounds to the variables to be optimized
	// Thrid, is to generate a kinematically feasible initial robot state

	Envi_Map_Defi();
	// Then it is to compute the normal and tangential vector of the map
	Envi_Map_Normal_Cal(Envi_Map);
	xlow_vec = dlib::zeros_matrix<double>(26,1);		xupp_vec = dlib::zeros_matrix<double>(26,1);
	ctrl_low_vec = dlib::zeros_matrix<double>(10,1);	ctrl_upp_vec = dlib::matrix<double>(10,1) ;
	xlow_vec(0) = rIxlow; 					xupp_vec(0) = rIxupp;
	xlow_vec(1) = rIylow; 					xupp_vec(1) = rIyupp;
	xlow_vec(2) = thetalow; 				xupp_vec(2) = thetaupp;
	xlow_vec(3) = q1low; 					xupp_vec(3) = q1upp;
	xlow_vec(4) = q2low; 					xupp_vec(4) = q2upp;
	xlow_vec(5) = q3low; 					xupp_vec(5) = q3upp;
	xlow_vec(6) = q4low; 					xupp_vec(6) = q4upp;
	xlow_vec(7) = q5low; 					xupp_vec(7) = q5upp;
	xlow_vec(8) = q6low; 					xupp_vec(8) = q6upp;
	xlow_vec(9) = q7low; 					xupp_vec(9) = q7upp;
	xlow_vec(10) = q8low; 					xupp_vec(10) = q8upp;
	xlow_vec(11) = q9low; 					xupp_vec(11) = q9upp;
	xlow_vec(12) = q10low; 					xupp_vec(12) = q10upp;
	xlow_vec(0+13) = rIxdotlow; 			xupp_vec(0+13) = rIxdotupp;
	xlow_vec(1+13) = rIydotlow; 			xupp_vec(1+13) = rIydotupp;
	xlow_vec(2+13) = thetadotlow; 			xupp_vec(2+13) = thetadotupp;
	xlow_vec(3+13) = q1dotlow; 				xupp_vec(3+13) = q1dotupp;
	xlow_vec(4+13) = q2dotlow; 				xupp_vec(4+13) = q2dotupp;
	xlow_vec(5+13) = q3dotlow; 				xupp_vec(5+13) = q3dotupp;
	xlow_vec(6+13) = q4dotlow; 				xupp_vec(6+13) = q4dotupp;
	xlow_vec(7+13) = q5dotlow; 				xupp_vec(7+13) = q5dotupp;
	xlow_vec(8+13) = q6dotlow; 				xupp_vec(8+13) = q6dotupp;
	xlow_vec(9+13) = q7dotlow; 				xupp_vec(9+13) = q7dotupp;
	xlow_vec(10+13) = q8dotlow; 			xupp_vec(10+13) = q8dotupp;
	xlow_vec(11+13) = q9dotlow; 			xupp_vec(11+13) = q9dotupp;
	xlow_vec(12+13) = q10dotlow; 			xupp_vec(12+13) = q10dotupp;

	ctrl_low_vec(0) = -tau1_max;			ctrl_upp_vec(0) = -ctrl_low_vec(0);
	ctrl_low_vec(1) = -tau2_max;			ctrl_upp_vec(1) = -ctrl_low_vec(1);
	ctrl_low_vec(2) = -tau3_max;			ctrl_upp_vec(2) = -ctrl_low_vec(2);
	ctrl_low_vec(3) = -tau4_max;			ctrl_upp_vec(3) = -ctrl_low_vec(3);
	ctrl_low_vec(4) = -tau5_max;			ctrl_upp_vec(4) = -ctrl_low_vec(4);
	ctrl_low_vec(5) = -tau6_max;			ctrl_upp_vec(5) = -ctrl_low_vec(5);
	ctrl_low_vec(6) = -tau7_max;			ctrl_upp_vec(6) = -ctrl_low_vec(6);
	ctrl_low_vec(7) = -tau8_max;			ctrl_upp_vec(7) = -ctrl_low_vec(7);
	ctrl_low_vec(8) = -tau9_max;			ctrl_upp_vec(8) = -ctrl_low_vec(8);
	ctrl_low_vec(9) = -tau10_max;			ctrl_upp_vec(9) = -ctrl_low_vec(9);

	vector<double> Robot_State_Init;
	ifstream Initial_Robot_State_File;              // This is to read the initial angle and angular velocities
	Initial_Robot_State_File.open("robot_angle_init.txt");
	if(Initial_Robot_State_File.is_open())
	{
		double data_each_line = 0.0;
		while(Initial_Robot_State_File>>data_each_line)
		{
			Robot_State_Init.push_back(data_each_line);
		}
		Initial_Robot_State_File.close();
	}
	else
	{
		printf("Unable to open robot_angle_init.txt file!\n");
	}
	Initial_Robot_State_File.open("robot_velocity_init.txt");
	if(Initial_Robot_State_File.is_open())
	{
		double data_each_line = 0.0;
		while(Initial_Robot_State_File>>data_each_line)
		{
			Robot_State_Init.push_back(data_each_line);
		}
		Initial_Robot_State_File.close();
	}
	else
	{
		printf("Unable to open robot_velocity_init.txt file!\n");
	}

	Tree_Node RootNode;
	RootNode.Node_StateNDot = StateVec2StateNDot(Robot_State_Init);
	RootNode.sigma = sigma_i;
	Structure_P.Node_i = RootNode;
	// Robot_Plot_fn(RootNode.Node_StateNDot);
	// // If the default configuration would like to be viewed
	// Robot_StateNDot Robot_StateNDot_init(Robot_State_Init);
	// std::string input_name = "init_given";
	// Robot_Plot_fn(Robot_StateNDot_init,input_name);

	snoptProblem Default_Init_Pr;                     // This is the name of the Optimization problem
	// Allocate and initialize
	std:vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
	integer n = Robot_State_Init.size();
	Default_Init_Pr_ObjNConstraint(Robot_State_Init, ObjNConstraint_Val, ObjNConstraint_Type);
	integer neF = ObjNConstraint_Val.size();     							  // 1 objective function
	integer lenA  =  n * neF;                         // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

	integer *iAfun = new integer[lenA];              //
	integer *jAvar = new integer[lenA];
	doublereal *A  = new doublereal[lenA];

	integer lenG   = lenA;
	integer *iGfun = new integer[lenG];
	integer *jGvar = new integer[lenG];

	doublereal *x      = new doublereal[n];
	doublereal *xlow   = new doublereal[n];
	doublereal *xupp   = new doublereal[n];
	doublereal *xmul   = new doublereal[n];
	integer    *xstate = new    integer[n];

	doublereal *F      = new doublereal[neF];
	doublereal *Flow   = new doublereal[neF];
	doublereal *Fupp   = new doublereal[neF];
	doublereal *Fmul   = new doublereal[neF];
	integer    *Fstate = new integer[neF];

	integer nxnames = 1;
	integer nFnames = 1;
	char *xnames = new char[nxnames*8];
	char *Fnames = new char[nFnames*8];

	integer    ObjRow = 0;
	doublereal ObjAdd = 0;

	// Set the upper and lower bounds.
	for (int i = 0; i < n; i++) {
		xlow[i] = xlow_vec(i);
		xupp[i] = xupp_vec(i);
		xstate[i] = 0.0;
		x[i] = Robot_State_Init[i];  	// Initial guess
	}

	for(int i = 0; i<neF; i++){
		// The lower bound is the same
		Flow[i] = 0.0;
		if(ObjNConstraint_Type[i]>0)	// Inequality constraint
		{	Fupp[i] = Inf;}
		else{
			Fupp[i] = 0.0;}
	}

	// Load the data for ToyProb ...
	Default_Init_Pr.setPrintFile  ( "Default_Init_Pr.out" );
	Default_Init_Pr.setProblemSize( n, neF );
	Default_Init_Pr.setObjective  ( ObjRow, ObjAdd );
	Default_Init_Pr.setA          ( lenA, iAfun, jAvar, A );
	Default_Init_Pr.setG          ( lenG, iGfun, jGvar );
	Default_Init_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Default_Init_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Default_Init_Pr.setXNames     ( xnames, nxnames );
	Default_Init_Pr.setFNames     ( Fnames, nFnames );
	Default_Init_Pr.setProbName   ( "Default_Init_Pr" );
	Default_Init_Pr.setUserFun    ( Default_Init_Pr_);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Default_Init_Pr.computeJac    ();
	Default_Init_Pr.setIntParameter( "Derivative option", 0 );
	integer Cold = 0, Basis = 1, Warm = 2;
	Default_Init_Pr.solve( Cold );

	// Take the value out from x
	for (int i = 0; i < n; i++)
	{
		Robot_State_Init[i] = x[i];
	}
	// Robot_StateNDot Init_Opt_vec(Robot_State_Init);
	// Robot_Plot_fn(Init_Opt_vec);

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;		 delete []Flow;	  delete []Fupp;
	delete []Fmul;	 delete []Fstate;

	delete []xnames; delete []Fnames;

	return Robot_State_Init;
}

int Default_Init_Pr_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
	     integer    iu[],    integer *leniu,
	     doublereal ru[],    integer *lenru )
{
	std::vector<double> ObjNConstraint_Val,ObjNConstraint_Type;

	// Initial guess of the robot configurations
	std::vector<double> Robot_State_Init = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
	std::vector<double> Robot_State_Opt;
	for (int i = 0; i < 26; i++)
	{
		Robot_State_Opt.push_back(x[i]);
	}

	Default_Init_Pr_ObjNConstraint(Robot_State_Opt, ObjNConstraint_Val,ObjNConstraint_Type);
	for (int i = 0; i < ObjNConstraint_Val.size(); i++){
		F[i] = ObjNConstraint_Val[i];}
	return 0;
}
void Default_Init_Pr_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);
	std::vector<double> Robostate_ref = StateNDot2StateVec(Structure_P.Node_i.Node_StateNDot);
	std::vector<double> Robostate_offset = Vec_Minus(Opt_Seed, Robostate_ref);
	double Robostate_offset_val = 0.0;
	for (int i = 0; i < Robostate_offset.size(); i++){
		Robostate_offset_val = Robostate_offset_val + Robostate_offset[i] * Robostate_offset[i];
	}
	std::vector<double> sigma = Structure_P.Node_i.sigma;
	ObjNConstraint_Val.push_back(Robostate_offset_val);
	ObjNConstraint_Type.push_back(1);

	dlib::matrix<double,6,1> End_Effector_Dist;
	std::vector<int> End_Effector_Obs(6);

	End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

	dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Matrix_result;
	std::vector<double> sigma_temp;
	sigma_temp = Sigma2Pos(sigma, 0);
	Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Pos(sigma, 1);
	Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Vel(sigma);
	Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

	// 1. Active constraints have to be satisfied: Position and Velocity
	Matrix_result = Eqn_Pos_Matrix * End_Effector_Dist;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	Matrix_result = Eqn_Vel_Matrix * End_Effector_Vel;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// 2. Inactive constraints have to be strictly away from the obstacle
	dlib::matrix<double> ones_vector, temp_matrix;
	ones_vector = ONES_VECTOR_fn(6);
	Matrix_result = Ineqn_Pos_Matrix * (End_Effector_Dist - ones_vector * mini);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

	// 3. Middle joints have to be strictly away from the obs
	temp_matrix = Middle_Joint_Obs_Dist_Fn(StateNDot_Init_i);
	Matrix_result = temp_matrix - ones_vector * mini;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
	return;
}
