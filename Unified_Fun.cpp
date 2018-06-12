#include <stdio.h>
#include <string.h>
#include <iostream>
#include "snopt.hh"
#include "Unified_Header.h"
#include "snoptProblem.hh"
#include <fstream>
#include <cmath>
#include <dlib/matrix.h>
#include <algorithm>

// Pre-define of the bounds
double Inf = 1.1e20;					double PI = 3.1415926535897932384626;					double pi = PI;
// There are three types of variables in this optimization problem: robot state, control torques and contact forces
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

dlib::matrix<double,26,1> xlow_vec;					dlib::matrix<double,26,1> xupp_vec;
dlib::matrix<double,10,1> ctrl_low_vec;				dlib::matrix<double,10,1> ctrl_upp_vec;
dlib::matrix<double>  Envi_Map;						dlib::matrix<double> Envi_Map_Normal, Envi_Map_Tange; // The Normal and tangential vector of the plane
/**
 * Some global values are defined
 * Description
 */
double mini = 0.05;			int Grids = 10;			double mu = 0.5;
std::vector<Tree_Node_Ptr> All_Nodes;				// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		    // The kinetic energy of each nodes

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
void Add_Node2Tree(Tree_Node &Current_Node)
{
	// This function will add the current node to the All_Nodes vector
	All_Nodes.push_back(&Current_Node);
	Frontier_Nodes.push_back(&Current_Node);
	Frontier_Nodes_Cost.push_back(Current_Node.KE);
}
int Minimum_Index(std::vector<double> &Given_vec)
{
	int index = 0;
	for(int i = 1; i < Given_vec.size(); i++)
	{
		if(Given_vec[i] < Given_vec[index])
		index = i;
	}
	return index;
}
Tree_Node Pop_Node()
{
	// This function will pop the node outfrom the current Frontier according to the kinetic energy
	int Min_Ind = Minimum_Index(Frontier_Nodes_Cost);
	Tree_Node Current_Node = *Frontier_Nodes[Min_Ind];
	Frontier_Nodes.erase(Frontier_Nodes.begin()+Min_Ind);
	Frontier_Nodes_Cost.erase(Frontier_Nodes_Cost.begin()+Min_Ind);
	return Current_Node;
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
std::vector<double> StateNDot2StateVec(Robot_StateNDot &Robot_StateNDot_i)
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
dlib::matrix<double> Middle_Joint_Obs_Dist_Fn(Robot_StateNDot &StateNDot_Init_i){
	dlib::matrix<double> Middle_Joint_Obs_Dist_Matrix;
	Middle_Joint_Obs_Dist_Matrix = dlib::ones_matrix<double>(6,1);
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");
	int Obs_Dist_Index;		double rH_Obs, rK_Obs, rM_Obs, rN_Obs, rI_Obs, rT_Obs;
	Obs_Dist_Fn(rH, rH_Obs, Obs_Dist_Index);
	Obs_Dist_Fn(rK, rK_Obs, Obs_Dist_Index);
	Obs_Dist_Fn(rM, rM_Obs, Obs_Dist_Index);
	Obs_Dist_Fn(rN, rN_Obs, Obs_Dist_Index);
	Obs_Dist_Fn(rI, rI_Obs, Obs_Dist_Index);
	Obs_Dist_Fn(rT, rT_Obs, Obs_Dist_Index);

	Middle_Joint_Obs_Dist_Matrix(0) = rH_Obs;
	Middle_Joint_Obs_Dist_Matrix(1) = rK_Obs;
	Middle_Joint_Obs_Dist_Matrix(2) = rM_Obs;
	Middle_Joint_Obs_Dist_Matrix(3) = rN_Obs;
	Middle_Joint_Obs_Dist_Matrix(4) = rI_Obs;
	Middle_Joint_Obs_Dist_Matrix(5) = rT_Obs;
	return Middle_Joint_Obs_Dist_Matrix;}
dlib::matrix<double> ONES_VECTOR_fn(int Dim){
	dlib::matrix<double> Ones_vector;
	Ones_vector = dlib::ones_matrix<double>(Dim,1);
	return Ones_vector;}
void ObjNConstraint_ValNType_Update(dlib::matrix<double> &Matrix_result, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type, int Constraint_Type)
{
	// This function is used to update the value of ObjNConstraint_Val and Type
	for (int i = 0; i < Matrix_result.nr(); i++)
	{
		ObjNConstraint_Val.push_back(Matrix_result(i));
		ObjNConstraint_Type.push_back(Constraint_Type);
	}
}
void Obs_Dist_Fn(std::vector<double> &r_Pos, double &Obs_Dist, int &Obs_Dist_Index)
{
	// 	This function is used to calculate the relative distance between the robot end effector and the nearby environment
	double Temp_Edge_x, Temp_Edge_y, Temp_offset_x, Temp_offset_y, Normal_vector_i_x, Normal_vector_i_y;
	std::vector<double> Obs_Dist_vec;
	for (int i = 0; i < Envi_Map.nr(); i++)
	{
		Temp_Edge_x = Envi_Map(i,0);
		Temp_Edge_y = Envi_Map(i,1);
		Temp_offset_x = r_Pos[0] - Temp_Edge_x;
		Temp_offset_y = r_Pos[1] - Temp_Edge_y;
		Normal_vector_i_x = Envi_Map_Normal(i,0);
		Normal_vector_i_y = Envi_Map_Normal(i,1);
		Obs_Dist_vec.push_back(Temp_offset_x * Normal_vector_i_x +Temp_offset_y * Normal_vector_i_y);
	}
	Obs_Dist = *std::min(Obs_Dist_vec.begin(), Obs_Dist_vec.begin() + Obs_Dist_vec.size());
	Obs_Dist_Index = *std::min_element(Obs_Dist_vec.begin(), Obs_Dist_vec.begin() + Obs_Dist_vec.size());
	return;
}
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

	double T;
	T = (q1dot*q1dot)*cos(q2)*2.9575E-1+(q1dot*q1dot)*cos(q3)*(1.3E1/5.0E2)+(q2dot*q2dot)*cos(q3)*(1.3E1/5.0E2)+(q4dot*q4dot)*cos(q5)*2.9575E-1+(q4dot*q4dot)*cos(q6)*(1.3E1/5.0E2)+(q5dot*q5dot)*cos(q6)*(1.3E1/5.0E2)+(q7dot*q7dot)*cos(q8)*(6.3E1/3.2E2)+(q9dot*q9dot)*cos(q10)*(6.3E1/3.2E2)+(thetadot*thetadot)*cos(q2)*2.9575E-1+(thetadot*thetadot)*cos(q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5)*2.9575E-1+(thetadot*thetadot)*cos(q6)*(1.3E1/5.0E2)-(thetadot*thetadot)*cos(q7)*6.25625E-1+(thetadot*thetadot)*cos(q8)*(6.3E1/3.2E2)-(thetadot*thetadot)*cos(q9)*6.25625E-1+(thetadot*thetadot)*cos(q10)*(6.3E1/3.2E2)-(q1dot*q1dot)*sin(q3)*(1.3E1/5.0E2)-(q2dot*q2dot)*sin(q3)*(1.3E1/5.0E2)-(q4dot*q4dot)*sin(q6)*(1.3E1/5.0E2)-(q5dot*q5dot)*sin(q6)*(1.3E1/5.0E2)-(thetadot*thetadot)*sin(q3)*(1.3E1/5.0E2)-(thetadot*thetadot)*sin(q6)*(1.3E1/5.0E2)+q10dot*q9dot*1.954166666666667E-1+q1dot*q2dot*2.785E-1+q1dot*q3dot*(1.0/4.0E1)+q2dot*q3dot*(1.0/4.0E1)+q4dot*q5dot*2.785E-1+q4dot*q6dot*(1.0/4.0E1)+q5dot*q6dot*(1.0/4.0E1)+q7dot*q8dot*1.954166666666667E-1+q10dot*thetadot*1.954166666666667E-1+q1dot*thetadot*8.418333333333333E-1+q2dot*thetadot*2.785E-1+q3dot*thetadot*(1.0/4.0E1)+q4dot*thetadot*8.418333333333333E-1+q5dot*thetadot*2.785E-1+q6dot*thetadot*(1.0/4.0E1)+q7dot*thetadot*4.824166666666667E-1+q8dot*thetadot*1.954166666666667E-1+q9dot*thetadot*4.824166666666667E-1+(q10dot*q10dot)*9.770833333333333E-2+(q1dot*q1dot)*4.209166666666667E-1+(q2dot*q2dot)*1.3925E-1+(q3dot*q3dot)*(1.0/8.0E1)+(q4dot*q4dot)*4.209166666666667E-1+(q5dot*q5dot)*1.3925E-1+(q6dot*q6dot)*(1.0/8.0E1)+(q7dot*q7dot)*2.412083333333333E-1+(q8dot*q8dot)*9.770833333333333E-2+(q9dot*q9dot)*2.412083333333333E-1+(rIxdot*rIxdot)*(2.71E2/1.0E1)+(rIydot*rIydot)*(2.71E2/1.0E1)+(thetadot*thetadot)*4.3795+(q1dot*q1dot)*cos(q2+q3)*(1.3E1/5.0E2)+(q4dot*q4dot)*cos(q5+q6)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q2+q3)*(1.3E1/5.0E2)+(thetadot*thetadot)*cos(q5+q6)*(1.3E1/5.0E2)-(thetadot*thetadot)*cos(q7+q8)*4.33125E-1-(thetadot*thetadot)*cos(q9+q10)*4.33125E-1-(q1dot*q1dot)*sin(q2+q3)*(1.3E1/5.0E2)-(q4dot*q4dot)*sin(q5+q6)*(1.3E1/5.0E2)-(thetadot*thetadot)*sin(q2+q3)*(1.3E1/5.0E2)-(thetadot*thetadot)*sin(q5+q6)*(1.3E1/5.0E2)+q10dot*q9dot*cos(q10)*(6.3E1/3.2E2)+q1dot*q2dot*cos(q2)*2.9575E-1+q1dot*q2dot*cos(q3)*(1.3E1/2.5E2)+q1dot*q3dot*cos(q3)*(1.3E1/5.0E2)+q2dot*q3dot*cos(q3)*(1.3E1/5.0E2)+q4dot*q5dot*cos(q5)*2.9575E-1+q4dot*q5dot*cos(q6)*(1.3E1/2.5E2)+q4dot*q6dot*cos(q6)*(1.3E1/5.0E2)+q5dot*q6dot*cos(q6)*(1.3E1/5.0E2)+q7dot*q8dot*cos(q8)*(6.3E1/3.2E2)+q10dot*thetadot*cos(q10)*(6.3E1/3.2E2)+q1dot*thetadot*cos(q2)*5.915E-1+q1dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q2dot*thetadot*cos(q2)*2.9575E-1+q2dot*thetadot*cos(q3)*(1.3E1/2.5E2)+q3dot*thetadot*cos(q3)*(1.3E1/5.0E2)+q4dot*thetadot*cos(q5)*5.915E-1+q4dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q5dot*thetadot*cos(q5)*2.9575E-1+q5dot*thetadot*cos(q6)*(1.3E1/2.5E2)+q6dot*thetadot*cos(q6)*(1.3E1/5.0E2)-q7dot*thetadot*cos(q7)*6.25625E-1+q7dot*thetadot*cos(q8)*(6.3E1/1.6E2)+q8dot*thetadot*cos(q8)*(6.3E1/3.2E2)-q9dot*thetadot*cos(q9)*6.25625E-1+q9dot*thetadot*cos(q10)*(6.3E1/1.6E2)+rIxdot*thetadot*cos(theta)*1.3585E1+sqrt(4.1E1)*(q1dot*q1dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(q2dot*q2dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(q3dot*q3dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(q4dot*q4dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(q5dot*q5dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(q6dot*q6dot)*cos(8.960553845713439E-1)*(1.0/1.0E3)+sqrt(4.1E1)*(thetadot*thetadot)*cos(8.960553845713439E-1)*(1.0/5.0E2)-q1dot*q2dot*sin(q3)*(1.3E1/2.5E2)-q1dot*q3dot*sin(q3)*(1.3E1/5.0E2)-q2dot*q3dot*sin(q3)*(1.3E1/5.0E2)-q4dot*q5dot*sin(q6)*(1.3E1/2.5E2)-q4dot*q6dot*sin(q6)*(1.3E1/5.0E2)-q5dot*q6dot*sin(q6)*(1.3E1/5.0E2)-q1dot*thetadot*sin(q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q3)*(1.3E1/2.5E2)-q3dot*thetadot*sin(q3)*(1.3E1/5.0E2)-q4dot*thetadot*sin(q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q6)*(1.3E1/2.5E2)-q6dot*thetadot*sin(q6)*(1.3E1/5.0E2)-rIydot*thetadot*sin(theta)*1.3585E1-sqrt(4.1E1)*(q1dot*q1dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(q2dot*q2dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(q3dot*q3dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(q4dot*q4dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(q5dot*q5dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(q6dot*q6dot)*sin(8.960553845713439E-1)*(1.0/1.0E3)-sqrt(4.1E1)*(thetadot*thetadot)*sin(8.960553845713439E-1)*(1.0/5.0E2)-q10dot*rIxdot*cos(q9+q10+theta)*(6.3E1/8.0E1)-q1dot*rIxdot*cos(q1+q2+theta)*(9.1E1/1.0E2)-q2dot*rIxdot*cos(q1+q2+theta)*(9.1E1/1.0E2)-q4dot*rIxdot*cos(q4+q5+theta)*(9.1E1/1.0E2)-q5dot*rIxdot*cos(q4+q5+theta)*(9.1E1/1.0E2)-q7dot*rIxdot*cos(q7+q8+theta)*(6.3E1/8.0E1)-q8dot*rIxdot*cos(q7+q8+theta)*(6.3E1/8.0E1)-q9dot*rIxdot*cos(q9+q10+theta)*(6.3E1/8.0E1)-rIxdot*thetadot*cos(q1+q2+theta)*(9.1E1/1.0E2)-rIxdot*thetadot*cos(q4+q5+theta)*(9.1E1/1.0E2)-rIxdot*thetadot*cos(q7+q8+theta)*(6.3E1/8.0E1)-rIxdot*thetadot*cos(q9+q10+theta)*(6.3E1/8.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q4dot*q4dot)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*cos(q5+q6-8.960553845713439E-1)*6.5E-3+q10dot*rIydot*sin(q9+q10+theta)*(6.3E1/8.0E1)+q1dot*rIydot*sin(q1+q2+theta)*(9.1E1/1.0E2)+q2dot*rIydot*sin(q1+q2+theta)*(9.1E1/1.0E2)+q4dot*rIydot*sin(q4+q5+theta)*(9.1E1/1.0E2)+q5dot*rIydot*sin(q4+q5+theta)*(9.1E1/1.0E2)+q7dot*rIydot*sin(q7+q8+theta)*(6.3E1/8.0E1)+q8dot*rIydot*sin(q7+q8+theta)*(6.3E1/8.0E1)+q9dot*rIydot*sin(q9+q10+theta)*(6.3E1/8.0E1)+rIydot*thetadot*sin(q1+q2+theta)*(9.1E1/1.0E2)+rIydot*thetadot*sin(q4+q5+theta)*(9.1E1/1.0E2)+rIydot*thetadot*sin(q7+q8+theta)*(6.3E1/8.0E1)+rIydot*thetadot*sin(q9+q10+theta)*(6.3E1/8.0E1)-q1dot*rIxdot*cos(q1+q2+q3+theta)*(2.0/2.5E1)-q2dot*rIxdot*cos(q1+q2+q3+theta)*(2.0/2.5E1)-q3dot*rIxdot*cos(q1+q2+q3+theta)*(2.0/2.5E1)-q4dot*rIxdot*cos(q4+q5+q6+theta)*(2.0/2.5E1)-q5dot*rIxdot*cos(q4+q5+q6+theta)*(2.0/2.5E1)-q6dot*rIxdot*cos(q4+q5+q6+theta)*(2.0/2.5E1)+q1dot*rIydot*cos(q1+q2+q3+theta)*(2.0/2.5E1)+q2dot*rIydot*cos(q1+q2+q3+theta)*(2.0/2.5E1)+q3dot*rIydot*cos(q1+q2+q3+theta)*(2.0/2.5E1)+q4dot*rIydot*cos(q4+q5+q6+theta)*(2.0/2.5E1)+q5dot*rIydot*cos(q4+q5+q6+theta)*(2.0/2.5E1)+q6dot*rIydot*cos(q4+q5+q6+theta)*(2.0/2.5E1)-rIxdot*thetadot*cos(q1+q2+q3+theta)*(2.0/2.5E1)-rIxdot*thetadot*cos(q4+q5+q6+theta)*(2.0/2.5E1)+rIydot*thetadot*cos(q1+q2+q3+theta)*(2.0/2.5E1)+rIydot*thetadot*cos(q4+q5+q6+theta)*(2.0/2.5E1)+q1dot*rIxdot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q2dot*rIxdot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q3dot*rIxdot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q4dot*rIxdot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q5dot*rIxdot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q6dot*rIxdot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q1dot*rIydot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q2dot*rIydot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q3dot*rIydot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+q4dot*rIydot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q5dot*rIydot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q6dot*rIydot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+rIxdot*thetadot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+rIxdot*thetadot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+rIydot*thetadot*sin(q1+q2+q3+theta)*(2.0/2.5E1)+rIydot*thetadot*sin(q4+q5+q6+theta)*(2.0/2.5E1)+q1dot*q2dot*cos(q2+q3)*(1.3E1/5.0E2)+q1dot*q3dot*cos(q2+q3)*(1.3E1/5.0E2)+q4dot*q5dot*cos(q5+q6)*(1.3E1/5.0E2)+q4dot*q6dot*cos(q5+q6)*(1.3E1/5.0E2)-q10dot*thetadot*cos(q9+q10)*4.33125E-1+q1dot*thetadot*cos(q2+q3)*(1.3E1/2.5E2)+q2dot*thetadot*cos(q2+q3)*(1.3E1/5.0E2)+q3dot*thetadot*cos(q2+q3)*(1.3E1/5.0E2)+q4dot*thetadot*cos(q5+q6)*(1.3E1/2.5E2)+q5dot*thetadot*cos(q5+q6)*(1.3E1/5.0E2)+q6dot*thetadot*cos(q5+q6)*(1.3E1/5.0E2)-q7dot*thetadot*cos(q7+q8)*4.33125E-1-q8dot*thetadot*cos(q7+q8)*4.33125E-1-q9dot*thetadot*cos(q9+q10)*4.33125E-1-q1dot*rIxdot*cos(q1+theta)*(3.9E1/2.0E1)-q4dot*rIxdot*cos(q4+theta)*(3.9E1/2.0E1)-q7dot*rIxdot*cos(q7+theta)*(9.1E1/8.0E1)-q9dot*rIxdot*cos(q9+theta)*(9.1E1/8.0E1)-rIxdot*thetadot*cos(q1+theta)*(3.9E1/2.0E1)-rIxdot*thetadot*cos(q4+theta)*(3.9E1/2.0E1)-rIxdot*thetadot*cos(q7+theta)*(9.1E1/8.0E1)-rIxdot*thetadot*cos(q9+theta)*(9.1E1/8.0E1)+sqrt(4.1E1)*(q1dot*q1dot)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q2dot*q2dot)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q4dot*q4dot)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(q5dot*q5dot)*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*(thetadot*thetadot)*cos(q6-8.960553845713439E-1)*6.5E-3-q1dot*q2dot*sin(q2+q3)*(1.3E1/5.0E2)-q1dot*q3dot*sin(q2+q3)*(1.3E1/5.0E2)-q4dot*q5dot*sin(q5+q6)*(1.3E1/5.0E2)-q4dot*q6dot*sin(q5+q6)*(1.3E1/5.0E2)-q1dot*thetadot*sin(q2+q3)*(1.3E1/2.5E2)-q2dot*thetadot*sin(q2+q3)*(1.3E1/5.0E2)-q3dot*thetadot*sin(q2+q3)*(1.3E1/5.0E2)-q4dot*thetadot*sin(q5+q6)*(1.3E1/2.5E2)-q5dot*thetadot*sin(q5+q6)*(1.3E1/5.0E2)-q6dot*thetadot*sin(q5+q6)*(1.3E1/5.0E2)+q1dot*rIydot*sin(q1+theta)*(3.9E1/2.0E1)+q4dot*rIydot*sin(q4+theta)*(3.9E1/2.0E1)+q7dot*rIydot*sin(q7+theta)*(9.1E1/8.0E1)+q9dot*rIydot*sin(q9+theta)*(9.1E1/8.0E1)+rIydot*thetadot*sin(q1+theta)*(3.9E1/2.0E1)+rIydot*thetadot*sin(q4+theta)*(3.9E1/2.0E1)+rIydot*thetadot*sin(q7+theta)*(9.1E1/8.0E1)+rIydot*thetadot*sin(q9+theta)*(9.1E1/8.0E1)+sqrt(4.1E1)*q1dot*q2dot*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q1dot*q3dot*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q2dot*q3dot*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q5dot*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q4dot*q6dot*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q5dot*q6dot*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*thetadot*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q2dot*thetadot*cos(q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q3dot*thetadot*cos(q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*thetadot*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q5dot*thetadot*cos(q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q6dot*thetadot*cos(q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*q2dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q1dot*q3dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q2dot*q3dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q4dot*q5dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q4dot*q6dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q5dot*q6dot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q1dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q2dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q3dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q4dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q5dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)+sqrt(4.1E1)*q6dot*thetadot*cos(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q1dot*q2dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q1dot*q3dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q2dot*q3dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q4dot*q5dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q4dot*q6dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q5dot*q6dot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q1dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q2dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q3dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q4dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q5dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q6dot*thetadot*sin(8.960553845713439E-1)*(1.0/5.0E2)-sqrt(4.1E1)*q1dot*rIxdot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*q2dot*rIxdot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*q3dot*rIxdot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*q4dot*rIxdot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*q5dot*rIxdot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*q6dot*rIxdot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*rIxdot*thetadot*cos(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)-sqrt(4.1E1)*rIxdot*thetadot*cos(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q1dot*rIydot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q2dot*rIydot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q3dot*rIydot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q4dot*rIydot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q5dot*rIydot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q6dot*rIydot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*rIydot*thetadot*sin(q1+q2+q3+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*rIydot*thetadot*sin(q4+q5+q6+theta-8.960553845713439E-1)*(1.0/5.0E1)+sqrt(4.1E1)*q1dot*q2dot*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*q3dot*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q5dot*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*q6dot*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q1dot*thetadot*cos(q2+q3-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q2dot*thetadot*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q3dot*thetadot*cos(q2+q3-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q4dot*thetadot*cos(q5+q6-8.960553845713439E-1)*(1.3E1/1.0E3)+sqrt(4.1E1)*q5dot*thetadot*cos(q5+q6-8.960553845713439E-1)*6.5E-3+sqrt(4.1E1)*q6dot*thetadot*cos(q5+q6-8.960553845713439E-1)*6.5E-3;
	return T;
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
void End_Effector_Obs_Dist_Fn(dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,6,1> &End_Effector_Dist, std::vector<int> &End_Effector_Obs)
{
	std::vector<double> r_Pos(2);		double Obs_Dist_i;			int Obs_Dist_Index;
	for (int i = 0; i < 6; i++){
		r_Pos[0] = End_Effector_Pos(2*i);
		r_Pos[1] = End_Effector_Pos(2*i+1);
		Obs_Dist_Fn(r_Pos, Obs_Dist_i, Obs_Dist_Index);
		End_Effector_Dist(i) = Obs_Dist_i;
		End_Effector_Obs[i] = Obs_Dist_Index;}
	return;
}
void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,12,1> &End_Effector_Pos, dlib::matrix<double,12,1> &End_Effector_Vel)
{
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");

	End_Effector_Pos(0) = rA[0];
	End_Effector_Pos(1) = rA[1];
	End_Effector_Pos(2) = rB[0];
	End_Effector_Pos(3) = rB[1];
	End_Effector_Pos(4) = rC[0];
	End_Effector_Pos(5) = rC[1];
	End_Effector_Pos(6) = rD[0];
	End_Effector_Pos(7) = rD[1];
	End_Effector_Pos(8) = rE[0];
	End_Effector_Pos(9) = rE[1];
	End_Effector_Pos(10) = rF[0];
	End_Effector_Pos(11) = rF[1];

	std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");
	std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
	std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");
	std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
	std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
	std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");

	End_Effector_Vel(0) = vA[0];
	End_Effector_Vel(1) = vA[1];
	End_Effector_Vel(2) = vB[0];
	End_Effector_Vel(3) = vB[1];
	End_Effector_Vel(4) = vC[0];
	End_Effector_Vel(5) = vC[1];
	End_Effector_Vel(6) = vD[0];
	End_Effector_Vel(7) = vD[1];
	End_Effector_Vel(8) = vE[0];
	End_Effector_Vel(9) = vE[1];
	End_Effector_Vel(10) = vF[0];
	End_Effector_Vel(11) = vF[1];
	return;
}
void Node_UpdateNCon(Tree_Node &Node_i, Robot_StateNDot &Node_StateNDot_i, std::vector<double> &sigma)
{	// This function is used to substitute the attribute and add it to the All_Nodes
	Node_i.Node_StateNDot = Node_StateNDot_i;
	Node_i.sigma = sigma;
	Node_i.KE = Kinetic_Energy_fn(Node_StateNDot_i);
	Node_i.Node_Index = All_Nodes.size();
	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	End_Effector_PosNVel(Node_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
	Node_i.End_Effector_Pos = End_Effector_Pos;
	Node_i.End_Effector_Vel = End_Effector_Vel;
	All_Nodes.push_back(&Node_i);
	Frontier_Nodes.push_back(&Node_i);
	Frontier_Nodes_Cost.push_back(Node_i.KE);
}
void Nodes_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	// This function is the main optimization constraint function
	double T, KE_i;		dlib::matrix<double> StateNDot_Coeff, Ctrl_Coeff, Contact_Force_Coeff; int Self_Opt_Flag;
	Opt_Seed_Unzip(Opt_Seed, T, StateNDot_Coeff, Ctrl_Coeff, Contact_Force_Coeff);
	std::vector<double> Robot_Config_init(13), Robot_Vel_init(13), Robot_VelfromPos_init(13);
	std::vector<double> Robot_Config_end(13),  Robot_Vel_end(13),  Robot_VelfromPos_end(13);
	dlib::matrix<double> Robot_Acc_init, Robot_Acc_end, Ctrl_init, Contact_Force_init, Ctrl_end, Contact_Force_end, Ctrl_ref, Contact_Force_ref;
	Robot_Acc_init = dlib::zeros_matrix<double>(13,1);			Robot_Acc_end = dlib::zeros_matrix<double>(13,1);
	Ctrl_init = dlib::zeros_matrix<double>(10,1);				Ctrl_end = dlib::zeros_matrix<double>(10,1);
	Contact_Force_init = dlib::zeros_matrix<double>(12,1);		Contact_Force_end = dlib::zeros_matrix<double>(12,1);
	Ctrl_ref = dlib::zeros_matrix<double>(10,1);;				Contact_Force_ref = dlib::zeros_matrix<double>(12,1);

	std::vector<double> Robot_Config_ref, Robot_Velocity_ref;
	Tree_Node Node_i, Node_i_child;
	Node_i = Structure_P.Node_i;		Node_i_child = Structure_P.Node_i_child;
	std::vector<double> Robot_State_init_vec = StateNDot2StateVec(Node_i.Node_StateNDot);
	std::vector<double> sigma_trans, sigma_goal;
	Sigma_TransNGoal(Node_i.sigma, Node_i_child.sigma, sigma_trans, sigma_goal, Self_Opt_Flag);

	std::vector<double> KE_tot; double KE_ref;

	if(Self_Opt_Flag ==1)
	{// In this way, it is the self_optimization`
		KE_ref = 0.1;
	}
	else{
		KE_ref = 0.2 * Node_i.KE;
	}

	// Define the continuity at grids
	for (int i = 0; i < Robot_State_init_vec.size()/2; i++) {
		Robot_Config_ref.push_back(Robot_State_init_vec[i]);
		Robot_Velocity_ref.push_back(Robot_State_init_vec[i+Robot_State_init_vec.size()/2]);
	}
	CtrlNContact_ForcefromCtrlNContact_Force_Coeff(Ctrl_Coeff,Contact_Force_Coeff, 0, 0, Ctrl_ref,  Contact_Force_ref);
	std::vector<double> StateNDot_vec_i;Robot_StateNDot Robot_StateNDot_i;
	dlib::matrix<double> StateNDot_Matrix, Matrix_result, Robot_Vel_Matrix, Robot_VelfromPos_Matrix, Robot_Vel_DlibMatrix, Robot_VelfromPos_DlibMatrix;

	// ObjNConst initialization
	ObjNConstraint_Val.push_back(0);
	ObjNConstraint_Type.push_back(1);

	for (int i = 0; i < Grids-1; i++) {
		Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, i, 0, Robot_Config_init, Robot_Vel_init, Robot_Acc_init, Robot_VelfromPos_init);
		CtrlNContact_ForcefromCtrlNContact_Force_Coeff(Ctrl_Coeff,Contact_Force_Coeff, i, 0, Ctrl_init,  Contact_Force_init);

		Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, i, 1, Robot_Config_end, Robot_Vel_end, Robot_Acc_end, Robot_VelfromPos_end);
		CtrlNContact_ForcefromCtrlNContact_Force_Coeff(Ctrl_Coeff,Contact_Force_Coeff, i, 1, Ctrl_end,  Contact_Force_end);

		// 1. Robotstate, control feasiblity constraint
		StateNDot_vec_i = PosNVel2StateVec(Robot_Config_init, Robot_Vel_init);
		Robot_StateNDot_i = StateVec2StateNDot(StateNDot_vec_i);
		KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
		KE_tot.push_back(KE_i);

		StateNDot_Matrix = StateVec2DlibMatrix_fn(StateNDot_vec_i);
		Matrix_result = StateNDot_Matrix - xlow_vec;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
		Matrix_result = xupp_vec - StateNDot_Matrix;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);

		Matrix_result = Ctrl_init - ctrl_low_vec;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
		Matrix_result = ctrl_upp_vec - Ctrl_init;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 1);
		for (int j = 0; j < Robot_Config_ref.size(); j++) {
			ObjNConstraint_Val.push_back(Robot_Config_init[j] - Robot_Config_ref[j]);
			ObjNConstraint_Type.push_back(0);
			ObjNConstraint_Val.push_back(Robot_Vel_init[j] - Robot_Velocity_ref[j]);
			ObjNConstraint_Type.push_back(0);}
		Matrix_result = Ctrl_init - Ctrl_ref;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

		Matrix_result = Contact_Force_init - Contact_Force_ref;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

		Robot_Vel_DlibMatrix = StateVec2DlibMatrix_fn(Robot_Vel_init);
		Robot_VelfromPos_DlibMatrix = StateVec2DlibMatrix_fn(Robot_VelfromPos_init);
		Matrix_result = Robot_Vel_DlibMatrix - Robot_VelfromPos_DlibMatrix;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

		// reference update
		Robot_Config_ref = Robot_Config_end;		Robot_Velocity_ref = Robot_Vel_end;		Ctrl_ref = Ctrl_end;		Contact_Force_ref = Contact_Force_end;
		Real_ObjNConstraint_Stage(Node_i, Node_i_child, sigma_trans, T, Robot_Config_init, Robot_Vel_init, Robot_Acc_init, Ctrl_init, Contact_Force_init, ObjNConstraint_Val, ObjNConstraint_Type);
		if (i == Grids - 2)
		{
			StateNDot_vec_i = PosNVel2StateVec(Robot_Config_end, Robot_Vel_end);
			Robot_StateNDot_i = StateVec2StateNDot(StateNDot_vec_i);
			KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
			KE_tot.push_back(KE_i);

			Real_ObjNConstraint_Stage(Node_i, Node_i_child,sigma_goal, T, Robot_Config_end, Robot_Vel_end, Robot_Acc_end, Ctrl_end, Contact_Force_end, ObjNConstraint_Val, ObjNConstraint_Type);
		}
	}
	ObjNConstraint_Val[0] = KE_Variation_fn(KE_tot);
	ObjNConstraint_Val.push_back(KE_ref - KE_i);
	ObjNConstraint_Type.push_back(1);
}
double KE_Variation_fn(std::vector<double> &KE_tot)
{	double KE_Variation = 0.0;
	for (int i = 0; i < KE_tot.size()-1; i++)
	{
		KE_Variation = KE_Variation + KE_tot[i] * KE_tot[i];
		// KE_Variation = KE_Variation + (KE_tot[i+1] - KE_tot[i]) * (KE_tot[i+1] - KE_tot[i]);
	}
	return KE_Variation;
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
void Real_ObjNConstraint_Stage(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double>& sigma, double T, std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Acc, dlib::matrix<double> &Ctrl, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	// 2. Dynamics Constraint
	Dynamics_Constraint(Pos, Vel, Acc, Ctrl, Contact_Force, ObjNConstraint_Val, ObjNConstraint_Type);

	// 3. Complementarity constraints: Distance!
	Distance_Velocity_Constraint(sigma, Pos, Vel, ObjNConstraint_Val, ObjNConstraint_Type);

	// 4. Complementarity constraints: Contact Force!
	Contact_Force_Complem_Constraint(Contact_Force, sigma, ObjNConstraint_Val, ObjNConstraint_Type);

	// 5. Contact force feasibility constraints: normal force should be positive and the friction cone constraint has to be satisfied
	Contact_Force_Feasibility_Constraint(Pos, Vel, Contact_Force, ObjNConstraint_Val, ObjNConstraint_Type);

	// 6. Contact maintenance constraint: the previous unchanged active constraint have to be satisfied
	Contact_Maintenance_Constraint(Node_i, Node_i_child, Pos, Vel, ObjNConstraint_Val, ObjNConstraint_Type);

}
void Contact_Maintenance_Constraint(Tree_Node &Node_i, Tree_Node &Node_i_child, std::vector<double> &Pos, std::vector<double> &Vel, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	std::vector<double> sigma_i = Node_i.sigma;
	std::vector<double> sigma_i_child = Node_i_child.sigma;

	dlib::matrix<double,12,1> End_Effector_Pos_ref, End_Effector_Vel_ref;
	End_Effector_Pos_ref = Node_i.End_Effector_Pos;
	End_Effector_Vel_ref = Node_i.End_Effector_Vel;

	std::vector<double> StateVec_i;	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	StateVec_i = PosNVel2StateVec(Pos, Vel);
	Robot_StateNDot StateNDot_Init_i(StateVec_i);
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);

	std::vector<double> sigma_maint(12);
	sigma_maint[0] = sigma_i[0] * sigma_i_child[0];		sigma_maint[1] = sigma_i[0] * sigma_i_child[0];
	sigma_maint[2] = sigma_i[0] * sigma_i_child[0];		sigma_maint[3] = sigma_i[0] * sigma_i_child[0];
	sigma_maint[4] = sigma_i[1] * sigma_i_child[1];		sigma_maint[5] = sigma_i[1] * sigma_i_child[1];
	sigma_maint[6] = sigma_i[1] * sigma_i_child[1];		sigma_maint[7] = sigma_i[1] * sigma_i_child[1];
	sigma_maint[8] = sigma_i[2] * sigma_i_child[2];		sigma_maint[9] = sigma_i[2] * sigma_i_child[2];
	sigma_maint[10] = sigma_i[3] * sigma_i_child[3];	sigma_maint[11] = sigma_i[3] * sigma_i_child[3];

	dlib::matrix<double> Maint_Matrix = Diag_Matrix_fn(sigma_maint);
	dlib::matrix<double> Matrix_result;

	Matrix_result = Maint_Matrix * (End_Effector_Pos - End_Effector_Pos_ref);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	Matrix_result = Maint_Matrix * (End_Effector_Vel - End_Effector_Vel_ref);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
}
void Contact_Force_Feasibility_Constraint(std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	std::vector<double> StateVec_i;	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	StateVec_i = PosNVel2StateVec(Pos, Vel);
	Robot_StateNDot StateNDot_Init_i(StateVec_i);
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);
	dlib::matrix<double,6,1> End_Effector_Dist; std::vector<int> End_Effector_Obs(6);
	End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);
	double Contact_Force_i_x, Contact_Force_i_y; dlib::matrix<double> Matrix_result;
	std:vector<double> Normal_Force, Tange_Force;
	for (int i = 0; i < Contact_Force.nr()/2; i++) {
		Contact_Force_i_x = Contact_Force(2*i);
		Contact_Force_i_y = Contact_Force(2*i+1);
		Normal_Force.push_back(Contact_Force_i_x * Envi_Map_Normal(End_Effector_Obs[i],0) + Contact_Force_i_y * Envi_Map_Normal(End_Effector_Obs[i],1));
		Tange_Force.push_back(Contact_Force_i_x * Envi_Map_Tange(End_Effector_Obs[i],0) + Contact_Force_i_y * Envi_Map_Tange(End_Effector_Obs[i],1));}

	for (int i = 0; i < Contact_Force.nr()/2; i++) {
		ObjNConstraint_Val.push_back(Normal_Force[i]);
		ObjNConstraint_Type.push_back(1);}

	double Normal_Force_1, Normal_Force_2, Normal_Force_3, Normal_Force_4;		double Tange_Force_1, Tange_Force_2, Tange_Force_3, Tange_Force_4;
	Normal_Force_1 = Normal_Force[0] + Normal_Force[1];		Normal_Force_2 = Normal_Force[2] + Normal_Force[3];		Normal_Force_3 = Normal_Force[4];	Normal_Force_4 = Normal_Force[5];
	Tange_Force_1 = Tange_Force[0] + Tange_Force[1];		Tange_Force_2 = Tange_Force[2] + Tange_Force[3];		Tange_Force_3 = Tange_Force[4];		Tange_Force_4 = Tange_Force[5];

	ObjNConstraint_Val.push_back(Normal_Force_1 * Normal_Force_1 * mu * mu - Tange_Force_1 * Tange_Force_1);
	ObjNConstraint_Type.push_back(1);
	ObjNConstraint_Val.push_back(Normal_Force_2 * Normal_Force_2 * mu * mu - Tange_Force_2 * Tange_Force_2);
	ObjNConstraint_Type.push_back(1);
	ObjNConstraint_Val.push_back(Normal_Force_3 * Normal_Force_3 * mu * mu - Tange_Force_3 * Tange_Force_3);
	ObjNConstraint_Type.push_back(1);
	ObjNConstraint_Val.push_back(Normal_Force_4 * Normal_Force_4 * mu * mu - Tange_Force_4 * Tange_Force_4);
	ObjNConstraint_Type.push_back(1);
	return;
}
void Contact_Force_Complem_Constraint(dlib::matrix<double> &Contact_Force, std::vector<double> &sigma, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	dlib::matrix<double> Contact_Force_Complem_Matrix, Matrix_result;
	std::vector<double> sigma_temp;
	sigma_temp.push_back(!sigma[0]);
	sigma_temp.push_back(!sigma[0]);
	sigma_temp.push_back(!sigma[0]);
	sigma_temp.push_back(!sigma[0]);
	sigma_temp.push_back(!sigma[1]);
	sigma_temp.push_back(!sigma[1]);
	sigma_temp.push_back(!sigma[1]);
	sigma_temp.push_back(!sigma[1]);
	sigma_temp.push_back(!sigma[2]);
	sigma_temp.push_back(!sigma[2]);
	sigma_temp.push_back(!sigma[3]);
	sigma_temp.push_back(!sigma[3]);
	Contact_Force_Complem_Matrix = Diag_Matrix_fn(sigma_temp);
	Matrix_result = Contact_Force_Complem_Matrix * Contact_Force;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	return;
}
void Distance_Velocity_Constraint(std::vector<double>& sigma, std::vector<double> &Pos, std::vector<double> &Vel, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	std::vector<double> Robotstate_Vec_i = PosNVel2StateVec(Pos, Vel);
	Robot_StateNDot StateNDot_Init_i(Robotstate_Vec_i);
	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);

	dlib::matrix<double,6,1> End_Effector_Dist;
	std::vector<int> End_Effector_Obs(6);

	End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

	dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Eqn_Maint_Matrix, Matrix_result;
	std::vector<double> sigma_temp;
	sigma_temp = Sigma2Pos(sigma, 0);			Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Pos(sigma, 1);			Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Vel(sigma);				Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

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

}
void Dynamics_Constraint(std::vector<double> &Pos, std::vector<double> &Vel, dlib::matrix<double> &Acc, dlib::matrix<double> &Ctrl, dlib::matrix<double> &Contact_Force, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	std::vector<double> Robotstate_Vec_i = PosNVel2StateVec(Pos, Vel);
	Robot_StateNDot Robot_StateNDot_i = StateVec2StateNDot(Robotstate_Vec_i);
	dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full, Jac_Full_Trans;
	Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
	Jac_Full_Trans = dlib::trans(Jac_Full);
	dlib::matrix<double> Matrix_result;
	Matrix_result = D_q * Acc + C_q_qdot - Jac_Full_Trans * Contact_Force - B_q * Ctrl;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
}
dlib::matrix<double> StateVec2DlibMatrix_fn(std::vector<double> &StateVec)
{	const int dim = StateVec.size();
	dlib::matrix<double> StateVec2DlibMatrix;
	StateVec2DlibMatrix = dlib::zeros_matrix<double>(dim,1);
	for (int i = 0; i < dim; i++) {
		StateVec2DlibMatrix(i) = StateVec[i];}
	return StateVec2DlibMatrix;
}
void CtrlNContact_ForcefromCtrlNContact_Force_Coeff(dlib::matrix<double> &Ctrl_Coeff,dlib::matrix<double> &Contact_Force_Coeff, int Grid_Ind, double s, dlib::matrix<double> &Ctrl_i,  dlib::matrix<double> &Contact_Force_i)
{	double x_a, x_b;
	for (int i = 0; i < 10; i++) {
		x_a = Ctrl_Coeff(2*i, Grid_Ind);
		x_b = Ctrl_Coeff(2*i+1, Grid_Ind);
		Ctrl_i(i) = x_a * s + x_b;}
	for (int i = 0; i < 12; i++) {
		x_a = Contact_Force_Coeff(2*i, Grid_Ind);
		x_b = Contact_Force_Coeff(2*i+1, Grid_Ind);
		Contact_Force_i(i) = x_a * s + x_b;}
}
std::vector<double> Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, int &Nodes_Opt_Flag)
{	// This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
	// However, the constraint will be set to use the direct collocation method
	Structure_P.Node_i = Node_i;		Structure_P.Node_i_child = Node_i_child;
	std::vector<double> Opt_Seed = Seed_Guess_Gene(Node_i, Node_i_child);
	std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
	Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	snoptProblem Nodes_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Opt_Seed.size();
	integer neF = ObjNConstraint_Val.size();
	integer lenA  =  n * neF;

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

	for (int i = 0; i < n; i++) {
		xlow[i] =-Inf;
		xupp[i] = Inf;
		xstate[i] = 0.0;
		x[i] = Opt_Seed[i];  	// Initial guess
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
	Nodes_Optimization_Pr.setPrintFile  ( "Nodes_Optimization_Pr.out" );
	Nodes_Optimization_Pr.setProblemSize( n, neF );
	Nodes_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
	Nodes_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
	Nodes_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
	Nodes_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Nodes_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Nodes_Optimization_Pr.setXNames     ( xnames, nxnames );
	Nodes_Optimization_Pr.setFNames     ( Fnames, nFnames );
	Nodes_Optimization_Pr.setProbName   ( "Nodes_Optimization_Pr" );
	Nodes_Optimization_Pr.setUserFun    ( Nodes_Optimization_Pr_fn);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Nodes_Optimization_Pr.computeJac    ();
	Nodes_Optimization_Pr.setIntParameter( "Derivative option", 0 );
	integer Cold = 0, Basis = 1, Warm = 2;
	Nodes_Optimization_Pr.solve( Cold );

	for (int i = 0; i < n; i++)
	{
		Opt_Seed[i] = x[i];
	}

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return Opt_Seed;
}
int Nodes_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru )
{	 std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
	for (int i = 0; i < *n; i++){
		Opt_Seed.push_back(x[i]);}
		Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	for (int i = 0; i < ObjNConstraint_Val.size(); i++){
		F[i] = ObjNConstraint_Val[i];}
	return 0;
}

std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child)
{	// This function will generate the spline coefficients needed for the further optimization
	double T = 0.5;
	// The first step is to generate a feasible configuration that can satisfy the contact mode in the node i child
	std::vector<double> Init_Config = StateNDot2StateVec(Node_i.Node_StateNDot);
	std::vector<double> Seed_Config = Seed_Guess_Gene_Robotstate(Node_i, Node_i_child);
	const int StateNDot_len = Init_Config.size();			const int Control_len = 10;			const int Contact_Force_len = 12;

	dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
	StateNDot_Traj = dlib::zeros_matrix<double>(StateNDot_len, Grids);				// each variable traj will be in a row fashion
	Ctrl_Traj = dlib::zeros_matrix<double>(Control_len, Grids);
	Contact_Force_Traj = dlib::zeros_matrix<double>(Contact_Force_len, Grids);

	dlib::matrix<double> StateNDot_Coeff, Ctrl_Coeff, Contact_Force_Coeff;
	StateNDot_Coeff = dlib::zeros_matrix<double>(4*StateNDot_len, Grids-1);
	Ctrl_Coeff = dlib::zeros_matrix<double>(2*Control_len, Grids-1);
	Contact_Force_Coeff = dlib::zeros_matrix<double>(2*Contact_Force_len, Grids-1);

	dlib::matrix<double> Robot_State_Interpol_i;
	for (int i = 0; i < StateNDot_len; i++){
		Robot_State_Interpol_i = dlib::linspace(Init_Config[i], Seed_Config[i], Grids);
		for (int j = 0; j < Grids; j++){
			StateNDot_Traj(i,j) = Robot_State_Interpol_i(j);}}

	// Here is the initialization of the spline coefficients for the stateNdot, control and contact force
	// First is to initialize: StateNDot
	double x_init, x_end, xdot_init, xdot_end;
	std::vector<double> CubicSpline_Coeff_Pos_vec, CubicSpline_Coeff_Vel_vec, PVA_init, PVA_end;
	for (int i = 0; i < Grids-1; i++)
	{	for (int j = 0; j < StateNDot_len/2; j++){
		    // Position first
			x_init = StateNDot_Traj(j,i);
			x_end = StateNDot_Traj(j,i+1);
			xdot_init = StateNDot_Traj(j+StateNDot_len/2,i);
			xdot_end = StateNDot_Traj(j+StateNDot_len/2,i+1);
			CubicSpline_Coeff_Pos_vec = CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end); // 4 by 1 vector: a, b, c, d
			// Velocity second
			PVA_init = CubicSpline_PosVelAcc4(T, CubicSpline_Coeff_Pos_vec[0], CubicSpline_Coeff_Pos_vec[1], CubicSpline_Coeff_Pos_vec[2], CubicSpline_Coeff_Pos_vec[3], 0);
			PVA_end  = CubicSpline_PosVelAcc4(T, CubicSpline_Coeff_Pos_vec[0], CubicSpline_Coeff_Pos_vec[1], CubicSpline_Coeff_Pos_vec[2], CubicSpline_Coeff_Pos_vec[3], 1);
			CubicSpline_Coeff_Vel_vec = CubicSpline_Coeff_fn(T, xdot_init,xdot_end, PVA_init[2], PVA_end[2]); // 4 by 1 vector: a, b, c, d
			StateNDot_Coeff(8*j,i)   = CubicSpline_Coeff_Pos_vec[0];
			StateNDot_Coeff(8*j+1,i) = CubicSpline_Coeff_Pos_vec[1];
			StateNDot_Coeff(8*j+2,i) = CubicSpline_Coeff_Pos_vec[2];
			StateNDot_Coeff(8*j+3,i) = CubicSpline_Coeff_Pos_vec[3];
			StateNDot_Coeff(8*j+4,i) = CubicSpline_Coeff_Vel_vec[0];
			StateNDot_Coeff(8*j+5,i) = CubicSpline_Coeff_Vel_vec[1];
			StateNDot_Coeff(8*j+6,i) = CubicSpline_Coeff_Vel_vec[2];
			StateNDot_Coeff(8*j+7,i) = CubicSpline_Coeff_Vel_vec[3];
		}
	}
	// cout<<StateNDot_Coeff<<endl;

	// Second is to initialize: Control and Contact Force
	std::vector<double> Robot_Pos(13), Robot_Vel(13), Robot_VelfromPos(13), Robotstate_Vec_i;		dlib::matrix<double> Robot_Acc; Robot_Acc = dlib::zeros_matrix<double>(13,1);
	dlib::matrix<double> D_q, B_q, C_q_qdot, Jac_Full, Jac_Full_Trans, Dynamics_LHS, Dynamics_RHS, Dynamics_RHS_Matrix;		Robot_StateNDot Robot_StateNDot_i;
	for (int i = 0; i < Grids-1; i++){
		Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, i, 0, Robot_Pos, Robot_Vel, Robot_Acc, Robot_VelfromPos);
		Robotstate_Vec_i = PosNVel2StateVec(Robot_Pos, Robot_Vel);
		Robot_StateNDot_i = StateVec2StateNDot(Robotstate_Vec_i);
		Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
		Dynamics_LHS = D_q * Robot_Acc + C_q_qdot;
		Dynamics_RHS_Matrix = Dynamics_RHS_Matrix_fn(Jac_Full, B_q);
		Dynamics_RHS = dlib::pinv(Dynamics_RHS_Matrix) * Dynamics_LHS;
		for (int j = 0; j < Dynamics_RHS.nr(); j++) {
			if(j<Contact_Force_Traj.nr())
			{	Contact_Force_Traj(j,i) = Dynamics_RHS(j);}
			else
			{
				Ctrl_Traj(j - Contact_Force_Traj.nr(),i) = Dynamics_RHS(j);
			}
		}
	}
	// cout<<StateNDot_Traj<<endl;			cout<<Contact_Force_Traj<<endl;			cout<<Ctrl_Traj<<endl;			cout<<StateNDot_Coeff<<endl;
	// THen the last column
	Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(T, StateNDot_Coeff, Grids-2, 1, Robot_Pos, Robot_Vel, Robot_Acc, Robot_VelfromPos);
	Robotstate_Vec_i = PosNVel2StateVec(Robot_Pos, Robot_Vel);
	Robot_StateNDot_i = StateVec2StateNDot(Robotstate_Vec_i);
	Dynamics_Matrices(Robot_StateNDot_i, D_q, B_q, C_q_qdot, Jac_Full);
	Dynamics_LHS = D_q * Robot_Acc + C_q_qdot;
	Dynamics_RHS_Matrix = Dynamics_RHS_Matrix_fn(Jac_Full, B_q);
	Dynamics_RHS = dlib::pinv(Dynamics_RHS_Matrix) * Dynamics_LHS;
	for (int j = 0; j < Dynamics_RHS.nr(); j++) {
		if(j<Contact_Force_Traj.nr())
		{	Contact_Force_Traj(j,Grids-1) = Dynamics_RHS(j);}
		else
		{
			Ctrl_Traj(j - Contact_Force_Traj.nr(),Grids-1) = Dynamics_RHS(j);
		}
	}
	// cout<<StateNDot_Traj<<endl;			cout<<Contact_Force_Traj<<endl;			cout<<Ctrl_Traj<<endl;
	// The calculation of the coefficients of the control and contact force is easier compared to the robot state due to the assumption of the linear equation
	Ctrl_Contact_Force_Coeff_fn(Ctrl_Traj, Contact_Force_Traj, Ctrl_Coeff, Contact_Force_Coeff);

	// The final task is to pile them into a single vector
	std::vector<double> Opt_Seed;
	Opt_Seed.push_back(T);
	Opt_Seed_Zip(Opt_Seed, StateNDot_Coeff, Ctrl_Coeff, Contact_Force_Coeff);

	// cout<<StateNDot_Coeff<<endl;
	// cout<<Ctrl_Coeff<<endl;
	// cout<<Contact_Force_Coeff<<endl;
	//
	return Opt_Seed;
}
void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T, dlib::matrix<double> & StateNDot_Coeff, dlib::matrix<double> & Ctrl_Coeff, dlib::matrix<double> & Contact_Force_Coeff)
{
	T = Opt_Seed[0];
	const int NumOfStateNDot_Coeff= 4 * 26;
	const int NumOfCtrl_Coeff= 2 * 10;
	const int NumOfContactForce_Coeff= 2 * 12;
	int Opt_Seed_Index = 1;
	StateNDot_Coeff = dlib::zeros_matrix<double>(NumOfStateNDot_Coeff, Grids-1);
	Ctrl_Coeff = dlib::zeros_matrix<double>(NumOfCtrl_Coeff, Grids-1);
	Contact_Force_Coeff = dlib::zeros_matrix<double>(NumOfContactForce_Coeff, Grids-1);
	// cout<<"To be compared!"<<endl;
	// 1. Retrieve the StateNDot_Coeff matrix
	for (int i = 0; i < Grids-1; i++) {
		for (int j = 0; j < NumOfStateNDot_Coeff; j++) {
			StateNDot_Coeff(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<StateNDot_Coeff<<endl;
	// 2. Retrieve the control matrix
	for (int i = 0; i < Grids-1; i++) {
		for (int j = 0; j < NumOfCtrl_Coeff; j++) {
			Ctrl_Coeff(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Ctrl_Coeff<<endl;
	// 3. Retrieve the contact force matrix
	for (int i = 0; i < Grids-1; i++) {
		for (int j = 0; j < NumOfContactForce_Coeff; j++) {
			Contact_Force_Coeff(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Contact_Force_Coeff<<endl;
}
void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Coeff, dlib::matrix<double> & Ctrl_Coeff, dlib::matrix<double> & Contact_Force_Coeff)
{	// This function is used to stack the coefficient matrices into a column vector
	for (int i = 0; i < StateNDot_Coeff.nc(); i++) {
		for (int j = 0; j < StateNDot_Coeff.nr(); j++) {
			Opt_Seed.push_back(StateNDot_Coeff(j,i));
		}
	}
	for (int i = 0; i < Ctrl_Coeff.nc(); i++) {
		for (int j = 0; j < Ctrl_Coeff.nr(); j++) {
			Opt_Seed.push_back(Ctrl_Coeff(j,i));
		}
	}
	for (int i = 0; i < Contact_Force_Coeff.nc(); i++) {
		for (int j = 0; j < Contact_Force_Coeff.nr(); j++) {
			Opt_Seed.push_back(Contact_Force_Coeff(j,i));
		}
	}
}
void Ctrl_Contact_Force_Coeff_fn(dlib::matrix<double> &Ctrl_Traj, dlib::matrix<double> &Contact_Force_Traj, dlib::matrix<double> &Ctrl_Coeff, dlib::matrix<double> &Contact_Force_Coeff)
{
	// y = a * s + b
	double Ctrl_init, Ctrl_end, Contact_Force_init, Contact_Force_end;
	for (int i = 0; i < Grids-1; i++){
		// Computation of the control coeff

		for (int j = 0; j < 10; j++) {
			Ctrl_init = Ctrl_Traj(j,i);
			Ctrl_end = Ctrl_Traj(j,i+1);
			Ctrl_Coeff(2*j,i) = Ctrl_end - Ctrl_init;
			Ctrl_Coeff(2*j+1,i) = Ctrl_init;
		}
		// Computation of the contact force coeff
		for (int j = 0; j < 12; j++) {
			Contact_Force_init = Contact_Force_Traj(j,i);
			Contact_Force_end = Contact_Force_Traj(j,i+1);
			Contact_Force_Coeff(2*j,i) = Contact_Force_end - Contact_Force_init;
			Contact_Force_Coeff(2*j+1,i) = Contact_Force_init;
		}
	}
	return;
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
std::vector<double> PosNVel2StateVec(std::vector<double> & Pos, std::vector<double> & Vel)
{
	std::vector<double> StateVec_i;
	for (int i = 0; i < Pos.size(); i++) {
		StateVec_i.push_back(Pos[i]);
	}
	for (int i = 0; i < Vel.size(); i++) {
		StateVec_i.push_back(Vel[i]);
	}
	return StateVec_i;
}
void Dynamics_Matrices(Robot_StateNDot &Node_StateNDot, dlib::matrix<double> &D_q, dlib::matrix<double> &B_q, dlib::matrix<double> &C_q_qdot, dlib::matrix<double> &Jac_Full)
{
	D_q = D_q_fn(Node_StateNDot);
	B_q = B_q_fn();
	C_q_qdot = C_q_qdot_fn(Node_StateNDot);
	Jac_Full = Jac_Full_fn(Node_StateNDot);
}
void Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(double T, dlib::matrix<double> &StateNDot_Coeff, int Grid_Ind, double s, std::vector<double> &Robot_Config,  std::vector<double> &Robot_Vel, dlib::matrix<double> &Robot_Acc, std::vector<double> &Robot_VelfromPos)
{
	//Here Grid_Ind denotes which Grid we are talking about and s is a value between 0 and 1
	std::vector<double> PVAVP_i;
	double x_a, x_b, x_c, x_d, xdot_a, xdot_b, xdot_c, xdot_d;
	for (int i = 0; i < 13; i++) {
		x_a = StateNDot_Coeff(8*i, Grid_Ind);
		x_b = StateNDot_Coeff(8*i+1, Grid_Ind);
		x_c = StateNDot_Coeff(8*i+2, Grid_Ind);
		x_d = StateNDot_Coeff(8*i+3, Grid_Ind);

		xdot_a = StateNDot_Coeff(8*i+4, Grid_Ind);
		xdot_b = StateNDot_Coeff(8*i+5, Grid_Ind);
		xdot_c = StateNDot_Coeff(8*i+6, Grid_Ind);
		xdot_d = StateNDot_Coeff(8*i+7, Grid_Ind);

		PVAVP_i = CubicSpline_PosVelAcc8(T, x_a, x_b, x_c, x_d, xdot_a, xdot_b, xdot_c, xdot_d, s);
		Robot_Config[i] = PVAVP_i[0];
		Robot_Vel[i] = PVAVP_i[1];
		Robot_Acc(i) = PVAVP_i[2];
		Robot_VelfromPos[i] = PVAVP_i[3];
	}
	return;
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
std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
	// This function is used to generate a configuration to initialize the optimization
	std::vector<double> Robot_State_Seed = StateNDot2StateVec(Node_i.Node_StateNDot);
	std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
	Seed_Conf_Optimization_ObjNConstraint(Robot_State_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	snoptProblem Seed_Conf_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Robot_State_Seed.size();
	integer neF = ObjNConstraint_Val.size();
	integer lenA  =  n * neF;

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

	for (int i = 0; i < n; i++) {
		xlow[i] = xlow_vec(i);
		xupp[i] = xupp_vec(i);
		xstate[i] = 0.0;
		x[i] = Robot_State_Seed[i];  	// Initial guess
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
	Seed_Conf_Optimization_Pr.setPrintFile  ( "Seed_Conf_Optimization_Pr.out" );
	Seed_Conf_Optimization_Pr.setProblemSize( n, neF );
	Seed_Conf_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
	Seed_Conf_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
	Seed_Conf_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
	Seed_Conf_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Seed_Conf_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Seed_Conf_Optimization_Pr.setXNames     ( xnames, nxnames );
	Seed_Conf_Optimization_Pr.setFNames     ( Fnames, nFnames );
	Seed_Conf_Optimization_Pr.setProbName   ( "Seed_Conf_Optimization_Pr" );
	Seed_Conf_Optimization_Pr.setUserFun    ( Seed_Conf_Optimization_Pr_fn_);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Seed_Conf_Optimization_Pr.computeJac    ();
	Seed_Conf_Optimization_Pr.setIntParameter( "Derivative option", 0 );
	integer Cold = 0, Basis = 1, Warm = 2;
	Seed_Conf_Optimization_Pr.solve( Cold );

	for (int i = 0; i < 26; i++)
	{
		Robot_State_Seed[i] = x[i];
	}

	Robot_StateNDot Robot_StateNDot_Seed(Robot_State_Seed);
	Robot_Plot_fn(Robot_StateNDot_Seed);

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return Robot_State_Seed;
}
int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru )
{	 std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
	for (int i = 0; i < 26; i++){
		Opt_Seed.push_back(x[i]);}
	Seed_Conf_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	for (int i = 0; i < ObjNConstraint_Val.size(); i++){
		F[i] = ObjNConstraint_Val[i];}
	return 0;
}
void Seed_Conf_Optimization_ObjNConstraint(std::vector<double> &Opt_Seed, std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	Robot_StateNDot StateNDot_Init_i(Opt_Seed);		dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel, End_Effector_Pos_ref, End_Effector_Vel_ref;
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);

	End_Effector_Pos_ref = Structure_P.Node_i.End_Effector_Pos;
	End_Effector_Vel_ref = Structure_P.Node_i.End_Effector_Vel;

	double KE_opt = Kinetic_Energy_fn(StateNDot_Init_i);
	ObjNConstraint_Val.push_back(KE_opt);
	ObjNConstraint_Type.push_back(1);

	std::vector<double> sigma_i = Structure_P.Node_i.sigma;
	std::vector<double> sigma_i_child = Structure_P.Node_i_child.sigma;
	dlib::matrix<double,6,1> End_Effector_Dist;
	std::vector<int> End_Effector_Obs(6);

	End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

	dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix, Eqn_Maint_Matrix, Matrix_result;
	std::vector<double> sigma_temp;
	sigma_temp = Sigma2Pos(sigma_i_child, 0);			Eqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Pos(sigma_i_child, 1);			Ineqn_Pos_Matrix = Diag_Matrix_fn(sigma_temp);
	sigma_temp = Sigma2Vel(sigma_i_child);				Eqn_Vel_Matrix = Diag_Matrix_fn(sigma_temp);

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

	// 4. One more constraint to be added is to maintain the active unchanged constraint
	Eqn_Maint_Matrix = Eqn_Maint_Matrix_fn(sigma_i, sigma_i_child);
	Matrix_result = Eqn_Maint_Matrix * (End_Effector_Pos_ref - End_Effector_Pos);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	return;
}
dlib::matrix<double> Eqn_Maint_Matrix_fn(std::vector<double> &sigma_i, std::vector<double> &sigma_i_child)
{
	// This function is used to generate the contact maintenance matrix
	std::vector<double> sigma_maint;
	sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
	sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
	sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
	sigma_maint.push_back(sigma_i[0] * sigma_i_child[0]);
	sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
	sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
	sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
	sigma_maint.push_back(sigma_i[1] * sigma_i_child[1]);
	sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
	sigma_maint.push_back(sigma_i[2] * sigma_i_child[2]);
	sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);
	sigma_maint.push_back(sigma_i[3] * sigma_i_child[3]);

	dlib::matrix<double> Eqn_Maint_Matrix;
	Eqn_Maint_Matrix = Diag_Matrix_fn(sigma_maint);
	return Eqn_Maint_Matrix;
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
