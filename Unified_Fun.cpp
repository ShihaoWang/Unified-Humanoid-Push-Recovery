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
double Inf = 1.1e20;

double PI = 3.1415926535897932384626;
double pi = PI;

double rIxlow = -Inf;                  double rIxupp = Inf;
double rIylow = -Inf;                  double rIyupp = Inf;
double thetalow = -pi;                 double thetaupp = pi;
double q1low = -2.1817;                double q1upp = 0.733;
double q2low = 0.0;                    double q2upp = 2.618;
double q3low = -1.3;                   double q3upp = 0.733;
double q4low = -2.1817;                double q4upp = 0.733;
double q5low = 0.0;                    double q5upp = 2.618;
double q6low = -1.3;                   double q6upp = 0.733;
double q7low = -3.14;                  double q7upp = 1.047;
double q8low = -2.391;                 double q8upp = 0.0;
double q9low = -3.14;                  double q9upp = 1.047;
double q10low = -2.391;                double q10upp = 0.0;

double AngRateLow = -3.0;              double AngRateHgh = 3.0;

double rIxdotlow = -Inf;               double rIxdotupp = Inf;
double rIydotlow = -Inf;               double rIydotupp = Inf;
double thetadotlow = -Inf;             double thetadotupp = Inf;
double q1dotlow = AngRateLow;          double q1dotupp = AngRateHgh;
double q2dotlow = AngRateLow;          double q2dotupp = AngRateHgh;
double q3dotlow = AngRateLow;          double q3dotupp = AngRateHgh;
double q4dotlow = AngRateLow;          double q4dotupp = AngRateHgh;
double q5dotlow = AngRateLow;          double q5dotupp = AngRateHgh;
double q6dotlow = AngRateLow;          double q6dotupp = AngRateHgh;
double q7dotlow = AngRateLow;          double q7dotupp = AngRateHgh;
double q8dotlow = AngRateLow;          double q8dotupp = AngRateHgh;
double q9dotlow = AngRateLow;          double q9dotupp = AngRateHgh;
double q10dotlow = AngRateLow;         double q10dotupp = AngRateHgh;

double tau1_max = 100;             		double tau2_max = 100;
double tau3_max = 100;					double tau4_max = 100;
double tau5_max = 100;             		double tau6_max = 100;
double tau7_max = 60;              		double tau8_max = 50;
double tau9_max = 60;             		double tau10_max = 50;

dlib::matrix<double,26,1> xlow_vec;
dlib::matrix<double,26,1> xupp_vec;
dlib::matrix<double,10,1> ctrl_low_vec;
dlib::matrix<double,10,1> ctrl_upp_vec;

double Step_Length_Max = 0.65 * (abs(sin(q1low)) + abs(sin(q1upp)));
int Nodes_Tot_Number;		dlib::matrix<double> Envi_Map;
double mini = 0.05;

std::vector<Tree_Node_Ptr> All_Nodes;				// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		// The kinetic energy of each nodes

void Add_Node(Tree_Node &Current_Node)
{
	// The nodes can only be added if
	All_Nodes.push_back(&Current_Node);
	Frontier_Nodes.push_back(&Current_Node);
	Frontier_Nodes_Cost.push_back(Kinetic_Energy_fn(Current_Node.StateNDot_Str));
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

int Add_Child2Par(Tree_Node &Child_Node, Tree_Node &Par_Node, int &Nodes_Tot_Number, std::vector<double> &sigma_i)
{

    Child_Node.Par_Node = &Par_Node;

    Child_Node.Node_Number = Nodes_Tot_Number + 1;

	Child_Node.sigma_i = sigma_i;

	Children_Nodes.push_back(&Child_Node);

    Nodes_Tot_Number = Nodes_Tot_Number + 1;

	Par_Node.Children_Nodes.push_back(&Child_Node);

    return 0;
}

std::vector<double> Default_Init(const std::vector<double> &sigma_i, Unified_Structure_P &P, int Flag)
{

	Envi_Map = dlib::ones_matrix<double>(2,4);
	Envi_Map(0,0) = -100.0;		Envi_Map(0,1) = 0.0;		Envi_Map(0,2) = 100.0;		Envi_Map(0,3) = 0.0;
	// Obs2: Vertical wall at some distance
	Envi_Map(1,0) = 5.0;			Envi_Map(1,1) = 0.0;		Envi_Map(1,2) = 5.0;		Envi_Map(1,3) = 10.0;

	// This function is used to initialize a given configuration
    Nodes_Tot_Number = 0;
	vector<double> Robot_State_Init;

	ifstream Initial_Robot_State_File;              // This is to read the initial angle and angular velocities
	Initial_Robot_State_File.open("init_robot_state.txt");
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
		printf("Unable to open file!\n");
	}

	if(Flag ==1)// This means that the given initial condition works for the constraint
	{
		return Robot_State_Init;
	}
	else
	{
		// Robot_StateNDot Robot_StateNDot_init(Robot_State_Init);
		// std::string input_name = "init_given";
		// Robot_Plot_fn(Robot_StateNDot_init,input_name);

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

		snoptProblem Default_Init_Pr;                     // This is the name of the Optimization problem
		// Allocate and initialize
		integer n = 14;
		integer neF = 17;     // 1 objective function
		integer lenA  =  n * neF;                              // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

		integer lenru = 14;             					  // This is used to pass the initial state into the usrfun
		doublereal *ru = new doublereal[lenru];
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
		// First set the lower and upper bounds for state
		xlow[0] = rIylow;		xupp[0] = rIyupp;
		xlow[1] = q1low;		xupp[1] = q1upp;
		xlow[2] = q2low;		xupp[2] = q2upp;
		xlow[3] = q3low;		xupp[3] = q3upp;
		xlow[4] = q4low;		xupp[4] = q4upp;
		xlow[5] = q5low;		xupp[5] = q5upp;
		xlow[6] = q6low;		xupp[6] = q6upp;
		xlow[7] = rIydotlow;	xupp[7] = rIydotupp;
		xlow[8] = q1dotlow;		xupp[8] = q1dotupp;
		xlow[9] = q2dotlow;		xupp[9] = q2dotupp;
		xlow[10] = q3dotlow;	xupp[10] = q3dotupp;
		xlow[11] = q4dotlow;	xupp[11] = q4dotupp;
		xlow[12] = q5dotlow;	xupp[12] = q5dotupp;
		xlow[13] = q6dotlow;	xupp[13] = q6dotupp;

		// Second set the lower and upper bounds for the objective function
		Flow[0] = 0;			Fupp[0] = Inf;		// This is the vertical position of the robot COM
		Flow[1] = 0;			Fupp[1] = 0;
		Flow[2] = 0;			Fupp[2] = 0;
		Flow[3] = 0;			Fupp[3] = 0;
		Flow[4] = 0;			Fupp[4] = 0;
		Flow[5] = 0;			Fupp[5] = 0;
		Flow[6] = 0;			Fupp[6] = 0;
		Flow[7] = 0;			Fupp[7] = 5.0;
		Flow[8] = 0;			Fupp[8] = 5.0;
		Flow[1+8] = 0;			Fupp[1+8] = 0;
		Flow[2+8] = 0;			Fupp[2+8] = 0;
		Flow[3+8] = 0;			Fupp[3+8] = 0;
		Flow[4+8] = 0;			Fupp[4+8] = 0;
		Flow[5+8] = 0;			Fupp[5+8] = 0;
		Flow[6+8] = 0;			Fupp[6+8] = 0;
		Flow[7+8] = 0;			Fupp[7+8] = 5.0;
		Flow[8+8] = 0;			Fupp[8+8] = 5.0;

		x[0] = Structure_P.rIy;				//rIy
		x[1] = Structure_P.q1;				//q1
		x[2] = Structure_P.q2;				//q2
		x[3] = Structure_P.q3;				//q3
		x[4] = Structure_P.q4;				//q4
		x[5] = Structure_P.q5;				//q5
		x[6] = Structure_P.q6;				//q6
		x[7] = Structure_P.rIydot;			//rIy
		x[8] = Structure_P.q1dot;			//q1
		x[9] = Structure_P.q2dot;			//q2
		x[10] = Structure_P.q3dot;			//q3
		x[11] = Structure_P.q4dot;			//q4
		x[12] = Structure_P.q5dot;			//q5
		x[13] = Structure_P.q6dot;			//q6

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

		Robot_State_Init[2-1] = x[0];				//rIy
		Robot_State_Init[4-1] = x[1];				//q1
		Robot_State_Init[5-1] = x[2];				//q2
		Robot_State_Init[6-1] = x[3];				//q3
		Robot_State_Init[7-1] = x[4];				//q4
		Robot_State_Init[8-1] = x[5];				//q5
		Robot_State_Init[9-1] = x[6];				//q6
		Robot_State_Init[2+13-1] = x[7];			//rIy
		Robot_State_Init[4+13-1] = x[8];			//q1
		Robot_State_Init[5+13-1] = x[9];			//q2
		Robot_State_Init[6+13-1] = x[10];			//q3
		Robot_State_Init[7+13-1] = x[11];			//q4
		Robot_State_Init[8+13-1] = x[12];			//q5
		Robot_State_Init[9+13-1] = x[13];			//q6

		delete []iAfun;  delete []jAvar;  delete []A;
		delete []iGfun;  delete []jGvar;

		delete []x;      delete []xlow;   delete []xupp;
		delete []xmul;   delete []xstate;

		delete []F;      delete []Flow;   delete []Fupp;
		delete []Fmul;   delete []Fstate;

		delete []xnames; delete []Fnames;

		return Robot_State_Init;
	}
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

	std::vector<double> Robot_State_Init_vec = Structure_P.Robot_State_Init;

	Robot_StateNDot StateNDot_Init_i(Robot_State_Init_vec);

	StateNDot_Init_i.rIy = x[0];			//rIy
	StateNDot_Init_i.q1 = x[1];				//q1
	StateNDot_Init_i.q2 = x[2];				//q2
	StateNDot_Init_i.q3 = x[3];				//q3
	StateNDot_Init_i.q4 = x[4];				//q4
	StateNDot_Init_i.q5 = x[5];				//q5
	StateNDot_Init_i.q6 = x[6];				//q6
	StateNDot_Init_i.rIydot = x[7];			//rIydot
	StateNDot_Init_i.q1dot = x[8];			//q1dot
	StateNDot_Init_i.q2dot = x[9];			//q2dot
	StateNDot_Init_i.q3dot = x[10];			//q3dot
	StateNDot_Init_i.q4dot = x[11];			//q4dot
	StateNDot_Init_i.q5dot = x[12];			//q5dot
	StateNDot_Init_i.q6dot = x[13];			//q6dot

	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_Init_i, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_Init_i, "rI");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");

	std::vector<double> rCOM = Ang_Pos_fn(StateNDot_Init_i, "rCOM");

	std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");
	std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
	std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");
	std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
	std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
	std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");

	std::vector<double> sigma_i = Structure_P.sigma_i;
	double sigma0_1 = 1.0*sigma_i[0];
	double sigma0_2 = 1.0*sigma_i[1];

	// for (int i = 0; i < 17; i++)
	// {
	// 	F[i] = 0.0;
	// }
	F[0] = -rT[1];
	// std::cout<<sigma0_1<<endl;
	// if(sigma0_1==1)
	// {
	// double val1 = -sigma0_1*(mini - rC[1]);
	// double val2 = -sigma0_1*(mini - rD[1]);

		F[1] = sigma0_1*rA[1];
		F[2] = sigma0_1*rB[1];
		F[3] = sigma0_1*vA[0];
		F[4] = sigma0_1*vA[1];
		F[5] = sigma0_1*vB[0];
		F[6] = sigma0_1*vB[1];
		F[7] = -sigma0_1*(sigma0_2==0)*(mini - rC[1]);
		F[8] = -sigma0_1*(sigma0_2==0)*(mini - rD[1]);
	// }
	// if(sigma0_2==1)
	// {
		F[9] =  sigma0_2*rC[1];
		F[10] = sigma0_2*rD[1];
		F[11] = sigma0_2*vC[0];
		F[12] = sigma0_2*vC[1];
		F[13] = sigma0_2*vD[0];
		F[14] = sigma0_2*vD[1];
		F[15] = -sigma0_2*(sigma0_1==0)*(mini - rA[1]);
		F[16] = -sigma0_2*(sigma0_1==0)*(mini - rB[1]);
	// }

	return 0;
}

Unified_Structure_P::Unified_Structure_P()
{
	// A default constructor
	// vector<double> Robot_State_Init;

	ifstream Initial_Robot_State_File;              // This is to read the initial angle and angular velocities
	Initial_Robot_State_File.open("init_robot_state.txt");
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
		printf("Unable to initialize anything yo!\n");
	}

	rIx = Robot_State_Init[0];
	rIy = Robot_State_Init[1];
	theta = Robot_State_Init[2];
	q1 = Robot_State_Init[3];
	q2 = Robot_State_Init[4];
	q3 = Robot_State_Init[5];
	q4 = Robot_State_Init[6];
	q5 = Robot_State_Init[7];
	q6 = Robot_State_Init[8];
	q7 = Robot_State_Init[9];
	q8 = Robot_State_Init[10];
	q9 = Robot_State_Init[11];
	q10 = Robot_State_Init[12];
	rIxdot = Robot_State_Init[13];
	rIydot = Robot_State_Init[14];
	thetadot = Robot_State_Init[15];
	q1dot = Robot_State_Init[16];
	q2dot = Robot_State_Init[17];
	q3dot = Robot_State_Init[18];
	q4dot = Robot_State_Init[19];
	q5dot = Robot_State_Init[20];
	q6dot = Robot_State_Init[21];
	q7dot = Robot_State_Init[22];
	q8dot = Robot_State_Init[23];
	q9dot = Robot_State_Init[24];
	q10dot = Robot_State_Init[25];

	for (int i = 0; i < 10; i++) {
		/* code */
		Opt_Ctrl_LowBd.push_back(0);
		Opt_Ctrl_UppBd.push_back(0);
	}
	for (int i = 0; i < 26; i++) {
		Opt_Conf_LowBd.push_back(0);
		Opt_Conf_UppBd.push_back(0);
	}

	Opt_Conf_LowBd[0] = rIxlow; 				Opt_Conf_UppBd[0] = rIxupp;
	Opt_Conf_LowBd[1] = rIylow; 				Opt_Conf_UppBd[1] = rIyupp;
	Opt_Conf_LowBd[2] = thetalow; 				Opt_Conf_UppBd[2] = thetaupp;
	Opt_Conf_LowBd[3] = q1low; 					Opt_Conf_UppBd[3] = q1upp;
	Opt_Conf_LowBd[4] = q2low; 					Opt_Conf_UppBd[4] = q2upp;
	Opt_Conf_LowBd[5] = q3low; 					Opt_Conf_UppBd[5] = q3upp;
	Opt_Conf_LowBd[6] = q4low; 					Opt_Conf_UppBd[6] = q4upp;
	Opt_Conf_LowBd[7] = q5low; 					Opt_Conf_UppBd[7] = q5upp;
	Opt_Conf_LowBd[8] = q6low; 					Opt_Conf_UppBd[8] = q6upp;
	Opt_Conf_LowBd[9] = q7low; 					Opt_Conf_UppBd[9] = q7upp;
	Opt_Conf_LowBd[10] = q8low; 				Opt_Conf_UppBd[10] = q8upp;
	Opt_Conf_LowBd[11] = q9low; 				Opt_Conf_UppBd[11] = q9upp;
	Opt_Conf_LowBd[12] = q10low; 				Opt_Conf_UppBd[12] = q10upp;
	Opt_Conf_LowBd[0+13] = rIxdotlow; 			Opt_Conf_UppBd[0+13] = rIxdotupp;
	Opt_Conf_LowBd[1+13] = rIydotlow; 			Opt_Conf_UppBd[1+13] = rIydotupp;
	Opt_Conf_LowBd[2+13] = thetadotlow; 		Opt_Conf_UppBd[2+13] = thetadotupp;
	Opt_Conf_LowBd[3+13] = q1dotlow; 			Opt_Conf_UppBd[3+13] = q1dotupp;
	Opt_Conf_LowBd[4+13] = q2dotlow; 			Opt_Conf_UppBd[4+13] = q2dotupp;
	Opt_Conf_LowBd[5+13] = q3dotlow; 			Opt_Conf_UppBd[5+13] = q3dotupp;
	Opt_Conf_LowBd[6+13] = q4dotlow; 			Opt_Conf_UppBd[6+13] = q4dotupp;
	Opt_Conf_LowBd[7+13] = q5dotlow; 			Opt_Conf_UppBd[7+13] = q5dotupp;
	Opt_Conf_LowBd[8+13] = q6dotlow; 			Opt_Conf_UppBd[8+13] = q6dotupp;
	Opt_Conf_LowBd[9+13] = q7dotlow; 			Opt_Conf_UppBd[9+13] = q7dotupp;
	Opt_Conf_LowBd[10+13] = q8dotlow; 			Opt_Conf_UppBd[10+13] = q8dotupp;
	Opt_Conf_LowBd[11+13] = q9dotlow; 			Opt_Conf_UppBd[11+13] = q9dotupp;
	Opt_Conf_LowBd[12+13] = q10dotlow; 			Opt_Conf_UppBd[12+13] = q10dotupp;

	Opt_Ctrl_LowBd[0] = -tau1_max;				Opt_Ctrl_UppBd[0] = tau1_max;
	Opt_Ctrl_LowBd[1] = -tau2_max;				Opt_Ctrl_UppBd[1] = tau2_max;
	Opt_Ctrl_LowBd[2] = -tau3_max;				Opt_Ctrl_UppBd[2] = tau3_max;
	Opt_Ctrl_LowBd[3] = -tau4_max;				Opt_Ctrl_UppBd[3] = tau4_max;
	Opt_Ctrl_LowBd[4] = -tau5_max;				Opt_Ctrl_UppBd[4] = tau5_max;
	Opt_Ctrl_LowBd[5] = -tau6_max;				Opt_Ctrl_UppBd[5] = tau6_max;
	Opt_Ctrl_LowBd[6] = -tau7_max;				Opt_Ctrl_UppBd[6] = tau7_max;
	Opt_Ctrl_LowBd[7] = -tau8_max;				Opt_Ctrl_UppBd[7] = tau8_max;
	Opt_Ctrl_LowBd[8] = -tau9_max;				Opt_Ctrl_UppBd[8] = tau9_max;
	Opt_Ctrl_LowBd[9] = -tau10_max;				Opt_Ctrl_UppBd[9] = tau10_max;
}

Robot_StateNDot::Robot_StateNDot()
{
	// A default constructor
	rIx = 0;
	rIy = 0.7230;
	theta = -0.0900;
	q1 = 0.3768;
	q2 = 0.0045;
	q3 = -0.2913;
	q4 = -1.0015;
	q5 = 0.1500;
	q6 = 0.2698;
	q7 = -0.6600;
	q8 = -0.6251;
	q9 = 0.6900;
	q10 = -0.2951;
	rIxdot = 0.2000;
	rIydot = -0.0605;
	thetadot = -0.2100;
	q1dot = -0.1239;
	q2dot = 1.3108;
	q3dot = -0.9768;
	q4dot = -1.4999;
	q5dot = 2.0000;
	q6dot = -1.2999;
	q7dot = 1.0000;
	q8dot = -2.0000;
	q9dot = -1.5708;
	q10dot = -1.5000;
}

Robot_StateNDot::Robot_StateNDot(std::vector<double> &Robot_AngleNRate)
{
    // An evaluated constructor
    rIx = Robot_AngleNRate[0];
    rIy = Robot_AngleNRate[1];
    theta = Robot_AngleNRate[2];
    q1 = Robot_AngleNRate[3];
    q2 = Robot_AngleNRate[4];
    q3 = Robot_AngleNRate[5];
    q4 = Robot_AngleNRate[6];
    q5 = Robot_AngleNRate[7];
    q6 = Robot_AngleNRate[8];
    q7 = Robot_AngleNRate[9];
    q8 = Robot_AngleNRate[10];
    q9 = Robot_AngleNRate[11];
    q10 = Robot_AngleNRate[12];

    rIxdot = Robot_AngleNRate[13];
    rIydot = Robot_AngleNRate[14];
    thetadot = Robot_AngleNRate[15];
    q1dot = Robot_AngleNRate[16];
    q2dot = Robot_AngleNRate[17];
    q3dot = Robot_AngleNRate[18];
    q4dot = Robot_AngleNRate[19];
    q5dot = Robot_AngleNRate[20];
    q6dot = Robot_AngleNRate[21];
    q7dot = Robot_AngleNRate[22];
    q8dot = Robot_AngleNRate[23];
    q9dot = Robot_AngleNRate[24];
    q10dot = Robot_AngleNRate[25];
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
// Here are the position vector functions
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

	dlib::matrix<double>  T(13,13);
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
	dlib::matrix<double> B_q(13,10);
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

	dlib::matrix<double> T(13,1);
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

	dlib::matrix<double> T(12,13);
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

	dlib::matrix<double> T(12,1);
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

double Obs_Dist_fn(std::vector<double> &r_pos, const char* s)
{
	// To make the problem easy to solve, here we only consider the flat/vertical environmental obstacles
	double x, y, Dist_temp;
	x = r_pos[0];
	y = r_pos[1];

	std::vector<double> Dist_Vec;

	for (int i = 0; i < Envi_Map.nr(); i++)
	{
		double x1 = Envi_Map(i,0);
		double y1 = Envi_Map(i,1);
		double x2 = Envi_Map(i,2);
		double y2 = Envi_Map(i,3);

		Dist_temp = 100.0;

		if(strcmp(s,"x")==0)
		{
			if(x1 == x2)
			{
				Dist_temp = x1 - x;
			}
		}
		else
		{
			if(strcmp(s,"y")==0)
			{
				if(y1 == y2)
				{
					Dist_temp = y - y1;
				}
			}
			else
			{
				Dist_temp = x1 - x;
				Dist_temp = min(Dist_temp, y - y1);
			}
		}
		Dist_Vec.push_back(Dist_temp);
	}
	return *std::min_element(Dist_Vec.begin(), Dist_Vec.end());
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
	if(strcmp(s,"vE")==0)
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
int Node_Expansion_fn(Tree_Node &Cur_Node, Unified_Structure_P Structure_P)
{
	// This function is used to expansion a certain node to its nearest children nodes
	/* This function is the main function used to expansion the given node to
	% its adjacent nodes without any connectivity test

	%%      Inputs;
	%           sigma_i:        the contact status at the time i
	%           x0:             the robot state at time i
	%           P:              the pre-load structure

	%%      Output:
	%           sigma_children: the updated queue after a node expansion

	%%      The main algorithm
	%
	%       Hand contact: 0-> Try the kinematical maximum step length to test the collision
	%                           if collision detected?
	%                               Expanded with adding one hand contact
	%                     1-> Since one hand is in contact,
	%                               Expanded with adding the other hand contact
	%                                             removing the current hand contact
	%                     2-> Now two hands are in contact,
	%                               Expanded with retracting either hand contact
	%       Foot contact: 0-> Add either foot contact point next
	%                     1-> Since one foot is in contact, adding one or
	%                     removing one
	%                     2-> Remove either foot contact*/

	std::vector<double> sigma_i = Cur_Node.sigma_i;

	double foot_AB_contas = sigma_i[1-1];
	double foot_CD_contas = sigma_i[2-1];
	double hand_E_contas = sigma_i[3-1];
	double hand_F_contas = sigma_i[4-1];
	// cout<<"foot_AB_contas value  is  "<<foot_AB_contas<<endl;
	// cout<<"foot_CD_contas value  is  "<<foot_CD_contas<<endl;
	// cout<<"hand_E_contas value  is  "<<hand_E_contas<<endl;
	// cout<<"hand_F_contas value  is  "<<hand_F_contas<<endl;


	Robot_StateNDot StateNDot_Init_i = Cur_Node.StateNDot_Str;

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

	std::vector<double> vA = Ang_Pos_fn(StateNDot_Init_i, "vA");
	std::vector<double> vB = Ang_Pos_fn(StateNDot_Init_i, "vB");
	std::vector<double> vC = Ang_Pos_fn(StateNDot_Init_i, "vC");
	std::vector<double> vD = Ang_Pos_fn(StateNDot_Init_i, "vD");
	std::vector<double> vE = Ang_Pos_fn(StateNDot_Init_i, "vE");
	std::vector<double> vF = Ang_Pos_fn(StateNDot_Init_i, "vF");

	// Hand contact expansion
	if ((hand_E_contas == 0)&&(hand_F_contas == 0))
	{
		// 1. No contact case
		std::vector<double> End_Hori_Array, End_Vert_Array;

		End_Hori_Array.push_back(rA[0]);		End_Vert_Array.push_back(rA[1]);
		End_Hori_Array.push_back(rB[0]);		End_Vert_Array.push_back(rB[1]);
		End_Hori_Array.push_back(rC[0]);		End_Vert_Array.push_back(rC[1]);
		End_Hori_Array.push_back(rD[0]);		End_Vert_Array.push_back(rD[1]);
		End_Hori_Array.push_back(rE[0]);		End_Vert_Array.push_back(rE[1]);
		End_Hori_Array.push_back(rF[0]);		End_Vert_Array.push_back(rF[1]);

		double End_Hori_Max = *std::max_element(End_Hori_Array.begin(), End_Hori_Array.end());
		double End_Vert_Min = *std::min_element(End_Vert_Array.begin(), End_Vert_Array.end());
		std::vector<double> r_pos(2); r_pos[0] = End_Hori_Max; r_pos[1] = End_Vert_Min;

		double Outreach = Obs_Dist_fn(r_pos, "x")- Step_Length_Max;
		if (Outreach<0)
		{
			// In this case, there could be a hand collision
			Tree_Node Node_Child1, Node_Child2;
			std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, 2, 1);
			std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, 3, 1);
			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
		}
	}
	else
	{
		//2. Two hand contact case
		if ((hand_E_contas == 1)&&(hand_F_contas == 1))
		{
			Tree_Node Node_Child1, Node_Child2;
			std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, 2, 0);
			std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, 3, 0);
			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
		}
		else
		{
			// 3. One hand contact case
			Tree_Node Node_Child1, Node_Child2, Node_Child3, Node_Child4;
			std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, 2, 1);
			std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, 3, 1);
			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
			std::vector<double> Node_Child3_sigma = sigma_modi(sigma_i, 2, 0);
			std::vector<double> Node_Child4_sigma = sigma_modi(sigma_i, 3, 0);
			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child3_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child4_sigma);

		}
	}

	// Foot contact expansion
	if ((foot_AB_contas == 0)&&(foot_CD_contas == 0))
	{
		//1. No contact case
		Tree_Node Node_Child1, Node_Child2;
		std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, 1, 1);
		std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, 2, 1);
		Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
		Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
	}
	else
	{
		// 2. Two foot contact case
		if ((foot_AB_contas == 1)&&(foot_CD_contas == 1))
		{
			Tree_Node Node_Child1, Node_Child2;
			std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, 1, 0);
			std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, 2, 0);
			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
		}
		else
		{
			int Active_Ind = 1, Inactive_Ind = 0;
			//3. One foot contact case
			if(sigma_i[0]>sigma_i[1])
			{
				Active_Ind = 0;
				Inactive_Ind = 1;
			}
			Tree_Node Node_Child1, Node_Child2;
			std::vector<double> Node_Child1_sigma = sigma_modi(sigma_i, Active_Ind, 0);
			std::vector<double> Node_Child2_sigma = sigma_modi(sigma_i, Inactive_Ind, 1);

			Add_Child2Par(Node_Child1,Cur_Node,Nodes_Tot_Number,Node_Child1_sigma);
			Add_Child2Par(Node_Child2,Cur_Node,Nodes_Tot_Number,Node_Child2_sigma);
		}
	}
	return 0;
}

std::vector<double> sigma_modi(std::vector<double> sigma_ref, int contas_ind, int AddOrRet)
{
	std::vector<double> sigma_child = sigma_ref;
	sigma_child[contas_ind] = 1.0*AddOrRet;
	return sigma_child;
}

int Nodes_Connectivity_Opt(Tree_Node &Node_i, Tree_Node &Node_i_child, Unified_Structure_P &Structure_P)
{
	/*% This function test the connectivity between a certain node and its child  node

	% The current idea is to make sure of the multiple shooting method since in
	% this case the dynamics constraints are satisfied automatically

	% The main idea to reach the sigma_child at the end step while minimizing the kinetic energy*/
	int Ctrl_No = 20;						// This is the default no of controls within a mode
	double Tme_Seed = 2;   					// This is the default time period for each segment



}

int Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child, Unified_Structure_P &Structure_P, int Ctrl_No, double Tme_Seed)
{
	int flag = 0;			// The default flag is set to be no-solution found
	// This function is used to initialize the seed guess for the optimization
	// The initialized Ctrl_Traj and StateNDot_Traj will be saved into  Node_i_child

	// Three stages need to be conducted to complete this whole initialization process
	vector<double> F_Constraint_vec;
	vector<double> Constraint_Status_vec;

	Structure_P.Node_i = Node_i;
	Structure_P.Node_i_child = Node_i_child;

	Seed_Conf_Constraint(Node_i, Node_i.StateNDot_Str, Node_i_child.sigma_i, F_Constraint_vec, Constraint_Status_vec);


	snoptProblem Seed_Conf_Constraint_Pr;                     // This is the name of the Optimization problem

	// Allocate and initialize
	// 1. Optimization for a feasible configuration that satisfies the desired mode

	integer n = 26;
	integer neF = Constraint_Status_vec.size();     // 1 objective function
	integer lenA  =  n * neF;                              // This is the number of nonzero elements in the linear part A    F(x) = f(x)+Ax

	integer lenru = 14;             					  // This is used to pass the initial state into the usrfun
	doublereal *ru = new doublereal[lenru];
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
	// First set the lower and upper bounds for state
	xlow[0] = rIxlow;				xupp[0] = rIxlow;
	xlow[1] = rIylow;				xupp[1] = rIylow;
	xlow[2] = thetalow;				xupp[2] = thetalow;
	xlow[3] = q1low;				xupp[3] = q1upp;
	xlow[4] = q2low;				xupp[4] = q2upp;
	xlow[5] = q3low;				xupp[5] = q3upp;
	xlow[6] = q4low;				xupp[6] = q4upp;
	xlow[7] = q5low;				xupp[7] = q5upp;
	xlow[8] = q6low;				xupp[8] = q6upp;
	xlow[9] = q7low;				xupp[9] = q7upp;
	xlow[10] = q8low;				xupp[10] = q8upp;
	xlow[11] = q9low;				xupp[11] = q9upp;
	xlow[12] = q10low;				xupp[12] = q10upp;
	xlow[13] = rIxdotlow;			xupp[13] = rIxdotupp;
	xlow[14] = rIydotlow;			xupp[14] = rIydotupp;
	xlow[15] = thetadotlow;			xupp[15] = thetadotupp;
	xlow[16] = q1dotlow;			xupp[16] = q1dotupp;
	xlow[17] = q2dotlow;			xupp[17] = q2dotupp;
	xlow[18] = q3dotlow;			xupp[18] = q3dotupp;
	xlow[19] = q4dotlow;			xupp[19] = q4dotupp;
	xlow[20] = q5dotlow;			xupp[20] = q5dotupp;
	xlow[21] = q6dotlow;			xupp[21] = q6dotupp;
	xlow[22] = q7dotlow;			xupp[22] = q7dotupp;
	xlow[23] = q8dotlow;			xupp[23] = q8dotupp;
	xlow[24] = q9dotlow;			xupp[24] = q9dotupp;
	xlow[25] = q10dotlow;			xupp[25] = q10dotupp;

	// Second set the lower and upper bounds for the objective function
	for (int i = 0; i < Constraint_Status_vec.size(); i++)
	{
		int Constraint_Status_vec_i = Constraint_Status_vec[i];
		if(Constraint_Status_vec_i==0)
		{
			Flow[i] = 0.0;			Fupp[i] = 0.0;
		}
		else
		{
			if(Constraint_Status_vec_i==1)
			{
				Flow[i] = -Inf;			Fupp[i] = 0.0;
			}
			else
			{
				Flow[i] = 0.0;			Fupp[i] = Inf;

			}
		}
	}

	// Load the data for ToyProb ...
	Seed_Conf_Constraint_Pr.setPrintFile  ( "Seed_Conf_Constraint_Pr.out" );
	Seed_Conf_Constraint_Pr.setProblemSize( n, neF );
	Seed_Conf_Constraint_Pr.setObjective  ( ObjRow, ObjAdd );
	Seed_Conf_Constraint_Pr.setA          ( lenA, iAfun, jAvar, A );
	Seed_Conf_Constraint_Pr.setG          ( lenG, iGfun, jGvar );
	Seed_Conf_Constraint_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Seed_Conf_Constraint_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Seed_Conf_Constraint_Pr.setXNames     ( xnames, nxnames );
	Seed_Conf_Constraint_Pr.setFNames     ( Fnames, nFnames );
	Seed_Conf_Constraint_Pr.setProbName   ( "Seed_Conf_Constraint_Pr" );
	Seed_Conf_Constraint_Pr.setUserFun    ( Seed_Conf_Constraint_Pr_fn_);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Seed_Conf_Constraint_Pr.computeJac    ();
	Seed_Conf_Constraint_Pr.setIntParameter( "Derivative option", 0 );
	integer Cold = 0, Basis = 1, Warm = 2;
	Seed_Conf_Constraint_Pr.solve( Cold );

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return 0;
}

void Seed_Conf_Constraint(Tree_Node &Node_i, Robot_StateNDot &StateNDot_i_child, vector<double> &sigma_i_child, vector<double> &F, vector<double> &Constraint_Status)
{
	// This function is used to generate the constraint value needed for the evaluation of the seed configuration initialization
	vector<double> sigma_i = Node_i.sigma_i;

	F.push_back(Kinetic_Energy_fn(StateNDot_i_child));		Constraint_Status.push_back(2);

	std::vector<double> rA = Ang_Pos_fn(StateNDot_i_child, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_i_child, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_i_child, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_i_child, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_i_child, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_i_child, "rF");
	std::vector<double> rG = Ang_Pos_fn(StateNDot_i_child, "rG");
	std::vector<double> rH = Ang_Pos_fn(StateNDot_i_child, "rH");
	std::vector<double> rI = Ang_Pos_fn(StateNDot_i_child, "rI");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_i_child, "rT");

	std::vector<double> rCOM = Ang_Pos_fn(StateNDot_i_child, "rCOM");

	std::vector<double> vA = Ang_Vel_fn(StateNDot_i_child, "vA");
	std::vector<double> vB = Ang_Vel_fn(StateNDot_i_child, "vB");
	std::vector<double> vC = Ang_Vel_fn(StateNDot_i_child, "vC");
	std::vector<double> vD = Ang_Vel_fn(StateNDot_i_child, "vD");
	std::vector<double> vE = Ang_Vel_fn(StateNDot_i_child, "vE");
	std::vector<double> vF = Ang_Vel_fn(StateNDot_i_child, "vF");

	double sigma_i_AB = sigma_i[1-1];
	double sigma_i_CD = sigma_i[2-1];
	double sigma_i_E = sigma_i[3-1];
	double sigma_i_F = sigma_i[4-1];

	double sigma_i_child_AB = sigma_i_child[1-1];
	double sigma_i_child_CD = sigma_i_child[2-1];
	double sigma_i_child_E = sigma_i_child[3-1];
	double sigma_i_child_F = sigma_i_child[4-1];

	vector<double> sigma_offset(4);
	for (int i = 0; i < 4; i++)
	{
		sigma_offset[i] = sigma_i_child[i] - sigma_i[i];
	}
	/* 1. Relative distance  constraints:
	a. all distancehave to be at least on the surface
	b. the desired mode has to be satisfied.
	*/
	int sigma_i_max = *max_element(sigma_i.begin(),sigma_i.end());
	int sigma_i_min = *min_element(sigma_i.begin(),sigma_i.end());

	int sigma_i_child_max = *max_element(sigma_i_child.begin(),sigma_i_child.end());
	int sigma_i_child_min = *min_element(sigma_i_child.begin(),sigma_i_child.end());
	int sigma_offset_max = *max_element(sigma_offset.begin(),sigma_offset.end());
	int sigma_offset_sum = sigma_offset[0] +  sigma_offset[1] +  sigma_offset[2] +  sigma_offset[3];

	if ((sigma_i_child_max == sigma_i_child_min)&&(sigma_i_child_max==0))
	{
		std::vector<double> Node_i_child_StateNDot_vec = StateNDot2StateVec(StateNDot_i_child);
		std::vector<double> Node_i_StateNDot_vec = StateNDot2StateVec(Node_i.StateNDot_Str);
		F.push_back(sqrt(Dot_Product(Node_i_child_StateNDot_vec, Node_i_child_StateNDot_vec))-sqrt(Dot_Product(Node_i_StateNDot_vec, Node_i_StateNDot_vec)) - 10*mini);
		Constraint_Status.push_back(1);
		F.push_back(sigma_offset[1-1] * vA[2-1]);		Constraint_Status.push_back(1);
		F.push_back(sigma_offset[1-1] * vB[2-1]);		Constraint_Status.push_back(1);
		F.push_back(sigma_offset[2-1] * vC[2-1]);		Constraint_Status.push_back(1);
		F.push_back(sigma_offset[2-1] * vD[2-1]);		Constraint_Status.push_back(1);
		F.push_back(-sigma_offset[3-1] * vE[1-1]);		Constraint_Status.push_back(1);
		F.push_back(-sigma_offset[4-1] * vF[1-1]);		Constraint_Status.push_back(1);
	}
	// The absolute relative distance between the robot end effector and the environment obstacles
	F.push_back(-Obs_Dist_fn(rA,"z") + mini * (!sigma_i_child_AB));		Constraint_Status.push_back(1);
	F.push_back(-Obs_Dist_fn(rB,"z") + mini * (!sigma_i_child_AB));		Constraint_Status.push_back(1);
	F.push_back(-Obs_Dist_fn(rC,"z") + mini * (!sigma_i_child_CD));		Constraint_Status.push_back(1);
	F.push_back(-Obs_Dist_fn(rD,"z") + mini * (!sigma_i_child_CD));		Constraint_Status.push_back(1);
	F.push_back(-Obs_Dist_fn(rE,"z") + mini * (!sigma_i_child_E));		Constraint_Status.push_back(1);
	F.push_back(-Obs_Dist_fn(rF,"z") + mini * (!sigma_i_child_F));		Constraint_Status.push_back(1);

	// The complementarity condition
	F.push_back(sigma_i_child_AB * rA[1]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_AB * vA[0]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_AB * vA[1]);		Constraint_Status.push_back(0);

	F.push_back(sigma_i_child_AB * rB[1]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_AB * vB[0]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_AB * vB[1]);		Constraint_Status.push_back(0);

	F.push_back(sigma_i_child_CD * rC[1]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_CD * vC[0]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_CD * vC[1]);		Constraint_Status.push_back(0);

	F.push_back(sigma_i_child_CD * rD[1]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_CD * vD[0]);		Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_CD * vD[1]);		Constraint_Status.push_back(0);

	F.push_back(sigma_i_child_E * Obs_Dist_fn(rE, "x"));	Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_E * vE[0]);	Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_E * vE[1]);	Constraint_Status.push_back(0);

	F.push_back(sigma_i_child_F * Obs_Dist_fn(rF, "x"));	Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_F * vF[0]);	Constraint_Status.push_back(0);
	F.push_back(sigma_i_child_F * vF[1]);	Constraint_Status.push_back(0);

	//  2. Contact Constraint Maintenance: the previous contacts have to be satisfied
	if(sigma_offset_max  ==0)
	{
		F.push_back((sigma_offset[1-1]==0) * sigma_i_AB * (rA[0] - Node_i.rA_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[1-1]==0) * sigma_i_AB * (rA[1] - Node_i.rA_ref[1]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[1-1]==0) * sigma_i_AB * (rB[0] - Node_i.rB_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[1-1]==0) * sigma_i_AB * (rB[1] - Node_i.rB_ref[1]));	Constraint_Status.push_back(0);

		F.push_back((sigma_offset[2-1]==0) * sigma_i_CD * (rC[0] - Node_i.rC_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[2-1]==0) * sigma_i_CD * (rC[1] - Node_i.rC_ref[1]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[2-1]==0) * sigma_i_CD * (rD[0] - Node_i.rD_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[2-1]==0) * sigma_i_CD * (rD[1] - Node_i.rD_ref[1]));	Constraint_Status.push_back(0);

		F.push_back((sigma_offset[3-1]==0) * sigma_i_E * (rE[0] - Node_i.rE_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[3-1]==0) * sigma_i_E * (rE[1] - Node_i.rE_ref[1]));	Constraint_Status.push_back(0);

		F.push_back((sigma_offset[4-1]==0) * sigma_i_F * (rF[0] - Node_i.rF_ref[0]));	Constraint_Status.push_back(0);
		F.push_back((sigma_offset[4-1]==0) * sigma_i_F * (rF[1] - Node_i.rF_ref[1]));	Constraint_Status.push_back(0);
	}
	else
	{
		F.push_back(sigma_i_AB * (rA[0] - Node_i.rA_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_AB * (rA[1] - Node_i.rA_ref[1]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_AB * (rB[0] - Node_i.rB_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_AB * (rB[1] - Node_i.rB_ref[1]));	Constraint_Status.push_back(0);

		F.push_back(sigma_i_CD * (rC[0] - Node_i.rC_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_CD * (rC[1] - Node_i.rC_ref[1]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_CD * (rD[0] - Node_i.rD_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_CD * (rD[1] - Node_i.rD_ref[1]));	Constraint_Status.push_back(0);

		F.push_back(sigma_i_E * (rE[0] - Node_i.rE_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_E * (rE[1] - Node_i.rE_ref[1]));	Constraint_Status.push_back(0);

		F.push_back(sigma_i_F * (rF[0] - Node_i.rF_ref[0]));	Constraint_Status.push_back(0);
		F.push_back(sigma_i_F * (rF[1] - Node_i.rF_ref[1]));	Constraint_Status.push_back(0);
	}
	// 3. Heuristic Stability Constraints: rI and rCOM have to be lied within the support polygon
	vector<double> r_Foot_Pos(4);	r_Foot_Pos.push_back(rA[0]);	r_Foot_Pos.push_back(rB[0]);	r_Foot_Pos.push_back(rC[0]);	r_Foot_Pos.push_back(rD[0]);
	double Temp_val1, Temp_val2;
	int Stance_Leg_Ind, Swing_Leg_Ind, Stance_Leg_Index;double Swing_Leg_Dis, Stance_Leg_Dis;
	if((sigma_i_max == sigma_i_min)&&(sigma_i_max==0))
	{
		F.push_back(rCOM[0] - rT[0]);			Constraint_Status.push_back(0);
	}
	else
	{
		if((sigma_i[0]==1)||(sigma_i[1]==1))
		{
			//        % At least foot contact is involved
			if(sigma_offset_max>0)
			{
				if((std::abs(sigma_offset[0])==1)||(std::abs(sigma_offset[1]==1)))
				{
					if(Node_i.vI_ref[0]>0)
					{
						Swing_Leg_Ind = 0;
						if(abs(sigma_offset[1]==1))
						{
							Swing_Leg_Ind = 1;
						}
						Stance_Leg_Ind = (Swing_Leg_Ind==0) + 1;
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						Swing_Leg_Dis = min(Temp_val1, Temp_val2);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						Stance_Leg_Dis = max(Temp_val1, Temp_val2);
						F.push_back(Stance_Leg_Dis - Swing_Leg_Dis);			Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						F.push_back(min(Temp_val1, Temp_val2) - rI[0]);			Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						F.push_back(rI[0] - max(Temp_val1, Temp_val2));			Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						F.push_back(min(Temp_val1, Temp_val2) - rCOM[0]);		Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						F.push_back(rCOM[0] - max(Temp_val1, Temp_val2));		Constraint_Status.push_back(1);
						F.push_back(min(rI[0], rCOM[0]) -rT[0]);		Constraint_Status.push_back(1);
						F.push_back(rT[0] - max(rI[0], rCOM[0]));		Constraint_Status.push_back(1);
					}
					else
					{
						Swing_Leg_Ind = 0;
						if(abs(sigma_offset[1]==1))
						{
							Swing_Leg_Ind = 1;
						}
						Stance_Leg_Ind = (Swing_Leg_Ind==0) + 1;
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						Swing_Leg_Dis = max(Temp_val1, Temp_val2);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						Stance_Leg_Dis = min(Temp_val1, Temp_val2);
						F.push_back(Stance_Leg_Dis - Swing_Leg_Dis);			Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						F.push_back(-max(Temp_val1, Temp_val2) + rI[0]);		Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						F.push_back(-rI[0] + min(Temp_val1, Temp_val2));		Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Ind+1];
						F.push_back(-max(Temp_val1, Temp_val2) + rCOM[0]);		Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Swing_Leg_Ind];
						Temp_val2 = r_Foot_Pos[2*Swing_Leg_Ind+1];
						F.push_back(-rCOM[0] + min(Temp_val1, Temp_val2));
						Constraint_Status.push_back(1);
						F.push_back(min(rI[0], rCOM[0]) -rT[0]);		Constraint_Status.push_back(1);
						F.push_back(rT[0] - max(rI[0], rCOM[0]));		Constraint_Status.push_back(1);
					}
				}
				else
				{
					F.push_back(rCOM[0]- rT[0]);			Constraint_Status.push_back(0);

				}
			}
			else
			{
				if(sigma_offset_sum==0)
				// In this case, we are doing the same mode optimization and we
				// would like the center of mass and rI to be within the support polygon
				{
					if(sigma_i[0] + sigma_i[1]==2)
					{
						Temp_val1 = *min_element(r_Foot_Pos.begin(),r_Foot_Pos.end());
						F.push_back(Temp_val1 - rI[0]);
						Constraint_Status.push_back(1);
						Temp_val1 =*max_element(r_Foot_Pos.begin(),r_Foot_Pos.end());
						F.push_back(rI[0] - Temp_val1);
						Constraint_Status.push_back(1);
						Temp_val1 = *min_element(r_Foot_Pos.begin(),r_Foot_Pos.end());
						F.push_back(Temp_val1 - rCOM[0]);
						Constraint_Status.push_back(1);
						Temp_val1 = *max_element(r_Foot_Pos.begin(),r_Foot_Pos.end());
						F.push_back(rCOM[0] - Temp_val1);
						Constraint_Status.push_back(1);
					}
					else
					{
						Stance_Leg_Index = 0;
						if(sigma_i[1]==1)
						{
							Stance_Leg_Index = 1;
						}
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Index];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Index+1];
						double Temp_val = min(Temp_val1, Temp_val2) - rI[0];
						F.push_back(Temp_val);
						Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Index];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Index+1];
						Temp_val = rI[0] - max(Temp_val1, Temp_val2);
						F.push_back(Temp_val);
						Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Index];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Index+1];
						Temp_val = min(Temp_val1, Temp_val2) - rCOM[0];
						F.push_back(Temp_val);
						Constraint_Status.push_back(1);
						Temp_val1 = r_Foot_Pos[2*Stance_Leg_Index];
						Temp_val2 = r_Foot_Pos[2*Stance_Leg_Index+1];
						Temp_val = rCOM[0] - max(Temp_val1, Temp_val2);
						F.push_back(Temp_val);
						Constraint_Status.push_back(1);
					}
				}
				else
				{
					F.push_back(sigma_offset[0] * vA[1]);			Constraint_Status.push_back(1);
					F.push_back(sigma_offset[0] * vB[1]);			Constraint_Status.push_back(1);
					F.push_back(sigma_offset[1] * vC[1]);			Constraint_Status.push_back(1);
					F.push_back(sigma_offset[1] * vD[1]);			Constraint_Status.push_back(1);
					F.push_back(sigma_offset[3] * vE[0]);			Constraint_Status.push_back(1);
					F.push_back(sigma_offset[4] * vF[0]);			Constraint_Status.push_back(1);
				}
			}
		}
	}
}

double Dot_Product(vector<double> &x1, vector<double> &x2)
{
	double result = 0.0;
	// Thisfunction is used to compute the doc product between two vectors
	for (int i = 0; i < x1.size(); i++)
	{
		result = result + x1[i] * x2[i];
	}
	return result;
}

void Reference_Dist_Vel_Update(Tree_Node &Node_i)
{
	// Thyis funciton is used to update the value of the reference in the given node
	Robot_StateNDot Node_i_StateNDot = Node_i.StateNDot_Str;

	Node_i.rA_ref = Ang_Pos_fn(Node_i_StateNDot, "rA");
	Node_i.rB_ref = Ang_Pos_fn(Node_i_StateNDot, "rB");
	Node_i.rC_ref = Ang_Pos_fn(Node_i_StateNDot, "rC");
	Node_i.rD_ref = Ang_Pos_fn(Node_i_StateNDot, "rD");
	Node_i.rE_ref = Ang_Pos_fn(Node_i_StateNDot, "rE");
	Node_i.rF_ref = Ang_Pos_fn(Node_i_StateNDot, "rF");

	Node_i.vA_ref = Ang_Pos_fn(Node_i_StateNDot, "vA");
	Node_i.vB_ref = Ang_Pos_fn(Node_i_StateNDot, "vB");
	Node_i.vC_ref = Ang_Pos_fn(Node_i_StateNDot, "vC");
	Node_i.vD_ref = Ang_Pos_fn(Node_i_StateNDot, "vD");
	Node_i.vE_ref = Ang_Pos_fn(Node_i_StateNDot, "vE");
	Node_i.vF_ref = Ang_Pos_fn(Node_i_StateNDot, "vF");
	Node_i.vI_ref = Ang_Pos_fn(Node_i_StateNDot, "vI");
	return;
}
int Seed_Conf_Constraint_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
	     integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru )
		 {

			 vector<double> x_i_child;

			 for (int i = 0; i < 26; i++)
			 {
				 x_i_child.push_back(x[i]);
			 }


			 vector<double> F_Constraint_vec;
			 vector<double> Constraint_Status_vec;

			 Robot_StateNDot x_i_child_StateNDot = StateVec2StateNDot(x_i_child);

			 Seed_Conf_Constraint(Structure_P.Node_i, x_i_child_StateNDot, Structure_P.Node_i_child.sigma_i, F_Constraint_vec, Constraint_Status_vec);

			 for (int i = 0; i < F_Constraint_vec.size(); i++)
			 {
				 F[i] = F_Constraint_vec[i];
			 }

			 return 0;
		 }
