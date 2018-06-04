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

// There are three types of variables in this optimization problem: robot state, control torques and contact forces

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

// 2. Control torques
double tau1_max = 100;             		double tau2_max = 100;
double tau3_max = 100;					double tau4_max = 100;
double tau5_max = 100;             		double tau6_max = 100;
double tau7_max = 60;              		double tau8_max = 50;
double tau9_max = 60;             		double tau10_max = 50;

// Save them into the bounds
dlib::matrix<double,26,1> xlow_vec;
dlib::matrix<double,26,1> xupp_vec;
dlib::matrix<double,10,1> ctrl_low_vec;
dlib::matrix<double,10,1> ctrl_upp_vec;
dlib::matrix<double,12,1> contact_force_low_vec;
dlib::matrix<double,12,1> contact_force_upp_vec;

int Total_Nodes_Number = 0; // The initializatin of the total nodes currently in the tree
dlib::matrix<double> Envi_Map;

/**
 * Some global values are defined
 * Description
 */

double mini = 0.05;			int Opt_Var_Per_Frame = 48;		int Constraints_Per_Frame = 59;
int Ctrl_No = 20;			double Tme_Seed = 0.5;			double mu = 0.5;

std::vector<Tree_Node_Ptr> All_Nodes;				// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		    // The kinetic energy of each nodes

std::vector<double> StateNDot_Seed;
std::vector<double> Ctrl_Seed;

dlib::matrix<double> Envi_Map_Defi()
{
	// 	This function is used to define the environment map for the simulation
	// Whenever this function gets called, it will return the array with the
	// environment obstacle information
	//
	// This is the default flat ground
	// This map is defined in a polyline manner with the first value denoting
	// the line length and the second value denoting the relative angle

	Envi_Map = dlib::ones_matrix<double>(2,2);
	Envi_Map(0,0) = 5.0;		Envi_Map(0,1) = 0.0;
	Envi_Map(1,0) = 3.0;		Envi_Map(1,1) = PI/2.0;
	return Envi_Map;
}

void Add_Node(Tree_Node &Current_Node)
{
	// This function will add the current node to the All_Nodes vector
	All_Nodes.push_back(&Current_Node);
	Frontier_Nodes.push_back(&Current_Node);
	Frontier_Nodes_Cost.push_back(Current_Node.Kinetic_Energy);
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

int Add_Child2Par(Tree_Node &Child_Node, Tree_Node &Par_Node, int &Total_Nodes_Number, std::vector<double> &sigma_i)
{

    Child_Node.Par_Node = &Par_Node;

    Child_Node.Node_Index_Number = Total_Nodes_Number + 1;

	Child_Node.sigma_i = sigma_i;

	Children_Nodes.push_back(&Child_Node);

    Total_Nodes_Number = Total_Nodes_Number + 1;

	Par_Node.Children_Nodes.push_back(&Child_Node);

    return 0;
}

std::vector<double> Default_Init(const std::vector<double> &sigma_i, int Flag)
{
	// This function is used to initialize the whole optimization process
	// First, is to substitute the map info into the Envi_Map matrix
	// Second, is to give the proper bounds to the variables to be optimized
	// Thrid, is to generate a kinematically feasible initial robot state

	// First job finished!
	Envi_Map = Envi_Map_Defi();

	// Second job finished!
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

	contact_force_low_vec = -dlib::ones_matrix<double>(12,1) * Inf;
	contact_force_upp_vec =  dlib::ones_matrix<double>(12,1) * Inf;


	// Third job started!
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

	if(Flag == 1)// This means that the given initial condition works for the constraint
	{
		return Robot_State_Init;
	}
	else
	{
		// If the default configuration would like to be viewed

		// Robot_StateNDot Robot_StateNDot_init(Robot_State_Init);
		// std::string input_name = "init_given";
		// Robot_Plot_fn(Robot_StateNDot_init,input_name);
		snoptProblem Default_Init_Pr;                     // This is the name of the Optimization problem
		// Allocate and initialize
		integer n = 26;
		integer neF = 25;     							  // 1 objective function
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
		// First set the lower and upper bounds for state
		xlow[0] = rIxlow;		xupp[0] = rIxupp;
		xlow[1] = rIylow;		xupp[1] = rIyupp;
		xlow[2] = thetalow;		xupp[2] = thetaupp;
		xlow[3] = q1low;		xupp[3] = q1upp;
		xlow[4] = q2low;		xupp[4] = q2upp;
		xlow[5] = q3low;		xupp[5] = q3upp;
		xlow[6] = q4low;		xupp[6] = q4upp;
		xlow[7] = q5low;		xupp[7] = q5upp;
		xlow[8] = q6low;		xupp[8] = q6upp;
		xlow[9] = q7low;		xupp[9] = q7upp;
		xlow[10] = q8low;		xupp[10] = q8upp;
		xlow[11] = q9low;		xupp[11] = q9upp;
		xlow[12] = q10low;		xupp[12] = q10upp;

		// Second set the lower and upper bounds for state
		xlow[0+13] = rIxdotlow;			xupp[0+13] = rIxdotupp;
		xlow[1+13] = rIydotlow;			xupp[1+13] = rIydotupp;
		xlow[2+13] = thetadotlow;		xupp[2+13] = thetadotupp;
		xlow[3+13] = q1dotlow;			xupp[3+13] = q1dotupp;
		xlow[4+13] = q2dotlow;			xupp[4+13] = q2dotupp;
		xlow[5+13] = q3dotlow;			xupp[5+13] = q3dotupp;
		xlow[6+13] = q4dotlow;			xupp[6+13] = q4dotupp;
		xlow[7+13] = q5dotlow;			xupp[7+13] = q5dotupp;
		xlow[8+13] = q6dotlow;			xupp[8+13] = q6dotupp;
		xlow[9+13] = q7dotlow;			xupp[9+13] = q7dotupp;
		xlow[10+13] = q8dotlow;			xupp[10+13] = q8dotupp;
		xlow[11+13] = q9dotlow;			xupp[11+13] = q9dotupp;
		xlow[12+13] = q10dotlow;		xupp[12+13] = q10dotupp;

		for (int i = 0; i < 26; i++)
		{
			xstate[i] = 0.0;
		}

		for(int i = 0; i<neF; i++)
		{
			// The lower bound is the same
			Flow[i] = 0.0;
		}

		// Second set the lower and upper bounds for the objective function
		Fupp[0] = Inf;		// This is the difference between the optimized configuration and the given configuration

		// The constraint bounds should be carefully defined:
		// They are actually 3 * 6 constraints.

		double sigma_i_j, Fupp_val, Flow_index, Fupp_index;

		for (int j = 0; j < 4; j++)
		{
			if(j < 2)
			{
				// The first two indicate the foot contact point
				sigma_i_j = sigma_i[j];
				if(sigma_i_j>0)	// In this case, the contact constraint should be active
				{
					Fupp_val = 0.0;
					Flow_index = 1 + j * 8;
					Fupp_index = Flow_index + 8;

					for (int i = Flow_index; i < Fupp_index; i++)
					{
						Fupp[i] = Fupp_val;
					}
				}
				else			// In this case, the contact constraints should be inactive
				{
					Fupp_val = Inf;
					Flow_index = 1 + j * 8;
					Fupp_index = Flow_index + 8;
					for (int i = Flow_index; i < Fupp_index; i++)
					{
						Fupp[i] = Fupp_val;
					}
				}
			}
			else
			{
				// The first two indicate the foot contact point
				sigma_i_j = sigma_i[j];
				if(sigma_i_j>0)	// In this case, the contact constraint should be active
				{
					Fupp_val = 0.0;
					Flow_index = 17 + (j-2) * 4;
					Fupp_index = Flow_index + 4;

					for (int i = Flow_index; i < Fupp_index; i++)
					{
						Fupp[i] = Fupp_val;
					}
				}
				else			// In this case, the contact constraints should be inactive
				{
					Fupp_val = Inf;
					Flow_index = 17 + (j-2) * 4;
					Fupp_index = Flow_index + 4;
					for (int i = Flow_index; i < Fupp_index; i++)
					{
						Fupp[i] = Fupp_val;
					}
				}
			}
		}
		// Initial guess
		for (int i = 0; i < 26; i++)
		{
			x[i] = Robot_State_Init[i];
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
		for (int i = 0; i < 26; i++)
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
	// Initial guess of the robot configurations
	std::vector<double> Robot_State_Init = Structure_P.Robot_State_Init;
	std::vector<double> Robot_State_Opt;
	for (int i = 0; i < 26; i++)
	{
		Robot_State_Opt.push_back(x[i]);
	}

	Robot_StateNDot StateNDot_Init_i(Robot_State_Opt);

	// These are the positions of the robot end effectors
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");

	// These are the positions of the intermediate joints and the head joint
	std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");

	std::vector<double> rCOM = Ang_Pos_fn(StateNDot_Init_i, "rCOM");

	std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");
	std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
	std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");
	std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
	std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
	std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");

	std::vector<double> sigma_i = Structure_P.sigma_i;

	double xG_offset = (Robot_State_Init[0]- Robot_State_Opt[0]) * (Robot_State_Init[0]- Robot_State_Opt[0]);
	double yG_offset = (Robot_State_Init[1]- Robot_State_Opt[1]) * (Robot_State_Init[1]- Robot_State_Opt[1]);

	// The optimized configuration should obey the initial xG and yG position
	F[0] = xG_offset + yG_offset;

	// 2 foot contacts
	F[1] = sigma_i[0]*rA[1];
	F[2] = sigma_i[0]*vA[0];
	F[3] = sigma_i[0]*vA[1];
	F[4] = (sigma_i[0]==0) * (rA[1] - mini);

	F[5] = sigma_i[0]*rB[1];
	F[6] = sigma_i[0]*vB[0];
	F[7] = sigma_i[0]*vB[1];
	F[8] = (sigma_i[0]==0) * (rB[1] - mini);

	F[9] = sigma_i[1]*rC[1];
	F[10] = sigma_i[1]*vC[0];
	F[11] = sigma_i[1]*vC[1];
	F[12] = (sigma_i[1]==0) * (rC[1] - mini);

	F[13] = sigma_i[1]*rD[1];
	F[14] = sigma_i[1]*vD[0];
	F[15] = sigma_i[1]*vD[1];
	F[16] = (sigma_i[1]==0) * (rD[1] - mini);

	// Two hand contacts
	int Obs_Choice_Ind; char Obs_sym = 'z';
	F[17] = sigma_i[2]*Obs_Dist_Fn(rE, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[18] = sigma_i[2]*Obs_Dist_Fn(rE, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[19] = sigma_i[2]*Obs_Dist_Fn(rE, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[20] = (sigma_i[2]==0) * (Obs_Dist_Fn(rE, Envi_Map, Obs_Choice_Ind, Obs_sym) - mini);

	F[21] = sigma_i[3]*Obs_Dist_Fn(rF, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[22] = sigma_i[3]*Obs_Dist_Fn(rF, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[23] = sigma_i[3]*Obs_Dist_Fn(rF, Envi_Map, Obs_Choice_Ind, Obs_sym);
	F[24] = (sigma_i[3]==0) * (Obs_Dist_Fn(rF, Envi_Map, Obs_Choice_Ind, Obs_sym) - mini);
	return 0;
}

double Obs_Dist_Fn(std::vector<double> &r_Pos, dlib::matrix<double> &Envi_Map, int &Obs_Choice_Ind, char char_sym)
{
	// 	This function is used to calculate the relative distance between the robot end effector and the nearby environment
	// To make it easier for the research purpose, here we only consider two environmental obstacles:  flat ground and a vertical wall
	double Flat_Grnd_Vert = 0.0;
	double Vert_Wall_Hori = Envi_Map(0,0);

	double r_Pos_x_offset = Vert_Wall_Hori - r_Pos[0];
	double r_Pos_y_offset = r_Pos[1] - Flat_Grnd_Vert;

	double Obs_Dist;

	// The difference between char_sym and Obs_Dist_Ind is that char_sym is specified by the user while Obs_Choice_Ind is determiend by the algorithm

	if(char_sym == 'x')
	{
		Obs_Dist = r_Pos_x_offset;
	}
	else
	{
		if( char_sym == 'y')
		{
			Obs_Dist = r_Pos_y_offset;
		}
		else
		{
			if ((r_Pos_x_offset * r_Pos_x_offset)>(r_Pos_y_offset * r_Pos_y_offset))
			{
				Obs_Choice_Ind = 0;
				Obs_Dist = r_Pos_y_offset;
			}
			else
			{
				Obs_Choice_Ind = 1;
				Obs_Dist = r_Pos_x_offset;
			}
		}
	}
	return Obs_Dist;
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
	B_q = dlib::zeros_matrix<double>(13,1);
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

}

std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
	int flag = 0;			// The default flag is set to be no-solution found
	std::vector<double> Robot_State_i, Robot_State_i_child, sigma_tran, sigma_goal;
	Robot_State_i = StateNDot2StateVec(Node_i.Node_StateNDot);

	Structure_P.Node_i = Node_i;
	Structure_P.Node_i_child = Node_i_child;

	// This function is used to initialize the seed guess for the optimization
	// The initialized Ctrl_Traj and StateNDot_Traj will be saved into  Node_i_child

	// Due to the elimination of the flying-in-air inertia shaping approach, here we will not consider it for now

	// Three stages need to be conducted to complete this whole initialization process

/*
	// Allocate and initialize
	// 1. Optimization for a feasible configuration that satisfies the desired mode
*/
	Robot_State_i_child = Seed_Guess_Gene_Robotstate(Node_i, Node_i_child);
	// Here the desired configuration has already been generated, then the state, control and control force need to be optimized
	dlib::matrix<double> Robot_State_Tot(26,Ctrl_No), Robot_State_Interpol_Array;
	dlib::matrix<double> Ctrl_Tot(10,Ctrl_No-1), Contact_Force_Tot(12,Ctrl_No-1);

	// Here the robot state has already been interpolated into Robot_State_Tot into the row fashion

	for (int i = 0; i < 26; i++)
	{
		Robot_State_Interpol_Array = dlib::linspace(Robot_State_i[i], Robot_State_i_child[i], Ctrl_No);

		for (int j = 0; j < Ctrl_No; j++)
		{
			Robot_State_Tot(i,j) = Robot_State_Interpol_Array(j);
		}
	}

	// Two more varaibles to be optimized: control and contact force and they will be initialied randomly

	Ctrl_Tot = dlib::randm(10,Ctrl_No-1);	Contact_Force_Tot = dlib::randm(12,Ctrl_No-1);

	// Finally is to concatenate them into a single column vector
	std::vector<double> Opt_Seed;
	Opt_Seed.push_back(Tme_Seed); // The first element is the time duration on each segment

	// Then it is the robot state
	for (int i = 0; i < Ctrl_No-1; i++)
	{
		for (int j = 0; j < 26; j++)
		{
			Opt_Seed.push_back(Robot_State_Tot(j,i+1));
		}
		for (int j = 0; j< 10; j++)
		{
			Opt_Seed.push_back(Ctrl_Tot(j,i));
		}
		for (int j = 0; j < 12; j++)
		{
			Opt_Seed.push_back(Contact_Force_Tot(j,i));
			/* code */
		}
	}

	return Opt_Seed;
}

std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child)
{
	// This function is used to generate a configuration to initialize the optimization

	std::vector<double> Robot_State_2BOpt = StateNDot2StateVec(Node_i.Node_StateNDot);

	snoptProblem Seed_Conf_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = 26;
	integer neF = 65;     							  // 1 objective function
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

	// Set the upper and lower bounds.
	// First set the lower and upper bounds for state
	xlow[0] = rIxlow;		xupp[0] = rIxupp;
	xlow[1] = rIylow;		xupp[1] = rIyupp;
	xlow[2] = thetalow;		xupp[2] = thetaupp;
	xlow[3] = q1low;		xupp[3] = q1upp;
	xlow[4] = q2low;		xupp[4] = q2upp;
	xlow[5] = q3low;		xupp[5] = q3upp;
	xlow[6] = q4low;		xupp[6] = q4upp;
	xlow[7] = q5low;		xupp[7] = q5upp;
	xlow[8] = q6low;		xupp[8] = q6upp;
	xlow[9] = q7low;		xupp[9] = q7upp;
	xlow[10] = q8low;		xupp[10] = q8upp;
	xlow[11] = q9low;		xupp[11] = q9upp;
	xlow[12] = q10low;		xupp[12] = q10upp;

	// Second set the lower and upper bounds for state
	xlow[0+13] = rIxdotlow;			xupp[0+13] = rIxdotupp;
	xlow[1+13] = rIydotlow;			xupp[1+13] = rIydotupp;
	xlow[2+13] = thetadotlow;		xupp[2+13] = thetadotupp;
	xlow[3+13] = q1dotlow;			xupp[3+13] = q1dotupp;
	xlow[4+13] = q2dotlow;			xupp[4+13] = q2dotupp;
	xlow[5+13] = q3dotlow;			xupp[5+13] = q3dotupp;
	xlow[6+13] = q4dotlow;			xupp[6+13] = q4dotupp;
	xlow[7+13] = q5dotlow;			xupp[7+13] = q5dotupp;
	xlow[8+13] = q6dotlow;			xupp[8+13] = q6dotupp;
	xlow[9+13] = q7dotlow;			xupp[9+13] = q7dotupp;
	xlow[10+13] = q8dotlow;			xupp[10+13] = q8dotupp;
	xlow[11+13] = q9dotlow;			xupp[11+13] = q9dotupp;
	xlow[12+13] = q10dotlow;		xupp[12+13] = q10dotupp;

	for (int i = 0; i < 26; i++)
	{
		xstate[i] = 0.0;
	}

	for(int i = 0; i<neF; i++)
	{
		// The lower bound is the same
		Flow[i] = 0.0;
		Fupp[i] = 0.0;
	}

	// Second set the lower and upper bounds for the objective function
	Fupp[0] = Inf;		// This is the difference between the optimized configuration and the given configuration

	// The constraint bounds should be carefully defined:
	// They are actually 3 * 6 constraints.

	for (int i = 25; i < 33; i++)
	{
		Fupp[i] = Inf;
	}

	// Initial guess
	for (int i = 0; i < 26; i++)
	{
		x[i] = Robot_State_2BOpt[i];
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
		Robot_State_2BOpt[i] = x[i];
	}

	Robot_StateNDot Init_Opt_vec(Robot_State_2BOpt);
	// Robot_Plot_fn(Init_Opt_vec);

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return Robot_State_2BOpt;
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
	Robot_StateNDot Node_i_StateNDot = Node_i.Node_StateNDot;
	End_Effector_PosNVel(Node_i_StateNDot, Node_i.End_Effector_Pos, Node_i.End_Effector_Vel);
	return;
}

int Node_Self_Opt(Tree_Node &Node_i)
{
	// This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
	int Optimization_Result = 0;	// 1 means it is possible to find a path that satisfies the constraints while minimizing the objective function value

	// Here there are two different methods to initialize the optimization seeds: 1. pure randomized 2. best guess interpolation

	std::vector<double> Opt_Seed = Seed_Guess_Gene(Node_i, Node_i);

	std::vector<double> sigma_i_change = Vec_Minus(Node_i.sigma_i, Node_i_child.sigma_i);
	double sigma_i_change_val = sigma_i_change[0] + sigma_i_change[1] + sigma_i_change[2] + sigma_i_change[3];
	if(sigma_i_change_val ==1)
	{
		// % In this case, it is making contact
		sigma_tran = Node_i.sigma_i;
		sigma_goal = Node_i_child.sigma_i;
	}
	else
	{
		// % In this case, it is retracting contact or maintaining the same	contact condition
		sigma_tran = Node_i_child.sigma_i;
		sigma_goal = Node_i_child.sigma_i;
	}
	Structure_P.sigma_tran = sigma_tran;
	Structure_P.sigma_goal = sigma_goal;

	// Now it is the pure



	return yes;

}

std::vector<double> Real_Optimization(std::vector<double> &Opt_Seed, std::vector<double> &sigma_tran, std::vector<double> &sigma_goal, int &Opt_Flag)
{
	// This function undertakes all the optimization computation burden

	snoptProblem Real_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Opt_Seed.size();
	integer neF = 65;     							  // 1 objective function
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

	// Set the upper and lower bounds.
	// First set the lower and upper bounds for state
	xlow[0] = rIxlow;		xupp[0] = rIxupp;
	xlow[1] = rIylow;		xupp[1] = rIyupp;
	xlow[2] = thetalow;		xupp[2] = thetaupp;
	xlow[3] = q1low;		xupp[3] = q1upp;
	xlow[4] = q2low;		xupp[4] = q2upp;
	xlow[5] = q3low;		xupp[5] = q3upp;
	xlow[6] = q4low;		xupp[6] = q4upp;
	xlow[7] = q5low;		xupp[7] = q5upp;
	xlow[8] = q6low;		xupp[8] = q6upp;
	xlow[9] = q7low;		xupp[9] = q7upp;
	xlow[10] = q8low;		xupp[10] = q8upp;
	xlow[11] = q9low;		xupp[11] = q9upp;
	xlow[12] = q10low;		xupp[12] = q10upp;

	// Second set the lower and upper bounds for state
	xlow[0+13] = rIxdotlow;			xupp[0+13] = rIxdotupp;
	xlow[1+13] = rIydotlow;			xupp[1+13] = rIydotupp;
	xlow[2+13] = thetadotlow;		xupp[2+13] = thetadotupp;
	xlow[3+13] = q1dotlow;			xupp[3+13] = q1dotupp;
	xlow[4+13] = q2dotlow;			xupp[4+13] = q2dotupp;
	xlow[5+13] = q3dotlow;			xupp[5+13] = q3dotupp;
	xlow[6+13] = q4dotlow;			xupp[6+13] = q4dotupp;
	xlow[7+13] = q5dotlow;			xupp[7+13] = q5dotupp;
	xlow[8+13] = q6dotlow;			xupp[8+13] = q6dotupp;
	xlow[9+13] = q7dotlow;			xupp[9+13] = q7dotupp;
	xlow[10+13] = q8dotlow;			xupp[10+13] = q8dotupp;
	xlow[11+13] = q9dotlow;			xupp[11+13] = q9dotupp;
	xlow[12+13] = q10dotlow;		xupp[12+13] = q10dotupp;

	for (int i = 0; i < 26; i++)
	{
		xstate[i] = 0.0;
	}

	for(int i = 0; i<neF; i++)
	{
		// The lower bound is the same
		Flow[i] = 0.0;
		Fupp[i] = 0.0;
	}

	// Second set the lower and upper bounds for the objective function
	Fupp[0] = Inf;		// This is the difference between the optimized configuration and the given configuration

	// The constraint bounds should be carefully defined:
	// They are actually 3 * 6 constraints.

	for (int i = 25; i < 33; i++)
	{
		Fupp[i] = Inf;
	}

	// Initial guess
	for (int i = 0; i < 26; i++)
	{
		x[i] = Robot_State_2BOpt[i];
	}

	// Load the data for ToyProb ...
	Real_Optimization_Pr.setPrintFile  ( "Real_Optimization_Pr.out" );
	Real_Optimization_Pr.setProblemSize( n, neF );
	Real_Optimization_Pr.setObjective  ( ObjRow, ObjAdd );
	Real_Optimization_Pr.setA          ( lenA, iAfun, jAvar, A );
	Real_Optimization_Pr.setG          ( lenG, iGfun, jGvar );
	Real_Optimization_Pr.setX          ( x, xlow, xupp, xmul, xstate );
	Real_Optimization_Pr.setF          ( F, Flow, Fupp, Fmul, Fstate );
	Real_Optimization_Pr.setXNames     ( xnames, nxnames );
	Real_Optimization_Pr.setFNames     ( Fnames, nFnames );
	Real_Optimization_Pr.setProbName   ( "Real_Optimization_Pr_" );
	Real_Optimization_Pr.setUserFun    ( Seed_Conf_Optimization_Pr_fn_);
	// snopta will compute the Jacobian by finite-differences.
	// The user has the option of calling  snJac  to define the
	// coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
	Real_Optimization_Pr.computeJac    ();
	Real_Optimization_Pr.setIntParameter( "Derivative option", 0 );
	integer Cold = 0, Basis = 1, Warm = 2;
	Real_Optimization_Pr.solve( Cold );

	for (int i = 0; i < 26; i++)
	{
		Robot_State_2BOpt[i] = x[i];
	}

	Robot_StateNDot Init_Opt_vec(Robot_State_2BOpt);
	// Robot_Plot_fn(Init_Opt_vec);

	delete []iAfun;  delete []jAvar;  delete []A;
	delete []iGfun;  delete []jGvar;

	delete []x;      delete []xlow;   delete []xupp;
	delete []xmul;   delete []xstate;

	delete []F;      delete []Flow;   delete []Fupp;
	delete []Fmul;   delete []Fstate;

	delete []xnames; delete []Fnames;

	return Robot_State_2BOpt;
}

std::vector<double> Opt_Constraint(double Tme_Seed, dlib:matrix<double> &Robot_State_Tot, dlib:matrix<double> &Ctrl_Tot, dlib:matrix<double> &Contact_Force_Tot, int Flag_Choice, std::vector<double> &Constraint_Type)
{
	// This function works to provide the constraint used for the optimization also help initalization the optimization problem

	// Here the Flag_Choice has three possibilities: 0------------------------> Node_Self_Opt
	//												 1------------------------> Nodes_Connectivity_Opt
	//												 2------------------------> Pure Constraint formulation

	double Obj_val = 0, Obj_val_i;		// Initialize this value to be 0 and  will be updated later
	double Not_Sigma_t_1, Not_Sigma_t_2, Not_Sigma_t_3, Not_Sigma_t_4, Kinetic_Energy_t, Kinetic_Energy_tp1;
	std::vector<double> Contact_Ind(8), sigma_tran, sigma_goal, sigma_temp;

	sigma_tran = Structure_P.sigma_tran;		sigma_goal = Structure_P.sigma_goal;

	Robot_StateNDot Robot_StateNDot_i;			double delta_t = Tme_Seed/(Ctrl_No * 1.0);

	Robot_State_init = Structure_P.Node_i.Node_StateNDot;

	std::vector<double> Opt_Constraint_vals, StateNDot_vec_i;
	Opt_Constraint_vals.push_back(Obj_val);
	Constraint_Type.push_back(1);			// The fist Constraint is an inequallity constraint

	dlib::matrix<double> D_q_i, B_q_i, C_q_qdot_i, Jac_Full_i, x0, x0p1, x0_state, x0_statedot, x0p1_state, x0p1_statedot, temp_result, qddot_i, Ctrl_i, Contact_Force_i, End_Effector_Obs_Dist;
	dlib::matrix<double> End_Eq_Pos_Matrix, Inq_Pos_Matrix, Eqn_Maint_Matrix, temp_matrix, Normal_Forces, Tang_Forces, End_Eq_Force_Matrix;
	dlib::matrix<double,16,1> End_Effector_Pos, End_Effector_Vel;
	for (int i = 0; i < Ctrl_No-1; i++)
	{
		// Take out the robot from the state vector
		if(i ==0){
			x0 = Robot_StateNDot2DlibMat(Robot_State_init);
			x0p1 = dlib::colm(Robot_State_Tot,0);}
		else{
			x0 = dlib::colm(Robot_State_Tot,i-1);
			x0p1 = dlib::colm(Robot_State_Tot,i);}
		StateNDot_vec_i = DlibMatCol2StdVec(x0p1, 26);
		Robot_StateNDot_i = StateVec2StateNDot(StateNDot_vec_i);

		D_q_i = D_q_fn(Robot_StateNDot_i)
		B_q_i  = B_q_fn();
		C_q_qdot_i = C_q_qdot_fn(Robot_StateNDot_i);
		Jac_Full_i = Jac_Full_fn(Robot_StateNDot_i);

		//
		//				1. First equality constraint is the dynamics update: qp1 	  = q    + qp1dot  * delta * t
		//																	   qdotp1 = qdot + qddotp1 * delta * t
		StateNStatedot_Distill(x0, x0_state, x0_statedot);
		StateNStatedot_Distill(x0p1, x0p1_state, x0p1_statedot);

		temp_result = x0_state + delta_t *x0p1_statedot - x0p1_state;
		Constraints2ValsNType(temp_result, Opt_Constraint_vals,Constraint_Type,Flag_Choice, 0);
		//
		//				2. Second equality constraint is the dynamics constraint
		//
		qddot_i = (x0p1_statedot - x0_statedot)/delta_t;
		Ctrl_i = dlib::colm(Ctrl_Tot, i);
		Contact_Force_i = dlib::colm(Contact_Force_Tot, i);
		temp_result = D_q_i * qddot_i + C_q_qdot_i - dlib::trans(Jac_Full_i) * Contact_Force_i - B_q_i * Ctrl_i;
		Constraints2ValsNType(temp_result, Opt_Constraint_vals,Constraint_Type,Flag_Choice, 0);
		//
		//				3. Distance constraints: certain end effectors have to be active/inactive: Two implications: active has to be zero and inactive has to be nonzero
		//
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
		End_Effector_Obs_Dist = End_Effector_Obs_Dist_Fn(End_Effector_Pos, Contact_Ind);

		Inq_Pos_Matrix = dlib::zeros_matrix<double>(8,8);
		temp_matrix = dlib::ones_matrix<double>(8,1);

		if(i == Ctrl_No-2){
			sigma_temp = sigma_goal;}
		else{
			sigma_temp = sigma_tran;}

		Inq_Pos_Matrix(0,0) = !sigma_temp[0];        Inq_Pos_Matrix(1,1) = !sigma_temp[0];
		Inq_Pos_Matrix(2,2) = !sigma_temp[1];        Inq_Pos_Matrix(3,3) = !sigma_temp[1];
		Inq_Pos_Matrix(4,4) = !sigma_temp[2];        Inq_Pos_Matrix(5,5) = !sigma_temp[3];

		temp_result = Inq_Pos_Matrix * (End_Effector_Obs_Dist - temp_matrix * mini);
		Constraints2ValsNType(temp_result, Opt_Constraint_vals,Constraint_Type,Flag_Choice, 1);

		if(i == Ctrl_No-2){
			End_Eq_Pos_Matrix = dlib::zeros_matrix<double>(8,8);
			End_Eq_Pos_Matrix(0,0) = sigma_temp[0];        End_Eq_Pos_Matrix(1,1) = sigma_temp[0];
	        End_Eq_Pos_Matrix(2,2) = sigma_temp[1];        End_Eq_Pos_Matrix(3,3) = sigma_temp[1];
	        End_Eq_Pos_Matrix(4,4) = sigma_temp[2];        End_Eq_Pos_Matrix(5,5) = sigma_temp[3];
			temp_result = End_Eq_Pos_Matrix * End_Effector_Obs_Dist;
			Constraints2ValsNType(temp_result, Opt_Constraint_vals,Constraint_Type,Flag_Choice, 0);}
		//
		//				4. Complementarity constraints: Contact Force!
		//
		dlib::matrix<double> Normal_Force_Ones_Matrix = dlib::ones_matrix<double>(6,1);
		Contact_Force_Proj(Contact_Force_i, Contact_Ind, Normal_Forces, Tang_Forces);
		temp_result = Normal_Forces - Normal_Force_Ones_Matrix * mini;	// Normal forces have to be nonnegative
		Constraints2ValsNType(temp_result, Opt_Constraint_vals,Constraint_Type,Flag_Choice, 1);

	    // The Complementarity condition requires the product of contact force with relative distance to be equal to zero
		End_Eq_Force_Matrix = dlib::zeros_matrix<double>(12,12);		// This matrix helps select the force that have to be zero
		End_Eq_Force_Matrix(0,0) = !sigma_temp[0];				End_Eq_Force_Matrix(1,1) = !sigma_temp[0];
		End_Eq_Force_Matrix(2,2) = !sigma_temp[0];				End_Eq_Force_Matrix(3,3) = !sigma_temp[0];

		End_Eq_Force_Matrix(4,4) = !sigma_temp[1];				End_Eq_Force_Matrix(5,5) = !sigma_temp[1];
		End_Eq_Force_Matrix(6,6) = !sigma_temp[1];				End_Eq_Force_Matrix(7,7) = !sigma_temp[1];

		End_Eq_Force_Matrix(8,8) = !sigma_temp[2];				End_Eq_Force_Matrix(9,9) = !sigma_temp[2];
		End_Eq_Force_Matrix(10,10) = !sigma_temp[3];			End_Eq_Force_Matrix(11,11) = !sigma_temp[3];
		temp_result = End_Eq_Force_Matrix * Contact_Force_i;	// Normal forces have to be nonnegative
		Constraints2ValsNType(temp_result, Opt_Constraint_vals, Constraint_Type, Flag_Choice, 0);
		//
		//				5. Friction cone constraints:
		//
		temp_result = Friction_Cone_Constraint(Normal_Forces, Tang_Forces);
		Constraints2ValsNType(temp_result, Opt_Constraint_vals, Constraint_Type, Flag_Choice, 1);
		//
		//				6. Contact Constraint maintenance: the previous contacts have to be maintained
		//
		Eqn_Maint_Matrix = dlib::matrix<double>(16,16);
		Eqn_Maint_Matrix(0,0) = sigma_tran[0] * sigma_goal[0];
		Eqn_Maint_Matrix(1,1) = sigma_tran[0] * sigma_goal[0];
		Eqn_Maint_Matrix(2,2) = sigma_tran[0] * sigma_goal[0];
		Eqn_Maint_Matrix(3,3) = sigma_tran[0] * sigma_goal[0];
		Eqn_Maint_Matrix(4,4) = sigma_tran[1] * sigma_goal[1];
		Eqn_Maint_Matrix(5,5) = sigma_tran[1] * sigma_goal[1];
		Eqn_Maint_Matrix(6,6) = sigma_tran[1] * sigma_goal[1];
		Eqn_Maint_Matrix(7,7) = sigma_tran[1] * sigma_goal[1];
		Eqn_Maint_Matrix(8,8) = sigma_tran[2] * sigma_goal[2];
		Eqn_Maint_Matrix(9,9) = sigma_tran[2] * sigma_goal[2];
		Eqn_Maint_Matrix(10,10) = sigma_tran[3] * sigma_goal[3];
		Eqn_Maint_Matrix(11,11) = sigma_tran[3] * sigma_goal[3];
		dlib::matrix<double,16,1> ref_End_Effector_Pos = Structure_P.Node_i.End_Effector_Pos;
		temp_result = Eqn_Maint_Matrix * (ref_End_Effector_Pos - End_Effector_Pos);
		Constraints2ValsNType(temp_result, Opt_Constraint_vals, Constraint_Type, Flag_Choice, 0);
		//
		//				7. Kinetic energy rate constraint
		//
		if(i == 0){
			Kinetic_Energy_t = Structure_P.Node_i.Kinetic_Energy;
			Kinetic_Energy_tp1 = Kinetic_Energy_fn(StateVec2StateNDot(DlibMatCol2StdVec(dlib::colm(Robot_State_Tot,0))));}
		else{
			Kinetic_Energy_t = Kinetic_Energy_fn(StateVec2StateNDot(DlibMatCol2StdVec(dlib::colm(Robot_State_Tot,i))));
			Kinetic_Energy_tp1 = Kinetic_Energy_fn(StateVec2StateNDot(DlibMatCol2StdVec(dlib::colm(Robot_State_Tot,i+1))));}
		Obj_val_i = (Kinetic_Energy_tp1 - Kinetic_Energy_t)*(Kinetic_Energy_tp1 - Kinetic_Energy_t)/(delta_t * delta_t);
		Obj_val = Obj_val + Obj_val_i;
	}
	// The last two constraints are objective function and the final kinetic energy constraint
	Opt_Constraint_vals[0] = Obj_val;
	dlib::matrix<double> Kinetic_Energy_constraint(1);
	if(Flag_Choice == 0){
		// In this case, it is the Node_Self_Opt
		Kinetic_Energy_constraint = 0.1 - Kinetic_Energy_tp1;}
	else{
		// In this case, it is the connectivity constraint
		Kinetic_Energy_constraint = 0.1*Structure_P.Node_i.Kinetic_Energy - Kinetic_Energy_tp1;}
	temp_result = Kinetic_Energy_constraint;
	Constraints2ValsNType(temp_result, Opt_Constraint_vals, Constraint_Type, Flag_Choice, 1);
	return Opt_Constraint_vals;
}

dlib::matrix<double> Friction_Cone_Constraint(dlib::matrix<double> &Normal_Forces, dlib::matrix<double> &Tang_Forces)
{
	// This function is used to compute the holonomic constraint
	dlib::matrix<double> Friction_Cone_Constraint_val(4);
	double Tang_Force_1, Tang_Force_2, Tang_Force_3, Tang_Force_4;
	double Norm_Force_1, Norm_Force_2, Norm_Force_3, Norm_Force_4;
	Tang_Force_1 = Tang_Force(0) + Tang_Force(1);
	Norm_Force_1 = Normal_Force(0) + Normal_Force(1);
	Tang_Force_2 = Tang_Force(2) + Tang_Force(3);
	Norm_Force_2 = Normal_Force(2) + Normal_Force(3);
	Tang_Force_3 = Tang_Force(4);
	Norm_Force_3 = Normal_Force(4);
	Tang_Force_4 = Tang_Force(5);
	Norm_Force_4 = Normal_Force(5);
	Friction_Cone_Constraint_val(0) = mu * mu * Norm_Force_1 * Norm_Force_1 - Tang_Force_1 * Tang_Force_1;
	Friction_Cone_Constraint_val(1) = mu * mu * Norm_Force_2 * Norm_Force_2 - Tang_Force_2 * Tang_Force_2;
	Friction_Cone_Constraint_val(2) = mu * mu * Norm_Force_3 * Norm_Force_3 - Tang_Force_3 * Tang_Force_3;
	Friction_Cone_Constraint_val(3) = mu * mu * Norm_Force_4 * Norm_Force_4 - Tang_Force_4 * Tang_Force_4;
	return Friction_Cone_Constraint_val;}

void Contact_Force_Proj(dlib::matrix<double> &Contact_Force_i, std::vector<double> & Contact_Ind, dlib::matrix<double> &Normal_Forces, dlib::matrix<double> &Tang_Force){
	// This function is used to project the contact force into the normal force and the tangential force
	for (int i = 0; i < 6; i++) {
		if(Contact_Ind[i]==0){
			Tang_Force[i] = Contact_Force_i[2*i];
			Normal_Forces[i] = Contact_Force_i[2*i+1];}
		else{
			Tang_Force[i] = Contact_Force_i[2*i+1];
			Normal_Forces[i] = Contact_Force_i[2*i];}
	}
}

void Constraints2ValsNType(dlib::matrix<double> &temp_result, std::vector<double> &Opt_Constraint_vals, std::vector<double> &Constraint_Type, int Flag_Choice, int Constraint_Type_val)
{
	// This function is used to retrieve the computational result back to Vals and Constraint type
	for (int i = 0; i < temp_result.nr(); i++)
	{
		Opt_Constraint_vals.push_back(temp_result(i));
		if(Flag_Choice==2)
		{
			Constraint_Type.push_back(Constraint_Type_val);
		}
	}
}

void StateNStatedot_Distill(dlib::matrix<double> & MatCol, dilb::matrix<double> &MatState, dilb::matrix<double> &MatStatedot)
{
	// This function isused to distill the state and statedot from the vector
	for (int i = 0; i < 13; i++)
	{
		MatState(i) = MatCol(i);
		MatStatedot(i) = MatCol(i+13);
	}
}

dlib::matrix<double> Robot_StateNDot2DlibMat(Robot_StateNDot &Robot_StateNDot_i)
{
	dlib::matrix<double> Dlib_StateNDot_i(26,1);

	std::vector<double> StateVec = StateNDot2StateVec(Robot_StateNDot_i);

	for (int i = 0; i < 26; i++)
	{
		Dlib_StateNDot_i(i) = StateVec[i]
	}
	return Dlib_StateNDot_i;
}

std::vector<double> DlibMatCol2StdVec(dlib::matrix<double> &One_Col, int Col_Length)
{
	// This function is used to transform the dlib matrix column to a double type vector

	std::vector<double> StdVec(Col_Length);
	for (int i = 0; i < Col_Length; i++)
	{
		StdVec[i] = One_Col(i);
	}
	return StdVec;
}

int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru )
{
	// Initial guess of the robot configurations
	Tree_Node Node_i,  Node_i_child;		Node_i = Structure_P.Node_i;		Node_i_child = Structure_P.Node_i_child;
	std::vector<double> sigma_i, sigma_i_child;
	sigma_i = Node_i.sigma_i;				sigma_i_child = Node_i_child.sigma_i;
	dlib::matrix<double,16,1> ref_End_Effector_Pos = Node_i.End_Effector_Pos;
	dlib::matrix<double,16,1> ref_End_Effector_Vel = Node_i.End_Effector_Vel;

	std::vector<double> Robot_State_Opt;
	for (int i = 0; i < 26; i++)
	{
		Robot_State_Opt.push_back(x[i]);
	}

	Robot_StateNDot StateNDot_Init_i(Robot_State_Opt);
	dlib::matrix<double,16,1> End_Effector_Pos, End_Effector_Vel;
	End_Effector_PosNVel(StateNDot_Init_i, End_Effector_Pos, End_Effector_Vel);
	// These are the positions of the robot end effectors
	dlib::matrix<double> End_Effector_Obs_Dist = End_Effector_Obs_Dist_Fn(End_Effector_Pos);

	dlib::matrix<double> End_Effector_Pos_Constraint, End_Effector_Vel_Constraint, Inq_Pos_Constraint, Eqn_Maint_Constraint;

	// The optimized configuration should consider the kinetic energy
	F[0] = Kinetic_Energy_fn(StateNDot_Init_i);
	//	1.Relative distance constraints: a.all distance have to be at least on the surface b.the desired motion hasto be satisfied

	dlib::matrix<double,8,8> Eqn_Pos_Matrix, Inq_Pos_Matrix;
	Eqn_Pos_Matrix = dlib::zeros_matrix<double>(8,8);
	Inq_Pos_Matrix = dlib::zeros_matrix<double>(8,8);
	dlib::matrix<double,16,16> Eqn_Vel_Matrix, Eqn_Maint_Matrix;
	Eqn_Vel_Matrix = dlib::zeros_matrix<double>(16,16);
	Eqn_Maint_Matrix = dlib::zeros_matrix<double>(16,16);

	// printf("Eqn_Pos_Matrix Number of rows: %ld\n",Eqn_Pos_Matrix.nr());
	// printf("Eqn_Pos_Matrix Number of columns: %ld\n",Eqn_Pos_Matrix.nc());
	//
	// printf("Inq_Pos_Matrix Number of rows: %ld\n",Inq_Pos_Matrix.nr());
	// printf("Inq_Pos_Matrix Number of columns: %ld\n",Inq_Pos_Matrix.nc());
	//
	// printf("Eqn_Vel_Matrix Number of rows: %ld\n",Eqn_Vel_Matrix.nr());
	// printf("Eqn_Vel_Matrix Number of columns: %ld\n",Eqn_Vel_Matrix.nc());
	//
	// printf("Eqn_Maint_Matrix Number of rows: %ld\n",Eqn_Maint_Matrix.nr());
	// printf("Eqn_Maint_Matrix Number of columns: %ld\n",Eqn_Maint_Matrix.nc());

	Eqn_Pos_Matrix(0,0) = sigma_i_child[0];
	Eqn_Pos_Matrix(1,1) = sigma_i_child[0];

	Eqn_Pos_Matrix(2,2) = sigma_i_child[1];
	Eqn_Pos_Matrix(3,3) = sigma_i_child[1];

	Eqn_Pos_Matrix(4,4) = sigma_i_child[2];
	Eqn_Pos_Matrix(5,5) = sigma_i_child[3];

	////////////////////////////////////////////////////////////////////////////////
	Eqn_Vel_Matrix(0,0) = sigma_i_child[0];
	Eqn_Vel_Matrix(1,1) = sigma_i_child[0];
	Eqn_Vel_Matrix(2,2) = sigma_i_child[0];
	Eqn_Vel_Matrix(3,3) = sigma_i_child[0];

	Eqn_Vel_Matrix(4,4) = sigma_i_child[1];
	Eqn_Vel_Matrix(5,5) = sigma_i_child[1];
	Eqn_Vel_Matrix(6,6) = sigma_i_child[1];
	Eqn_Vel_Matrix(7,7) = sigma_i_child[1];

	Eqn_Vel_Matrix(8,8) = sigma_i_child[2];
	Eqn_Vel_Matrix(9,9) = sigma_i_child[2];
	Eqn_Vel_Matrix(10,10) = sigma_i_child[3];
	Eqn_Vel_Matrix(11,11) = sigma_i_child[3];

	End_Effector_Pos_Constraint = Eqn_Pos_Matrix * End_Effector_Obs_Dist;	 // 8 by 1
	End_Effector_Vel_Constraint = Eqn_Vel_Matrix * End_Effector_Vel; 		 // 16 by 1

	for (int i = 1; i < 9; i++)
	{
		F[i] = End_Effector_Pos_Constraint(i);
		// printf("End_Effector_Pos_Constraint(%d): %f\n", i, F[i]);
	}
	for (int i = 9; i < 25; i++)
	{
		F[i] = End_Effector_Vel_Constraint(i-9);
		// printf("End_Effector_Vel_Constraint(%d): %f\n", i - 9, F[i]);
	}

	Inq_Pos_Matrix(0,0) = !sigma_i_child[0];
	Inq_Pos_Matrix(1,1) = !sigma_i_child[0];
	Inq_Pos_Matrix(2,2) = !sigma_i_child[1];
	Inq_Pos_Matrix(3,3) = !sigma_i_child[1];
	Inq_Pos_Matrix(4,4) = !sigma_i_child[3];
	Inq_Pos_Matrix(5,5) = !sigma_i_child[4];

	dlib::matrix<double> End_Effector_Obs_vec = dlib::ones_matrix<double>(8,1);
	Inq_Pos_Constraint = Inq_Pos_Matrix * (End_Effector_Obs_Dist - End_Effector_Obs_vec * mini);
	for (int i = 25; i < 33; i++)
	{
		F[i] = Inq_Pos_Constraint(i-25);
		// printf("Inq_Pos_Constraint(%d): %f\n", i-25, F[i]);
	}

	/*
	**	2. Constraint maintenance constraints
	*/
	Eqn_Maint_Matrix(0,0) = sigma_i[0] * sigma_i_child[0];
	Eqn_Maint_Matrix(1,1) = sigma_i[0] * sigma_i_child[0];
	Eqn_Maint_Matrix(2,2) = sigma_i[0] * sigma_i_child[0];
	Eqn_Maint_Matrix(3,3) = sigma_i[0] * sigma_i_child[0];

	Eqn_Maint_Matrix(4,4) = sigma_i[1] * sigma_i_child[1];
	Eqn_Maint_Matrix(5,5) = sigma_i[1] * sigma_i_child[1];
	Eqn_Maint_Matrix(6,6) = sigma_i[1] * sigma_i_child[1];
	Eqn_Maint_Matrix(7,7) = sigma_i[1] * sigma_i_child[1];

	Eqn_Maint_Matrix(8,8) = sigma_i[2] * sigma_i_child[2];
	Eqn_Maint_Matrix(9,9) = sigma_i[2] * sigma_i_child[2];
	Eqn_Maint_Matrix(10,10) = sigma_i[3] * sigma_i_child[3];
	Eqn_Maint_Matrix(11,11) = sigma_i[3] * sigma_i_child[3];

	dlib::matrix<double,16,1> offset_End_Effector_Pos, offset_End_Effector_Vel;
	for (int i = 0; i < 16; i++)
	{
		offset_End_Effector_Pos(i) = ref_End_Effector_Pos(i) - End_Effector_Pos(i);
		offset_End_Effector_Vel(i) = ref_End_Effector_Vel(i) - End_Effector_Vel(i);
	}
	//
	// std::cout<<offset_End_Effector_Pos<<endl;
	// std::cout<<offset_End_Effector_Vel<<endl;
	// std::cout << "Eqn_Maint_Matrix"<< endl;
	// std::cout << Eqn_Maint_Matrix << endl;

	offset_End_Effector_Pos = Eqn_Maint_Matrix * offset_End_Effector_Pos;
	offset_End_Effector_Vel = Eqn_Maint_Matrix * offset_End_Effector_Vel;
	// std::cout<<offset_End_Effector_Pos<<endl;
	// std::cout<<offset_End_Effector_Vel<<endl;
	for (int i = 33; i < 49; i++)
	{
		F[i] = offset_End_Effector_Pos(i-33);
		// printf("offset_End_Effector_Pos(%d): %f\n", i - 33, F[i]);

	}
	for (int i = 49; i < 65; i++)
	{
		F[i] = offset_End_Effector_Vel(i-49);
		// printf("offset_End_Effector_Vel(%d): %f\n", i - 49, F[i]);
	}

	return 0;
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

dlib::matrix<double> End_Effector_Obs_Dist_Fn(dlib::matrix<double,16,1> &End_Effector_Pos, std::vector<double> &Contact_Ind)
{
	dlib::matrix<double,8,1> End_Effector_Dist;
	int Obs_Choice_Ind; char char_sym;
	for (int i = 0; i < 8; i++)
	{
		std::vector<double> r_Pos;
		r_Pos.push_back(End_Effector_Pos(2*i));
		r_Pos.push_back(End_Effector_Pos(2*i+1));

		if (i<4)
		{
			char_sym = 'y';
		}
		else
		{
			char_sym= 'z';
		}

		End_Effector_Dist(i) = Obs_Dist_Fn(r_Pos, Envi_Map, Obs_Choice_Ind, char_sym);
		Contact_Ind[i] = Obs_Choice_Ind;
	}
	return End_Effector_Dist;
}

dlib::matrix<double> End_Effector_Obs_Dist_Fn(dlib::matrix<double,16,1> &End_Effector_Pos)
{
	dlib::matrix<double,8,1> End_Effector_Dist;
	int Obs_Choice_Ind; char char_sym;
	for (int i = 0; i < 8; i++)
	{
		std::vector<double> r_Pos;
		r_Pos.push_back(End_Effector_Pos(2*i));
		r_Pos.push_back(End_Effector_Pos(2*i+1));

		if (i<4)
		{
			char_sym = 'y';
		}
		else
		{
			char_sym= 'z';
		}

		End_Effector_Dist(i) = Obs_Dist_Fn(r_Pos, Envi_Map, Obs_Choice_Ind, char_sym);
	}
	return End_Effector_Dist;
}

void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,16,1> &End_Effector_Pos, dlib::matrix<double,16,1> &End_Effector_Vel)
{
	std::vector<double> rA = Ang_Pos_fn(StateNDot_Init_i, "rA");
	std::vector<double> rB = Ang_Pos_fn(StateNDot_Init_i, "rB");
	std::vector<double> rC = Ang_Pos_fn(StateNDot_Init_i, "rC");
	std::vector<double> rD = Ang_Pos_fn(StateNDot_Init_i, "rD");
	std::vector<double> rE = Ang_Pos_fn(StateNDot_Init_i, "rE");
	std::vector<double> rF = Ang_Pos_fn(StateNDot_Init_i, "rF");
	std::vector<double> rT = Ang_Pos_fn(StateNDot_Init_i, "rT");
	std::vector<double> rCOM = Ang_Pos_fn(StateNDot_Init_i, "rCOM");

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
	End_Effector_Pos(12) = rT[0];
	End_Effector_Pos(13) = rT[1];
	End_Effector_Pos(14) = rCOM[0];
	End_Effector_Pos(15) = rCOM[1];

	// These are the positions of the intermediate joints and the head joint
	// std::vector<double> rH = Ang_Pos_fn(StateNDot_Init_i, "rH");
	// std::vector<double> rK = Ang_Pos_fn(StateNDot_Init_i, "rK");
	// std::vector<double> rM = Ang_Pos_fn(StateNDot_Init_i, "rM");
	// std::vector<double> rN = Ang_Pos_fn(StateNDot_Init_i, "rN");

	// std::vector<double> rCOM = Ang_Pos_fn(StateNDot_Init_i, "rCOM");

	std::vector<double> vA = Ang_Vel_fn(StateNDot_Init_i, "vA");
	std::vector<double> vB = Ang_Vel_fn(StateNDot_Init_i, "vB");
	std::vector<double> vC = Ang_Vel_fn(StateNDot_Init_i, "vC");
	std::vector<double> vD = Ang_Vel_fn(StateNDot_Init_i, "vD");
	std::vector<double> vE = Ang_Vel_fn(StateNDot_Init_i, "vE");
	std::vector<double> vF = Ang_Vel_fn(StateNDot_Init_i, "vF");
	std::vector<double> vT = Ang_Vel_fn(StateNDot_Init_i, "vT");
	std::vector<double> vCOM = Ang_Vel_fn(StateNDot_Init_i, "vCOM");

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
	End_Effector_Vel(12) = vT[0];
	End_Effector_Vel(13) = vT[1];
	End_Effector_Vel(14) = vCOM[0];
	End_Effector_Vel(15) = vCOM[1];
}
