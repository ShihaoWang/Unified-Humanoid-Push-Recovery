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


void Opt_Seed_Unzip(std::vector<double> &Opt_Seed, double &T_tot, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{
	const int NumOfStateNDot_Traj = 26;
	const int NumOfCtrl_Traj = 10;
	const int NumOfContactForce_Traj = 12;

	T_tot = Opt_Seed[0];			int Opt_Seed_Index = 1;

	StateNDot_Traj = dlib::zeros_matrix<double>(NumOfStateNDot_Traj, Grids);
	Ctrl_Traj = dlib::zeros_matrix<double>(NumOfCtrl_Traj, Grids);
	Contact_Force_Traj = dlib::zeros_matrix<double>(NumOfContactForce_Traj, Grids);

	// 1. Retrieve the StateNDot_Traj matrix
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < NumOfStateNDot_Traj; j++) {
			StateNDot_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<StateNDot_Traj<<endl;
	// 2. Retrieve the control matrix
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < NumOfCtrl_Traj; j++) {
			Ctrl_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Ctrl_Traj<<endl;
	// 3. Retrieve the contact force matrix
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < NumOfContactForce_Traj; j++) {
			Contact_Force_Traj(j,i) = Opt_Seed[Opt_Seed_Index];
			Opt_Seed_Index = Opt_Seed_Index + 1;}}
	// cout<<Contact_Force_Traj<<endl;
	return ;
}
void Opt_Seed_Zip(std::vector<double> &Opt_Seed, dlib::matrix<double> & StateNDot_Traj, dlib::matrix<double> & Ctrl_Traj, dlib::matrix<double> & Contact_Force_Traj)
{	// This function is used to stack the coefficient matrices into a column vector
	for (int i = 0; i < StateNDot_Traj.nc(); i++) {
		for (int j = 0; j < StateNDot_Traj.nr(); j++) {
			Opt_Seed.push_back(StateNDot_Traj(j,i));
		}
	}
	for (int i = 0; i < Ctrl_Traj.nc(); i++) {
		for (int j = 0; j < Ctrl_Traj.nr(); j++) {
			Opt_Seed.push_back(Ctrl_Traj(j,i));
		}
	}
	for (int i = 0; i < Contact_Force_Traj.nc(); i++) {
		for (int j = 0; j < Contact_Force_Traj.nr(); j++) {
			Opt_Seed.push_back(Contact_Force_Traj(j,i));
		}
	}
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

std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child)
{	// This function will generate the spline coefficients needed for the further optimization
	// double T = 0.5;
	double T = Time_Seed;
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
	// cout<<StateNDot_Traj<<endl;			cout<<Ctrl_Traj<<endl;				cout<<Contact_Force_Traj<<endl;
	// The calculation of the coefficients of the control and contact force is easier compared to the robot state due to the assumption of the linear equation
	// Ctrl_Contact_Force_Coeff_fn(Ctrl_Traj, Contact_Force_Traj, Ctrl_Coeff, Contact_Force_Coeff);

	// The final task is to pile them into a single vector
	std::vector<double> Opt_Seed;
	Opt_Seed.push_back(T * (Grids - 1));
	Opt_Seed_Zip(Opt_Seed, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);

	// cout<<StateNDot_Coeff<<endl;				cout<<Ctrl_Coeff<<endl;			cout<<Contact_Force_Coeff<<endl;
	return Opt_Seed;
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
