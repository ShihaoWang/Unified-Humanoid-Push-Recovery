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

dlib::matrix<double> xlow_vec;						dlib::matrix<double> xupp_vec;
dlib::matrix<double> ctrl_low_vec;					dlib::matrix<double> ctrl_upp_vec;
dlib::matrix<double>  Envi_Map;						dlib::matrix<double> Envi_Map_Normal, Envi_Map_Tange; // The Normal and tangential vector of the plane
/**
 * Some global values are defined
 * Description
 */
double mini = 0.05;			int Grids = 10;			double mu = 0.35;
double Time_Seed = 0.5; 							// This value will be adaptively changed to formulate an optimal solution
std::vector<Tree_Node_Ptr> All_Nodes;				// All nodes are here!
std::vector<Tree_Node_Ptr> Children_Nodes;			// All children nodes!
std::vector<Tree_Node_Ptr> Frontier_Nodes;			// Only Frontier ndoes!
std::vector<double> Frontier_Nodes_Cost;		    // The kinetic energy of each nodes


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
	double T_tot, T;										dlib::matrix<double> StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj;
	// Opt_Seed = Opt_Soln_Load();
	Opt_Seed_Unzip(Opt_Seed, T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj);
	T = T_tot/(Grids - 1)*1.0;
	Tree_Node Node_i, Node_i_child;							Node_i = Structure_P.Node_i;								Node_i_child = Structure_P.Node_i_child;
	std::vector<double> Robotstate_Initial_Vec = StateNDot2StateVec(Node_i.Node_StateNDot);
	int Self_Opt_Flag;										std::vector<double> sigma, sigma_i, sigma_i_child;
	sigma_i = Node_i.sigma;									sigma_i_child = Node_i_child.sigma;
	std::vector<double> sigma_trans, sigma_goal;			Sigma_TransNGoal(sigma_i, sigma_i_child, sigma_trans, sigma_goal, Self_Opt_Flag);
	double KE_ref;
	if(Self_Opt_Flag ==1)
	{// In this way, it is the self_optimization`
		 KE_ref = 0.01;
	}
	else{	KE_ref = 0.25 * Node_i.KE;}

	// Objective function initialization
	ObjNConstraint_Val.push_back(0);
	ObjNConstraint_Type.push_back(1);

	// 1. The first constraint is to make sure that the initial condition matches the given initial condition
	dlib::matrix<double> StateNDot_Traj_1st_colm, Robotstate_Initial_VecDlib, Matrix_result;
	StateNDot_Traj_1st_colm = dlib::colm(StateNDot_Traj, 0);			Robotstate_Initial_VecDlib = StateVec2DlibMatrix_fn(Robotstate_Initial_Vec);
	// cout<<StateNDot_Traj_1st_colm<<endl;								cout<<Robotstate_Initial_VecDlib<<endl;
	Matrix_result = Quadratic_Minus(StateNDot_Traj_1st_colm, Robotstate_Initial_VecDlib);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	// THe constraint will be added term by term
	dlib::matrix<double> Robostate_Dlib_Front, Robostate_Dlib_Back, Ctrl_Front, Ctrl_Back, Contact_Force_Front, Contact_Force_Back, Robotstate_Mid_Acc;

	// 2. The second constraint is Dynamics constraints
	dlib::matrix<double> D_q_Front, B_q_Front, 	C_q_qdot_Front, Jac_Full_Front, Jac_Full_Trans_Front;
	dlib::matrix<double> D_q_Mid, 	B_q_Mid, 	C_q_qdot_Mid, 	Jac_Full_Mid, 	Jac_Full_Trans_Mid;
	dlib::matrix<double> D_q_Back, 	B_q_Back, 	C_q_qdot_Back, 	Jac_Full_Back, 	Jac_Full_Trans_Back;
	dlib::matrix<double> Dynamics_LHS, Dynamics_RHS;
	Robot_StateNDot Robot_StateNDot_Front, Robot_StateNDot_Mid, Robot_StateNDot_Back;
	for (int i = 0; i < Grids-1; i++)
	{	// Get the robot state, ctrl, and contact force at the front and end of each segment
		Robostate_Dlib_Front = dlib::colm(StateNDot_Traj, i);						Robostate_Dlib_Back = dlib::colm(StateNDot_Traj, i+1);
		Ctrl_Front = dlib::colm(Ctrl_Traj,i);										Ctrl_Back = dlib::colm(Ctrl_Traj,i+1);
		Contact_Force_Front = dlib::colm(Contact_Force_Traj,i);						Contact_Force_Back = dlib::colm(Contact_Force_Traj,i+1);

		// cout<<Robostate_Dlib_Front<<endl;											cout<<Robostate_Dlib_Back<<endl;
		// cout<<Ctrl_Front<<endl;														cout<<Ctrl_Back<<endl;
		// cout<<Contact_Force_Front<<endl;											cout<<Contact_Force_Back<<endl;

		// Compute the Dynamics matrices at Front and Back
		Robot_StateNDot_Front = DlibRobotstate2StateNDot(Robostate_Dlib_Front);		Dynamics_Matrices(Robot_StateNDot_Front, D_q_Front, B_q_Front, C_q_qdot_Front, Jac_Full_Front);
		Robot_StateNDot_Back = DlibRobotstate2StateNDot(Robostate_Dlib_Back);		Dynamics_Matrices(Robot_StateNDot_Back, D_q_Back, B_q_Back, C_q_qdot_Back, Jac_Full_Back);

		Robot_StateNDot_MidNAcc(T, Robot_StateNDot_Front, Robot_StateNDot_Back, Ctrl_Front, Ctrl_Back, Contact_Force_Front, Contact_Force_Back, Robot_StateNDot_Mid, Robotstate_Mid_Acc, ObjNConstraint_Val, ObjNConstraint_Type);
		Dynamics_Matrices(Robot_StateNDot_Mid, D_q_Mid, B_q_Mid, C_q_qdot_Mid, Jac_Full_Mid);		Jac_Full_Trans_Mid = dlib::trans(Jac_Full_Mid);
		Dynamics_LHS = D_q_Mid * Robotstate_Mid_Acc + C_q_qdot_Mid;
		Dynamics_RHS = Jac_Full_Trans_Mid * (0.5 * Contact_Force_Front + 0.5 * Contact_Force_Back) + B_q_Mid * (0.5 * Ctrl_Front + 0.5 * Ctrl_Back);
		// cout<<Dynamics_LHS<<endl;													cout<<Dynamics_RHS<<endl;
		// Matrix_result = Quadratic_Minus(Dynamics_LHS, Dynamics_RHS);
		Matrix_result = Dynamics_LHS - Dynamics_RHS;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	}
	// 3. Complementarity constraints: Distance!
	dlib::matrix<double,12,1> End_Effector_Pos, End_Effector_Vel;					dlib::matrix<double,6,1> End_Effector_Dist;
	dlib::matrix<double> Robostate_Dlib_i;											Robot_StateNDot Robot_StateNDot_i;
	std::vector<int> End_Effector_Obs(6);
	dlib::matrix<double> Eqn_Pos_Matrix, Ineqn_Pos_Matrix, Eqn_Vel_Matrix;
	std::vector<double> sigma_temp;
	for (int i = 0; i < Grids; i++) {
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);							Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
		End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);

		if(i<Grids-1){	sigma = sigma_trans;}
		else{			sigma = sigma_goal;}

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
	}

	// 4. Complementarity constraints: Contact Force!
	dlib::matrix<double> Contact_Force_Complem_Matrix, Contact_Force_i;
	// Since up to here, the constraint version works, now we will try the variable bounds version
	for (int i = 0; i < Grids; i++) {
		std::vector<double> sigma_temp;				Contact_Force_i = dlib::colm(Contact_Force_Traj,i);
		sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);			sigma_temp.push_back(!sigma[0]);				sigma_temp.push_back(!sigma[0]);
		sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);			sigma_temp.push_back(!sigma[1]);				sigma_temp.push_back(!sigma[1]);
		sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[2]);			sigma_temp.push_back(!sigma[3]);				sigma_temp.push_back(!sigma[3]);
		Contact_Force_Complem_Matrix = Diag_Matrix_fn(sigma_temp);								Matrix_result = Contact_Force_Complem_Matrix * Contact_Force_i;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);}

	// 5. Contact force feasibility constraints: 1. Normal force should be positive and 2. the friction cone constraint has to be satisfied
	double Contact_Force_i_x, Contact_Force_i_y;												std:vector<double> Normal_Force, Tange_Force;
	double Normal_Force_1, Normal_Force_2, Normal_Force_3, Normal_Force_4;						double Tange_Force_1, Tange_Force_2, Tange_Force_3, Tange_Force_4;
	for (int i = 0; i < Grids; i++) {
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);										Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);			End_Effector_Obs_Dist_Fn(End_Effector_Pos, End_Effector_Dist, End_Effector_Obs);
		Contact_Force_i = dlib::colm(Contact_Force_Traj,i);
		for (int j = 0; j < Contact_Force_i.nr()/2; j++) {
			Contact_Force_i_x = Contact_Force_i(2*i);											Contact_Force_i_y = Contact_Force_i(2*i+1);
			Normal_Force.push_back(Contact_Force_i_x * Envi_Map_Normal(End_Effector_Obs[j],0) + Contact_Force_i_y * Envi_Map_Normal(End_Effector_Obs[j],1));
			Tange_Force.push_back(Contact_Force_i_x * Envi_Map_Tange(End_Effector_Obs[j],0) + Contact_Force_i_y * Envi_Map_Tange(End_Effector_Obs[j],1));}
		for (int j = 0; j < Contact_Force_i.nr()/2; j++) {
			ObjNConstraint_Val.push_back(Normal_Force[j]);										ObjNConstraint_Type.push_back(1);}
		Normal_Force_1 = Normal_Force[0] + Normal_Force[1];		Normal_Force_2 = Normal_Force[2] + Normal_Force[3];		Normal_Force_3 = Normal_Force[4];	Normal_Force_4 = Normal_Force[5];
		Tange_Force_1 = Tange_Force[0] + Tange_Force[1];		Tange_Force_2 = Tange_Force[2] + Tange_Force[3];		Tange_Force_3 = Tange_Force[4];		Tange_Force_4 = Tange_Force[5];
		ObjNConstraint_Val.push_back(Normal_Force_1 * Normal_Force_1 * mu * mu - Tange_Force_1 * Tange_Force_1);		ObjNConstraint_Type.push_back(1);
		ObjNConstraint_Val.push_back(Normal_Force_2 * Normal_Force_2 * mu * mu - Tange_Force_2 * Tange_Force_2);		ObjNConstraint_Type.push_back(1);
		ObjNConstraint_Val.push_back(Normal_Force_3 * Normal_Force_3 * mu * mu - Tange_Force_3 * Tange_Force_3);		ObjNConstraint_Type.push_back(1);
		ObjNConstraint_Val.push_back(Normal_Force_4 * Normal_Force_4 * mu * mu - Tange_Force_4 * Tange_Force_4);		ObjNConstraint_Type.push_back(1);}
	//
	// 6. Contact maintenance constraint: the previous unchanged active constraint have to be satisfied
	dlib::matrix<double> End_Effector_Pos_ref;													End_Effector_Pos_ref = Node_i.End_Effector_Pos;
	std::vector<double> sigma_maint(12);														dlib::matrix<double> Matrix_Minus_result, End_Effector_PosDlib;
	for (int i = 0; i < Grids; i++) {
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i);										Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		End_Effector_PosNVel(Robot_StateNDot_i, End_Effector_Pos, End_Effector_Vel);
		sigma_maint[0] = sigma_i[0] * sigma_i_child[0];		sigma_maint[1] = sigma_i[0] * sigma_i_child[0];
		sigma_maint[2] = sigma_i[0] * sigma_i_child[0];		sigma_maint[3] = sigma_i[0] * sigma_i_child[0];
		sigma_maint[4] = sigma_i[1] * sigma_i_child[1];		sigma_maint[5] = sigma_i[1] * sigma_i_child[1];
		sigma_maint[6] = sigma_i[1] * sigma_i_child[1];		sigma_maint[7] = sigma_i[1] * sigma_i_child[1];
		sigma_maint[8] = sigma_i[2] * sigma_i_child[2];		sigma_maint[9] = sigma_i[2] * sigma_i_child[2];
		sigma_maint[10] = sigma_i[3] * sigma_i_child[3];	sigma_maint[11] = sigma_i[3] * sigma_i_child[3];
		dlib::matrix<double> Maint_Matrix = Diag_Matrix_fn(sigma_maint);
		End_Effector_PosDlib = End_Effector_Pos;
		Matrix_Minus_result = End_Effector_PosDlib - End_Effector_Pos_ref;
		Matrix_result = Maint_Matrix * Matrix_Minus_result;
		ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	}
	// 7. Objective function value: the first value is fixed
	double KE_i;					std::vector<double> KE_tot;
	for (int i = 0; i < Grids-1; i++) {
		Robostate_Dlib_i = dlib::colm(StateNDot_Traj, i+1);										Robot_StateNDot_i = DlibRobotstate2StateNDot(Robostate_Dlib_i);
		KE_i = Kinetic_Energy_fn(Robot_StateNDot_i);
		KE_tot.push_back(KE_i);
	}
	ObjNConstraint_Val[0] = KE_Variation_fn(KE_tot, T);
	// ObjNConstraint_Val[0] = Traj_Variation(StateNDot_Traj);
	// ObjNConstraint_Val[0] = Torque_Sum(Ctrl_Traj);
	// ObjNConstraint_Val[0] = Joint_Velocity_Sum(StateNDot_Traj);
	ObjNConstraint_Val.push_back(KE_ref - KE_i);
	ObjNConstraint_Type.push_back(1);
}
double Torque_Sum(dlib::matrix<double> &Ctrl_Traj)
{
	double Torque_Sum_Val = 0.0;
	dlib::matrix<double> Ctrl_Traj_i;
	for (int i = 0; i < Ctrl_Traj.nc(); i++) {
		Ctrl_Traj_i = dlib::colm(Ctrl_Traj, i);
		for (int j = 0; j < Ctrl_Traj_i.nr(); j++) {
			Torque_Sum_Val = Torque_Sum_Val + Ctrl_Traj_i(j) * Ctrl_Traj_i(j);
		}
	}
	return Torque_Sum_Val;
}
double Joint_Velocity_Sum(dlib::matrix<double> &StateNDot_Traj)
{
	// This function is used to calcualte the variation of the state trajectory
	dlib::matrix<double> StateNDot_Traj_i;
	double Traj_Variation_Val = 0.0;
	for (int i = 0; i < StateNDot_Traj.nc(); i++) {
		StateNDot_Traj_i = dlib::colm(StateNDot_Traj, i);
		for (int j = 0; j < StateNDot_Traj_i.nr()/2; j++) {
			Traj_Variation_Val = Traj_Variation_Val +  StateNDot_Traj_i(j+13) * StateNDot_Traj_i(j+13);
		}
	}
	// Traj_Variation_Val = Traj_Variation_Val * Traj_Variation_Val;
	return Traj_Variation_Val;
}
double Traj_Variation(dlib::matrix<double> &StateNDot_Traj)
{
	// This function is used to calcualte the variation of the state trajectory
	dlib::matrix<double> StateNDot_Traj_Front, StateNDot_Traj_Back, Matrix_result;
	double Traj_Variation_Val = 0.0;
	for (int i = 0; i < StateNDot_Traj.nc()-1; i++) {
		StateNDot_Traj_Front = dlib::colm(StateNDot_Traj, i);
		StateNDot_Traj_Back = dlib::colm(StateNDot_Traj, i+1);
		Matrix_result = StateNDot_Traj_Front - StateNDot_Traj_Back;
		for (int j = 0; j < Matrix_result.nr()/2; j++) {
			Traj_Variation_Val = Traj_Variation_Val +  Matrix_result(j) * Matrix_result(j);
		}
	}
	Traj_Variation_Val = Traj_Variation_Val * Traj_Variation_Val;
	return Traj_Variation_Val;
}
Robot_StateNDot DlibRobotstate2StateNDot(dlib::matrix<double> &DlibRobotstate)
{
	// This function is used to convert the dlib matrix robot state to Robot_StateNDot type
	std::vector<double> Robot_StateNDot_vec;
	for (int i = 0; i < DlibRobotstate.nr(); i++) {
		Robot_StateNDot_vec.push_back(DlibRobotstate(i));}
	Robot_StateNDot Robot_StateNDot_i(Robot_StateNDot_vec);
	return Robot_StateNDot_i;
}
void Robot_StateNDot_MidNAcc(double T, const Robot_StateNDot &Robot_StateNDot_Front, const Robot_StateNDot &Robot_StateNDot_Back, const dlib::matrix<double> &Ctrl_Front, const dlib::matrix<double> &Ctrl_Back, const dlib::matrix<double> &Contact_Force_Front, const dlib::matrix<double> &Contact_Force_Back, Robot_StateNDot &Robot_StateNDot_Mid, dlib::matrix<double> &Robotstate_Mid_Acc,std::vector<double> &ObjNConstraint_Val, std::vector<double> &ObjNConstraint_Type)
{
	std::vector<double> Robotstate_Vec_Front, Robotstate_Vec_Back;
	dlib::matrix<double> D_q_Front, B_q_Front, 	C_q_qdot_Front, Jac_Full_Front, Jac_Full_Trans_Front, Acc_Front;
	dlib::matrix<double> D_q_Back, 	B_q_Back, 	C_q_qdot_Back, 	Jac_Full_Back, 	Jac_Full_Trans_Back, Acc_Back;

	Dynamics_Matrices(Robot_StateNDot_Front, D_q_Front, B_q_Front, C_q_qdot_Front, Jac_Full_Front);
	Jac_Full_Trans_Front = dlib::trans(Jac_Full_Front);
	Dynamics_Matrices(Robot_StateNDot_Back, D_q_Back, B_q_Back, C_q_qdot_Back, Jac_Full_Back);
	Jac_Full_Trans_Back = dlib::trans(Jac_Full_Back);

	Acc_Front = dlib::inv(D_q_Front) * (Jac_Full_Trans_Front * Contact_Force_Front + B_q_Front * Ctrl_Front - C_q_qdot_Front);
	Acc_Back =  dlib::inv(D_q_Back) *  (Jac_Full_Trans_Back *  Contact_Force_Back +  B_q_Back *  Ctrl_Back -  C_q_qdot_Back);

	// For each variable in the state, we will calculate the Cubic spline coefficients and the middle state N Acc
	std::vector<double> Robotstate_Vec_Mid(26);										Robotstate_Mid_Acc = dlib::zeros_matrix<double>(13,1);
	Robotstate_Vec_Front = StateNDot2StateVec(Robot_StateNDot_Front);
	Robotstate_Vec_Back = StateNDot2StateVec(Robot_StateNDot_Back);

	double x_init, x_end, xdot_init, xdot_end, xddot_init, xddot_end;				std::vector<double> CubicSpline_Coeff;

	// 1. Calculate the Robotstate_Mid_Acc Pos
	for (int i = 0; i < 13; i++) {
		x_init = Robotstate_Vec_Front[i];				x_end = Robotstate_Vec_Back[i];
		xdot_init = Robotstate_Vec_Front[i+13];			xdot_end = Robotstate_Vec_Back[i+13];
		CubicSpline_Coeff = CubicSpline_Coeff_fn(T, x_init, x_end, xdot_init, xdot_end);
		Robotstate_Vec_Mid[i] = CubicSpline_Evaluation_fn(CubicSpline_Coeff, 0.5);
	}
	// 2. Calculate the Robotstate_Mid_Acc Vel and Acc
	for (int i = 0; i < 13; i++) {
		xdot_init = Robotstate_Vec_Front[i+13];			xdot_end = Robotstate_Vec_Back[i+13];
		xddot_init = Acc_Front(i);						xddot_end = Acc_Back(i);
		CubicSpline_Coeff = CubicSpline_Coeff_fn(T, xdot_init, xdot_end, xddot_init, xddot_end);
		Robotstate_Vec_Mid[i+13] = CubicSpline_Evaluation_fn(CubicSpline_Coeff, 0.5);
		Robotstate_Mid_Acc(i) = CubicSpline_1stOrder_Evaluation_fn(CubicSpline_Coeff, 0.5, T);
	}
	// cout<<Acc_Front<<endl;			cout<<Acc_Back<<endl;			cout<<Robotstate_Mid_Acc<<endl;
	Robot_StateNDot_Mid = StateVec2StateNDot(Robotstate_Vec_Mid);
	for (int i = 1; i < 10; i++) {
		ObjNConstraint_Val.push_back(Acc_Front(i+3));
		ObjNConstraint_Type.push_back(2);
		ObjNConstraint_Val.push_back(Robotstate_Mid_Acc(i+3));
		ObjNConstraint_Type.push_back(2);
		ObjNConstraint_Val.push_back(Acc_Back(i+3));
		ObjNConstraint_Type.push_back(2);
	}
	// ObjNConstraint_ValNType_Update(Robotstate_Mid_Acc, ObjNConstraint_Val, ObjNConstraint_Type, 2);
}
void Quadratic_Angular_Sum_Cal(std::vector<double> &Robot_Vel,double &Quadratic_Angular_Sum)
{
	for (int i = 0; i < Robot_Vel.size(); i++) {
		Quadratic_Angular_Sum = Quadratic_Angular_Sum + Robot_Vel[i] * Robot_Vel[i];
	}
}
dlib::matrix<double> StateNDot_ref_fn(std::vector<double> &Robot_Config_i, std::vector<double> &Robot_Velocity_i)
{
	dlib::matrix<double> StateNDot_ref_i;
	StateNDot_ref_i = dlib::zeros_matrix<double>(26,1);
	for (int i = 0; i < 13; i++) {
		StateNDot_ref_i(i) = Robot_Config_i[i];
		StateNDot_ref_i(i+13) = Robot_Velocity_i[i];}
	return StateNDot_ref_i;
}
double KE_Variation_fn(std::vector<double> &KE_tot, double T)
{	double KE_Variation = 0.0;
	// KE_Variation = KE_tot[KE_tot.size()-1];
	for (int i = 0; i < KE_tot.size(); i++)
	{
	// 	// KE_Variation = KE_Variation + 0.0*KE_tot[i];
	// 	// KE_Variation = KE_Variation + KE_tot[KE_tot.size()-1];
		KE_Variation = KE_Variation + KE_tot[i];
	}
	return KE_Variation;
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

	dlib::matrix<double> End_Effector_Pos_ref, End_Effector_Vel_ref;
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
	dlib::matrix<double> Matrix_result, Matrix_Minus_result;
	dlib::matrix<double> End_Effector_PosDlib;
	End_Effector_PosDlib = End_Effector_Pos;

	Matrix_Minus_result = Quadratic_Minus(End_Effector_PosDlib, End_Effector_Pos_ref);
	Matrix_result = Maint_Matrix * End_Effector_Pos_ref;
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);

	Matrix_result = Maint_Matrix * (End_Effector_Vel);
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
	dlib::matrix<double> Matrix_result, Dynamics_LHS, Dynamics_RHS;
	Dynamics_LHS = D_q * Acc + C_q_qdot;
	Dynamics_RHS = Jac_Full_Trans * Contact_Force + B_q * Ctrl;
	Matrix_result = Quadratic_Minus(Dynamics_LHS, Dynamics_RHS);
	ObjNConstraint_ValNType_Update(Matrix_result, ObjNConstraint_Val, ObjNConstraint_Type, 0);
	for (int i = 4; i < Acc.nr(); i++) {
		ObjNConstraint_Val.push_back(9 - Acc(i) * Acc(i));
		ObjNConstraint_Type.push_back(1);
	}
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
	Ctrl_i = dlib::zeros_matrix<double>(10,1);				Contact_Force_i = dlib::zeros_matrix<double>(12,1);
	for (int i = 0; i < 10; i++) {
		x_a = Ctrl_Coeff(2*i, Grid_Ind);
		x_b = Ctrl_Coeff(2*i+1, Grid_Ind);
		Ctrl_i(i) = x_a * s + x_b;}
	for (int i = 0; i < 12; i++) {
		x_a = Contact_Force_Coeff(2*i, Grid_Ind);
		x_b = Contact_Force_Coeff(2*i+1, Grid_Ind);
		Contact_Force_i(i) = x_a * s + x_b;}
}
int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child)
{	// This function will optimize the joint trajectories to minimize the robot kinetic energy while maintaining a smooth transition fashion
	// However, the constraint will be set to use the direct collocation method
	Structure_P.Node_i = Node_i;		Structure_P.Node_i_child = Node_i_child;
	// A solution will be generated iteratively
	int Opt_Flag = 0;
	std::vector<double> Opt_Seed;		double Constraint_Vio;
	for (int i = 0; i < 15; i++) {
		Time_Seed = 0.25 + 1.0 * i * 0.05;
		Opt_Seed = Nodes_Optimization_fn_inner(Node_i, Node_i_child, Constraint_Vio);
		if (Constraint_Vio<0.1)
		{
			Opt_Flag = 1;
			ofstream output_file;
			std::string filename;
			filename.append("Node_");
			filename.append(to_string(Node_i.Node_Index));
			filename.append("to");
			filename.append("Node_");
			filename.append(to_string(Node_i_child.Node_Index));
			filename.append("_Opt_Soln.txt");
			output_file.open(filename, std::ofstream::out);
			for (int i = 0; i < Opt_Seed.size(); i++)
			{
				output_file<<Opt_Seed[i]<<endl;
			}
			output_file.close();
			break;
		}
	}
}
std::vector<double> Nodes_Optimization_fn_inner(Tree_Node &Node_i, Tree_Node &Node_i_child, double &Constraint_Vio)
{
	std::vector<double> Opt_Seed = Seed_Guess_Gene(Node_i, Node_i_child);
	std::vector<double> ObjNConstraint_Val, ObjNConstraint_Type;
	Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	snoptProblem Nodes_Optimization_Pr;                     // This is the name of the Optimization problem for the robot configuration

	integer n = Opt_Seed.size();
	integer neF = ObjNConstraint_Val.size();
	integer lenA  =  n * neF;

	Structure_P.Opt_Val_No = Opt_Seed.size();		Structure_P.ObjNConst_No = ObjNConstraint_Val.size();

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
	}
	xlow[0] = 0.5;		xupp[0] = 1;
	int Index_Count = 1;
	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < 26; j++) {
			xlow[Index_Count] = xlow_vec(j);
			xupp[Index_Count] = xupp_vec(j);
			Index_Count = Index_Count + 1;}}

	for (int i = 0; i < Grids; i++) {
		for (int j = 0; j < 10; j++) {
			xlow[Index_Count] = ctrl_low_vec(j);
			xupp[Index_Count] = ctrl_upp_vec(j);
			Index_Count = Index_Count + 1;}}

	for (int i = 0; i < n; i++) {
		xstate[i] = 0.0;
		x[i] = Opt_Seed[i];  	// Initial guess
	}

	for(int i = 0; i<neF; i++)
	{
		if (ObjNConstraint_Type[i]>1.0)
		{
			Flow[i] = -3.0;
			Fupp[i] = 3.0;
		}
		else
		{
			Flow[i] = 0.0;
			if(ObjNConstraint_Type[i]>0.0)
			{
				Fupp[i] = Inf;
			}
			else
			{
				Fupp[i] = 0.0;
			}
		}
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
	Nodes_Optimization_Pr.setIntParameter("Major iterations limit", 150);
	Nodes_Optimization_Pr.setIntParameter("Minor iterations limit", 2000000);
	Nodes_Optimization_Pr.setIntParameter("Iterations limit", 2000000);
	Nodes_Optimization_Pr.setIntParameter("setFeaTol", 1e-4);
	Nodes_Optimization_Pr.setIntParameter("setOptTol", 1e-3);
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
	std::vector<double> ObjNConstraint_Final_Val, ObjNConstraint_Final_Type;
	Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Final_Val, ObjNConstraint_Final_Type);
	double Constraint_Vio_Temp;
	Constraint_Vio = 0.0;
	for (int i = 0; i < ObjNConstraint_Final_Val.size(); i++) {
		Constraint_Vio_Temp = abs(ObjNConstraint_Final_Val[i] * (!ObjNConstraint_Final_Type[i]));
		if (Constraint_Vio_Temp>Constraint_Vio)
		{Constraint_Vio = Constraint_Vio_Temp;}}
	return Opt_Seed;
}
void Contact_Force_Bounds(std::vector<double> &sigma, std::vector<double> &Contact_Force_Status_i)
{
	Contact_Force_Status_i[0] = !sigma[0];				Contact_Force_Status_i[1] = !sigma[0];
	Contact_Force_Status_i[2] = !sigma[0];				Contact_Force_Status_i[3] = !sigma[0];

	Contact_Force_Status_i[4] = !sigma[1];				Contact_Force_Status_i[5] = !sigma[1];
	Contact_Force_Status_i[6] = !sigma[1];				Contact_Force_Status_i[7] = !sigma[1];

	Contact_Force_Status_i[8] = !sigma[2];				Contact_Force_Status_i[9] = !sigma[2];
	Contact_Force_Status_i[10] = !sigma[3];				Contact_Force_Status_i[11] = !sigma[3];
}
int Nodes_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
		 integer    iu[],    integer *leniu,
		 doublereal ru[],    integer *lenru )
{	 std::vector<double> Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type;
	const int Opt_Val_No = Structure_P.Opt_Val_No ;				const int ObjNConst_No = Structure_P.ObjNConst_No;
	for (int i = 0; i < Opt_Val_No; i++)
	{
		Opt_Seed.push_back(x[i]);
	}

	Nodes_Optimization_ObjNConstraint(Opt_Seed, ObjNConstraint_Val, ObjNConstraint_Type);
	for (int i = 0; i < ObjNConst_No; i++)
	{
		F[i] = ObjNConstraint_Val[i];
	}

	return 0;
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
void Pos_Vel_Acc_VelfromPos_fromStateNdot_Coeff(double T, dlib::matrix<double> &StateNDot_Coeff, int Grid_Ind, double s, std::vector<double> &Robot_Config,  std::vector<double> &Robot_Vel, dlib::matrix<double> &Robot_Acc, std::vector<double> &Robot_VelfromPos)
{
	//Here Grid_Ind denotes which Grid we are talking about and s is a value between 0 and 1
	Robot_Acc = dlib::zeros_matrix<double>(13,1);
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
