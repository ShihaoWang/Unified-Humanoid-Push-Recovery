#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <numeric>

#include "Unified_Header.h"
#include "snfilewrapper.hh"
#include "snopt.hh"
#include "snoptProblem.hh"

using namespace std;
Unified_Structure_P Structure_P;

int main( int argc, char **argv)
{
	// This function is used to generate the contact graph
	std::vector<double> sigma_init(4);
	// sigma_init is used initialize the robot contact at the initial time
	sigma_init[0] = 1; 	sigma_init[1] = 0;  	sigma_init[2] = 0;  	sigma_init[3] = 0;
	Structure_P.sigma_i = sigma_init;
	std::vector<double> Robot_Init_vec = Default_Init(sigma_init, 0);
	Robot_StateNDot StateNDot_Init_Opt(Robot_Init_vec);
	// Robot_Plot_fn(StateNDot_Init_Opt);

	// After the robot state initialization, the next job is to conduct the multi-contact staiblization strategy

	// The root node initialization
	Tree_Node Root_Node;
	Root_Node.time = 0.0;
	Root_Node.Kinetic_Energy = Kinetic_Energy_fn(StateNDot_Init_Opt);
	Root_Node.Node_Index_Number = 0;

	Root_Node.Node_StateNDot = StateVec2StateNDot(Robot_Init_vec);
	Root_Node.Par_Node = NULL;
	Root_Node.sigma_i = sigma_init;

	Reference_Dist_Vel_Update(Root_Node);
	Add_Node(Root_Node);

	Tree_Node Current_Node;
	std::vector<double> Current_Node_sigma_i; double Current_Node_sigma_sum;
	while(Frontier_Nodes.size()>0)
	{
		/**
		* For the current node, first is the Node_Self_Opt to optimize a motion while maintain the current mode (the fly in the air case will not be considered)
		* 						if this does not work, then expand the current node into the adjacent nodes and figure out whether there exists a path to reach that contact mode
		* Description
		*/
		int Soln_Flag = 0;
		Current_Node = Pop_Node();
		Current_Node_sigma_i = Current_Node.sigma_i;
		Current_Node_sigma_sum = Current_Node_sigma_i[0] + Current_Node_sigma_i[1] + Current_Node_sigma_i[2] + Current_Node_sigma_i[3];
		if (Current_Node_sigma_sum>0)
		{
			// At least one contact is active in this case
			Soln_Flag = Nodes_Optimization_fn(Current_Node, Current_Node, 0);






			/* code */
		}
		// else
		// {
		// 	// This is the flying in the air case
		//
		// }





	}

	//
	// while(1)
	// {
	// 	// The first step is to do the Node_Expansion_fn
	// 	// The node chosen to be expanded is the node from Frontier with the
	// 	// minimum kinetic energy
	// 	Tree_Node  Node_i, Node_i_child;
	// 	Node_i = Pop_Node();
	// 	Node_Expansion_fn(Node_i,Structure_P);
	//
	// 	while(Node_i.Children_Nodes.size())
	// 	{
	// 		// First is to retrieve the child node from the parent node
	// 		// Node_i_child = *(Node_i.Children_Nodes[Node_i.Children_Nodes.size()-1]);
	// 		Node_i.Children_Nodes.pop_back();
	//
	// 	}
	//
	//
	// 	// The children nodes have been attached into the children ptr
	//
	//
	//
	//
	// }




	return 0;
}
