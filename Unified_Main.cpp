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
	std::vector<double> sigma_init;
	// sigma_init is used initialize the robot contact at the initial time

	ifstream Sigma_Init_File;              // This is to read the initial angle and angular velocities
	Sigma_Init_File.open("sigma_init.txt");
	if(Sigma_Init_File.is_open())
	{
		double data_each_line = 0.0;

		while(Sigma_Init_File>>data_each_line)
		{
			sigma_init.push_back(data_each_line);
		}
		Sigma_Init_File.close();
	}
	else
	{
		printf("Unable to open sigma_init.txt file!\n");
	}
	std::vector<double> Robot_Init_vec = Default_Init(sigma_init);
	Robot_StateNDot StateNDot_Init_Opt(Robot_Init_vec);
	Robot_Plot_fn(StateNDot_Init_Opt);

	// After the robot state initialization, the next job is to conduct the multi-contact staiblization strategy: the root node initialization
	Tree_Node Root_Node;
	Node_UpdateNCon(Root_Node, StateNDot_Init_Opt, sigma_init);
	Tree_Node Node_i, Node_i_child;
	int Self_Opt_Flag, Nodes_Opt_Flag;
	while(Frontier_Nodes.size()>0)
	{
		/**
		* For the current node, first is the Node_Self_Opt to optimize a motion while maintain the current mode
		* if this does not work, then expand the current node into the adjacent nodes then do the Nodes_Connectivity_Opt
		*/
		Node_i = Pop_Node(); 
		int Flag;
		std::vector<double> test;
		test = Nodes_Optimization_fn(Node_i, Node_i, Flag);



	}
	return 0;
}
