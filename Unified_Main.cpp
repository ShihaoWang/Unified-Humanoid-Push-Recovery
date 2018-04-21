#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include "Unified_Header.h"
#include "snfilewrapper.hh"
#include "snopt.hh"
#include "snoptProblem.hh"

using namespace std;
int Path_Find_Flag =0;
Unified_Structure_P Structure_P;
int main( int argc, char **argv)
{
	// This function is used to generate the contact graph
	std::vector<int> sigma_init(4);
	cout<<sigma_init[0]<<endl;
	// sigma_init is used initialize the robot configuration at the initial time
	sigma_init[0] = 1; 	sigma_init[1] = 0;  	sigma_init[2] = 0;  	sigma_init[3] = 0;
	Structure_P.sigma_i = sigma_init;
	std::vector<double> x0_init = Default_Init(sigma_init, Structure_P,0);
	Robot_StateNDot Init_Opt_vec(x0_init);
	Robot_Plot_fn(Init_Opt_vec);

	// The root node initialization
	Tree_Node Root_Node;
	Root_Node.time = 0.0;
	Root_Node.Kinetic_Energy = Kinetic_Energy_fn(Init_Opt_vec);
	Root_Node.Node_Number = 0;

	Root_Node.StateNDot_Str = StateVec2StateNDot(x0_init);
	Root_Node.Par_Node = NULL;
	Root_Node.sigma_i = sigma_init;
	Add_Node(Root_Node);

	while(1)
	{
		// The first step is to do the Node_Expansion_fn
		// The node chosen to be expanded is the node from Frontier with the
		// minimum kinetic energy

		Tree_Node  Current_Node;
		Current_Node = Pop_Node();
		cout<<"Haha"<<endl;

	}




	return 0;
}
