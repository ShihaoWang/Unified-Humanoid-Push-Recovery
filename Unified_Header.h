#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <dlib/matrix.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include "../../matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std;
// using namespace dlib;

#ifndef USERFUNCTION
#define USERFUNCTION

#include "snoptProblem.hh"
#ifdef __cplusplus
extern "C" {
#endif

    int Seed_Conf_Constraint_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

    int Default_Init_Pr_( integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

#ifdef __cplusplus
}
#endif

#endif

class Robot_StateNDot {
public:
	double rIx, rIy, theta, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;
	double rIxdot, rIydot, thetadot, q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot, q10dot;
    Robot_StateNDot();
    Robot_StateNDot(std::vector<double> &);
};
typedef struct Tree_Node *Tree_Node_Ptr;    // This is a pointer to the Tree_Node
struct Tree_Node
{
    //  The running time to this node from parents
    //  The kinetic energy and the unique index number in the All_Nodes vector
    double time;    double Kinetic_Energy;    int Node_Number;

    Robot_StateNDot StateNDot_Str;
    Tree_Node_Ptr Par_Node;
    std::vector<Tree_Node_Ptr> Children_Nodes;

    std::vector<double> sigma_i;
    std::vector<double> StateNDot_Traj;
    std::vector<double> Ctrl_Traj;

    std::vector<double> rA_ref;         std::vector<double> rB_ref;
    std::vector<double> rC_ref;         std::vector<double> rD_ref;
    std::vector<double> rE_ref;         std::vector<double> rF_ref;
    std::vector<double> vA_ref;         std::vector<double> vB_ref;
    std::vector<double> vC_ref;         std::vector<double> vD_ref;
    std::vector<double> vE_ref;         std::vector<double> vF_ref;
    std::vector<double> vI_ref;
};


class Unified_Structure_P{
public:
    double rIx, rIy, theta, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;
	double rIxdot, rIydot, thetadot, q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot, q10dot;
    double u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;
    Tree_Node Node_i;
    Tree_Node Node_i_child;
    std::vector<double> sigma_i;///        well thisis theonly for thye initialization
    std::vector<double> Opt_Ctrl_LowBd;
    std::vector<double> Opt_Ctrl_UppBd;
    std::vector<double> Opt_Conf_LowBd;
    std::vector<double> Opt_Conf_UppBd;
    std::vector<double> Robot_State_Init;
    int Ctrl_No;    double Tme_Seed;    double eps;
    Unified_Structure_P();
}; 
extern Unified_Structure_P Structure_P;

void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i);
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name);
std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  B_q_fn();
dlib::matrix<double>  C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
double Obs_Dist_fn(std::vector<double> &r_pos, const char* s);
std::vector<double> Default_Init(const std::vector<double> &sigma_i, Unified_Structure_P &P, int Flag);
std::vector<double> StateNDot2StateVec(Robot_StateNDot &Robot_StateNDot_i);
Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec);
double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i);
std::vector<double> sigma_modi(std::vector<double> sigma_ref, int contas_ind, int AddOrRet);
int Node_Expansion_fn(Tree_Node &Cur_Node, Unified_Structure_P Structure_P);
int Minimum_Index(std::vector<double> &Given_vec);
void Add_Node(Tree_Node &Cur_Node);
Tree_Node Pop_Node();
double Dot_Product(std::vector<double> &x1, std::vector<double> &x2);
void Seed_Conf_Constraint(Tree_Node &Node_i, Robot_StateNDot &x_i_child, vector<double> &sigma_i_child, vector<double> &F, vector<double> &Constraint_Status);
void Reference_Dist_Vel_Update(Tree_Node &Node_i);
