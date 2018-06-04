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
    int Default_Init_Pr_( integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

    int Seed_Conf_Optimization_Pr_fn_(integer    *Status, integer *n,    doublereal x[],
        integer    *needF,  integer *neF,  doublereal F[],
        integer    *needG,  integer *neG,  doublereal G[],
        char       *cu,     integer *lencu,
        integer    iu[],    integer *leniu,
        doublereal ru[],    integer *lenru );

    int Real_Optimization_Pr_fn(integer    *Status, integer *n,    doublereal x[],
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
    double time;    double Kinetic_Energy;    int Node_Index_Number;    // This means the index of the current nodes in the whole node array

    Robot_StateNDot Node_StateNDot;
    Tree_Node_Ptr Par_Node;
    std::vector<Tree_Node_Ptr> Children_Nodes;

    std::vector<double> sigma_i;
    std::vector<double> StateNDot_Traj;
    std::vector<double> Ctrl_Traj;
    std::vector<double> Contact_Force_Traj;

    // The following three are reserved for inertia shaping method
    std::vector<double> self_StateNDot_Traj;
    std::vector<double> self_Ctrl_Traj;
    std::vector<double> self_Contact_Force_Traj;
    dlib::matrix<double,16,1> End_Effector_Pos, End_Effector_Vel;
};


class Unified_Structure_P{
public:
    double rIx, rIy, theta, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;
	double rIxdot, rIydot, thetadot, q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot, q10dot;
    double u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;
    Tree_Node Node_i;
    Tree_Node Node_i_child;
    std::vector<double> sigma_i;                ///        well this is the only for thye initialization
    std::vector<double> sigma_tran, sigma_goal;
    std::vector<double> Opt_Ctrl_LowBd;
    std::vector<double> Opt_Ctrl_UppBd;
    std::vector<double> Opt_Conf_LowBd;
    std::vector<double> Opt_Conf_UppBd;
    std::vector<double> Robot_State_Init;
    int Ctrl_No;    double Tme_Seed;    double eps;
    Unified_Structure_P();
};
extern Unified_Structure_P Structure_P;
extern std::vector<Tree_Node_Ptr> All_Nodes, Children_Nodes, Frontier_Nodes;

/**
 * Functions that have been successfully tested!
 * Description
 */

double Obs_Dist_Fn(std::vector<double> &r_Pos, dlib::matrix<double> &Envi_Map, int &Obs_Choice_Ind, char char_sym);
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i);
void Robot_Plot_fn(const Robot_StateNDot &StateNDot_Init_i, std::string &name);
std::vector<double> Ang_Pos_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
std::vector<double> Ang_Vel_fn(const Robot_StateNDot &Robot_StateNDot_i, const char* s);
std::vector<double> StateNDot2StateVec(Robot_StateNDot &Robot_StateNDot_i);
Robot_StateNDot StateVec2StateNDot(std::vector<double> &StateVec);
std::vector<double> Default_Init(const std::vector<double> &sigma_i, int Flag);


dlib::matrix<double> D_q_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  B_q_fn();
dlib::matrix<double>  C_q_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double>  Jac_Full_fn(const Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double> Jacdot_qdot_fn(const Robot_StateNDot &Robot_StateNDot_i);
double Kinetic_Energy_fn(Robot_StateNDot &Robot_StateNDot_i);

int Minimum_Index(std::vector<double> &Given_vec);
void Add_Node(Tree_Node &Cur_Node);
Tree_Node Pop_Node();
void Reference_Dist_Vel_Update(Tree_Node &Node_i);
double Dot_Product(std::vector<double> &x1, std::vector<double> &x2);
void End_Effector_PosNVel(Robot_StateNDot &StateNDot_Init_i, dlib::matrix<double,16,1> &End_Effector_Pos, dlib::matrix<double,16,1> &End_Effector_Vel);
dlib::matrix<double> End_Effector_Obs_Dist_Fn(dlib::matrix<double,16,1> &End_Effector_Pos);


int Node_Self_Opt(Tree_Node &Node_i);
std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child);
std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child);
std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2);
dlib::matrix<double> End_Effector_Obs_Dist_Fn(dlib::matrix<double,16,1> &End_Effector_Pos, std::vector<double> &Contact_Ind);


std::vector<double> sigma_modi(std::vector<double> sigma_ref, int contas_ind, int AddOrRet);

void Seed_Conf_Constraint(Tree_Node &Node_i, Robot_StateNDot &x_i_child, vector<double> &sigma_i_child, vector<double> &F, vector<double> &Constraint_Status);
