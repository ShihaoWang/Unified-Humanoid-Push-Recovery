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
    double time, self_time;    double Kinetic_Energy;    int Node_Index_Number;    // This means the index of the current nodes in the whole node array

    Robot_StateNDot Node_StateNDot;
    Tree_Node_Ptr Par_Node;
    std::vector<Tree_Node_Ptr> Children_Nodes;

    std::vector<double> sigma_i;
    dlib::matrix<double> StateNDot_Traj;
    dlib::matrix<double> Ctrl_Traj;
    dlib::matrix<double>Contact_Force_Traj;

    // The following three are reserved for inertia shaping method
    dlib::matrix<double> self_StateNDot_Traj;
    dlib::matrix<double> self_Ctrl_Traj;
    dlib::matrix<double> self_Contact_Force_Traj;
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
    int Opt_Prob_Flag, n_variables, n_constraints;
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


int Nodes_Optimization_fn(Tree_Node &Node_i, Tree_Node &Node_i_child, int Opt_Flag);
std::vector<double> Seed_Guess_Gene(Tree_Node &Node_i, Tree_Node &Node_i_child);
std::vector<double> Seed_Guess_Gene_Robotstate(Tree_Node &Node_i, Tree_Node &Node_i_child);
std::vector<double> Vec_Minus(std::vector<double> &vec1, std::vector<double> &vec2);
dlib::matrix<double> End_Effector_Obs_Dist_Fn(dlib::matrix<double,16,1> &End_Effector_Pos, std::vector<double> &Contact_Ind);
int Real_Optimization(std::vector<double> &Opt_Seed, std::vector<double> &sigma_tran, std::vector<double> &sigma_goal, int Opt_Flag);
void OptSeed2DlibMat(std::vector<double> &Opt_Seed, double &Tme_Seed, dlib::matrix<double> &Robot_State_Tot, dlib::matrix<double> &Ctrl_Tot, dlib::matrix<double> &Contact_Force_Tot);

std::vector<double> DlibMatCol2StdVec(dlib::matrix<double> &One_Col, int Col_Length);
std::vector<double> Opt_Constraint(double Tme_Seed, dlib::matrix<double> &Robot_State_Tot, dlib::matrix<double> &Ctrl_Tot, dlib::matrix<double> &Contact_Force_Tot, int Flag_Choice, std::vector<double> &Constraint_Type);
dlib::matrix<double> Robot_StateNDot2DlibMat(Robot_StateNDot &Robot_StateNDot_i);
dlib::matrix<double> Friction_Cone_Constraint(dlib::matrix<double> &Normal_Force, dlib::matrix<double> &Tang_Force);
void Contact_Force_Proj(dlib::matrix<double> &Contact_Force_i, std::vector<double> & Contact_Ind, dlib::matrix<double> &Normal_Forces, dlib::matrix<double> &Tang_Force);
void Constraints2ValsNType(dlib::matrix<double> &temp_result, std::vector<double> &Opt_Constraint_vals, std::vector<double> &Constraint_Type, int Flag_Choice, int Constraint_Type_val);
void StateNStatedot_Distill(dlib::matrix<double> & MatCol, dlib::matrix<double> &MatState, dlib::matrix<double> &MatStatedot);
