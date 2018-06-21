#! /usr/bin/env python

import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
import pickle as pkl
from klampt import WorldModel
from klampt import vis
from klampt.math import vectorops
from klampt.model.trajectory import Trajectory

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

Grids = 0
pi = 3.1415926535897932384626
Aux_Link_Ind = [1, 3, 5, 6, 7, 11, 12, 13, 17, 18, 19, 20, 21, 23, 24, 26, 27, 28, 30, 31, 33, 34, 35]  # This is the auxilliary joints causing lateral motion
Act_Link_Ind = [0, 2, 4, 8, 9, 10, 14, 15, 16, 22, 25, 29, 32]                                          # This is the active joints to be considered
Ctrl_Link_Ind = Act_Link_Ind[3:]
Tot_Link_No = len(Aux_Link_Ind) + len(Act_Link_Ind)

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.nextone = False

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed",c)
        if c == 'q':
            self.quit = True
            return True
        if c == 'c':
            self.nextone = True
            return True
        return False

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]


def getPath(data):
    """Load configurations and return a traj, i.e. list of list

    :param data: a dictionary structure from parsing the solution. It has 't', 'x', 'u', etc
    :return t, traj, force: time stamp, configuration, and support force

    """
    t = data['t']
    x = data['x'][:, :7]  # only use q
    force = data['p']
    x_offset = 0
    # parse into traj
    traj = []
    for i in range(len(t)):
        qT, q1, q2, p1, p2, qx, qy = x[i]
        q = np.zeros(10)
        q[0] = -qx  # x-loc
        q[2] = qy  # z-loc
        q[4] = qT
        q[6] = q1
        q[7] = q2
        q[8] = p1
        q[9] = p2
        traj.append(q.tolist())
    return t, traj, force


def copy_copy(q, switched):
    qq = copy.copy(q)
    if switched:
        qq[6], qq[8] = qq[8], qq[6]
        qq[7], qq[9] = qq[9], qq[7]
    return qq

def Dimension_Recovery(low_dim_obj):
    high_dim_obj = np.zeros(Tot_Link_No)
    for i in range(0,len(Act_Link_Ind)):
        high_dim_obj[Act_Link_Ind[i]] = low_dim_obj[i]
    return high_dim_obj

def main():
    #creates a world and loads all the items on the command line
    world = WorldModel()
    res = world.readFile("/home/shihao/Multi-contat Push Recovery-Python/HRP2_Robot.xml")
    if not res:
        raise RuntimeError("Unable to load model")
    robot = world.robot(0)
    # coordinates.setWorldModel(world)
    plugin = MyGLPlugin(world)
    vis.pushPlugin(plugin)
    #add the world to the visualizer
    # vis.add("world", world)
    # vis.add("robot",world.robot(0))
    # vis.show()

    Kinetic_Energy = []  # This list is used to save the value of the kinetic energy

    T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj = Path_Loader()

    Time = np.linspace(0, T_tot, num=StateNDot_Traj.shape[1])

    for i in range(0, StateNDot_Traj.shape[1]):
        StateNDot_Traj_i = StateNDot_Traj[:,i]
        KE_i = KE_fn(robot, StateNDot_Traj_i)
        Kinetic_Energy.append(KE_i)

    fig, ax = plt.subplots()
    ax.plot(Time, Kinetic_Energy)
    ax.set(xlabel='time (s)', ylabel='voltage (mV)',title='About as simple as it gets, folks')
    ax.grid()

    fig.savefig("test.png")
    plt.show()


    # while vis.shown():
    #     pass

    # parse arguments and decide what we can do
    # args = getOnOffArgs('one', 'two')
    # if args.one:
    #     showOneStep(plugin, world, robot, ground1, ground2, ground3, link_foot, link_foot_other)
    # if args.two:
    #     showTwoStep(plugin, world, robot, ground1, ground2, ground3, ground4, link_foot, link_foot_other)
def Robot_ConfigNVel_Update(robot, x):
    OptConfig_Low = x[0:len(x)/2]
    OptVelocity_Low = x[len(x)/2:]
    OptConfig_High = Dimension_Recovery(OptConfig_Low)
    OptVelocity_High = Dimension_Recovery(OptVelocity_Low)
    robot.setConfig(OptConfig_High)
    robot.setVelocity(OptVelocity_High)

def KE_fn(robot, dataArray):
    #First we have to set the robot to be the corresponding configuration and angular velocities
    Robot_ConfigNVel_Update(robot, dataArray)
    D_q = np.asarray(robot.getMassMatrix())    # Here D_q is the effective inertia matrix
    qdot_i = dataArray[13:None]
    qdot_i = np.reshape(qdot_i,[13,1])
    qdot_i = Dimension_Recovery(qdot_i)
    qdot_i_trans = np.transpose(qdot_i)
    T = 0.5 * qdot_i_trans.dot(D_q.dot(qdot_i))
    return T

def Path_Loader():
    # This function is used to read the Opt_Soln txt file
    global Grids
    with open("Opt_Soln3.txt",'r') as robot_soln_file:
        robot_soln_i = robot_soln_file.readlines()
        robot_soln = [x.strip() for x in robot_soln_i]
        robot_soln = np.array(robot_soln, dtype = float)
    Grids = (robot_soln.shape[0]-1)/48
    T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj = Seed_Guess_Unzip(robot_soln)
    return T_tot, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj

def Seed_Guess_Unzip(Opt_Seed):
    # This function is used to unzip the Opt_Seed back into the coefficients
    T = Opt_Seed[0]
    # First is to retrieve the stateNdot coefficients from the Opt_Seed: 1 to 1 + (Grids - 1) * 26 * 4
    Var_Init = 1
    Var_Length = (Grids) * 26
    Var_End = Var_Init + Var_Length
    StateNDot_Traj = Opt_Seed[Var_Init:Var_End]
    StateNDot_Traj = np.reshape(StateNDot_Traj, (Grids, 26)).transpose()

    Var_Init = Var_End
    Var_Length = (Grids) * 10
    Var_End = Var_Init + Var_Length
    Ctrl_Traj = Opt_Seed[Var_Init:Var_End]
    Ctrl_Traj = np.reshape(Ctrl_Traj, (Grids, 10)).transpose()

    Var_Init = Var_End
    Var_Length = (Grids) * 12
    Var_End = Var_Init + Var_Length
    Contact_Force_Traj = Opt_Seed[Var_Init:Var_End]
    Contact_Force_Traj = np.reshape(Contact_Force_Traj, (Grids, 12)).transpose()
    # ipdb.set_trace()
    return T, StateNDot_Traj, Ctrl_Traj, Contact_Force_Traj

def showTwoStep(plugin, world, robot, ground1, ground2, ground3, ground4, link_foot, link_foot_other):
    """Read the data file and select trajectories to show, they are composed of two steps"""
    data = np.load('data/twoStepBunch.npz')['output']
    ndata = len(data)
    N = 20
    force_len = 0.5
    vis.show()
    for j in range(ndata):
        i = j
        print('Entering traj %d' % i)
        l0, l1, h0, h1, vel, phase0, phase1 = data[i]
        # set those grounds
        ground1.setConfig([0, 0, 0, 0, 0, 0])
        ground2.setConfig([l0, 0, h0, 0, 0, 0])
        ground3.setConfig([l0 + l1, 0, h0 + h1, 0, 0, 0])
        ground4.setConfig([2 * l0 + l1, 0, 2 * h0 + h1, 0, 0, 0])
        while vis.shown() and not plugin.quit:
            if plugin.nextone:  # check if we want next one
                plugin.nextone = False
                break
            # show phase0
            t, q, force = getTraj(phase0)
            h_ = t[1] - t[0]
            if h_ < 0:
                break
            nPoint = len(t)
            for k in range(nPoint):
                vis.lock()
                useq = copy_copy(q[k], False)
                robot.setConfig(useq)
                footpos = link_foot.getWorldPosition([0, 0, 0.5])
                support = np.array([force[k, 0], 0, force[k, 1]])
                use_support = support / np.linalg.norm(support) * force_len
                force_end = vectorops.add(footpos, use_support.tolist())
                vis.add("force", Trajectory([0, 1], [footpos, force_end]))
                vis.unlock()
                time.sleep(h_ * 5.0)
                vis.remove('force')
            # phase 1
            t, q, force = getTraj(phase1)
            h_ = t[1] - t[0]
            if h_ < 0:
                break
            print('h_ = %f' % h_)
            for k in range(nPoint):
                vis.lock()
                useq = copy_copy(q[k], True)
                robot.setConfig(useq)
                footpos = link_foot_other.getWorldPosition([0, 0, 0.5])
                support = np.array([force[k, 0], 0, force[k, 1]])
                use_support = support / np.linalg.norm(support) * force_len
                force_end = vectorops.add(footpos, use_support.tolist())
                vis.add("force", Trajectory([0, 1], [footpos, force_end]))
                vis.unlock()
                time.sleep(h_ * 5.0)
                vis.remove('force')
    while vis.shown() and not plugin.quit:
        vis.lock()
        vis.unlock()
        time.sleep(0.05)
    print("Ending visualization.")
    vis.kill()

def showOneStep(plugin, world, robot, ground1, ground2, ground3, link_foot, link_foot_other):
    """Read the data file and select a few trajectories to show"""
    # calculate the trajectory to show
    oneStepDBName = 'data/oneStepDatabase.npz'
    data = np.load(oneStepDBName)['arr_0']
    ndata = len(data)
    force_len = 0.5
    repeatN = 3
    vis.show()
    for i in range(ndata):
        print('Entering traj %d' % i)
        nowN = 0
        datai = data[i]
        h, l, sol = datai['h'], datai['l'], datai['sol']
        t, q, force = getTraj(sol)  # parse into useable format
        nPoint = len(t)
        h_ = t[1] - t[0]  # actual time step for trajectory
        # vis.animate(("world","marlo"), traj, speed=10)
        while vis.shown() and not plugin.quit:
            if nowN == 0:
                switched = False
                if plugin.nextone:  # try to exit once animation of one mode is finished
                    plugin.nextone = False
                    break
            for k in range(nPoint):
                vis.lock()
                if k == 0:  # manipulate the two steps if necessary
                    offset = nowN * np.array([l, 0, h])
                    ground1.setConfig([-(0 + offset[0]), 0 + offset[1], 0 + offset[2], 0, 0, 0])
                    ground2.setConfig([-(l + offset[0]), 0 + offset[1], h + offset[2], 0, 0, 0])
                    ground3.setConfig([-(2*l + offset[0]), 0 + offset[1], 2*h + offset[2], 0, 0, 0])
                    nowN += 1
                    if nowN == repeatN:
                        nowN = 0
                useq = copy_copy(q[k], switched)
                useq[0] -= offset[0]
                useq[2] += offset[2]
                robot.setConfig(useq)
                footpos = link_foot.getWorldPosition([0, 0, 0.4])
                otherfootpos = link_foot_other.getWorldPosition([0, 0, 0.4])
                if not switched:
                    footpos, otherfootpos = otherfootpos, footpos
                if k == 0:
                    print('offset {}'.format(offset))
                    print('useq {}'.format(useq))
                    print('foot1 {} foot2 {}'.format(footpos, otherfootpos))
                support = np.array([force[k, 0], 0, force[k, 1]])
                use_support = support / np.linalg.norm(support) * force_len
                force_end = vectorops.add(footpos, use_support.tolist())
                vis.add("force", Trajectory([0, 1], [footpos, force_end]))
                vis.unlock()
                #changes to the visualization must be done outside the lock
                time.sleep(h_ * 2.0)  # wait for actual time
                vis.remove('force')
            switched = not switched
        print('show next one')
    while vis.shown() and not plugin.quit:
        vis.lock()
        vis.unlock()
        time.sleep(0.05)
    print("Ending visualization.")
    vis.kill()


if __name__ == '__main__':
    main()
