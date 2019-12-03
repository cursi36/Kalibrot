# Kalibrot
Algorithm for Robot Kinematic Calibration

- Kalibrot is an optimization algorithm for solving the problem of finding the optimal DH parametres for correct robot kinematc calibration.
- The algorithm uses derivatives of the Cartesian position and orientation (computed through quaternions) which are retrieved analytically, thus speeding up the computations.
- Two different methods can be used: 
    1) traditional **pseudoinverse**
    2) **constarined quadratic programming** problem (solved using quadprog from matlab).

## Functions

The main functions are:
- `getModel.m` which solves the optimization problem
- `RobotKinematics.m` is a class building the Robot object which allows to compute the forward kinematics and the derivatives wrt DH parameters. It just takes as input the number of joints, their types, and the initial transfromation matrix.
- `Mycombvec.m` and `rsqrd.m` are just auxiliary functions.


### Tests
The `tests` folder contains two examples:
- calibration  of a 3R robot
- calibration of a Stanford manipulator (6DOF)

