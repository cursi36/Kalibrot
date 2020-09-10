# Kalibrot
Algorithm for Robot Kinematic Calibration

- Kalibrot is an optimization algorithm for solving the problem of finding the optimal DH parametres for correct robot kinematc calibration.
- The algorithm uses derivatives of the Cartesian position and orientation (computed through quaternions) which are retrieved analytically, thus speeding up the computations.
- Two different methods can be used: 
    1) traditional **pseudoinverse**;
    2) **constarined quadratic programming** problem (solved using quadprog from matlab).

## Functions
The main functions are:
- `getModel.m` which solves the optimization problem. It finds the optimal DH parameters and also outputs additional information such as the identifiable parameters or the observability measures;
- `RobotKinematics.m` is a class building the Robot object which allows to compute the forward kinematics and the derivatives wrt DH parameters. It just takes as input the number of joints, their types, and the initial transfromation matrix.
- `VisualizeResults.m` generates three plots for showing the calibrate DH parameters, the calibrated robot kinematic structure, and the identification matrix;
- The others are just auxiliary functions.

### Tests
The `tests` folder contains two examples:
- calibration  of a 3R robot;
- calibration of a Stanford manipulator (6DOF);
- calibration of a KUKA manipulator (7DOF).

- `getData_3R.m`, `getData_Stanford.m`, `getData_Kuka.m` are used to generate the calibration data for the three simulated robots.

### Tutorial
`tutorial.pdf` explains Kalibrot functionalities and how to use it.

### Contact
If you have any question and want to contribute to the project, please email `cursifrancesco@gmail.com` .
We thank you for your interest and help.

