function Calibration_3r
clear all
close all
clc

addpath("../")
%% Initiaizations
T_init = eye(4,4); %%matrix transformtaion for robot base
n_joints = 3; %%number of joints
types = 'rrr'; %type of joints (r = revolute, p = prismatic)

Robot = RobotKinematics(n_joints, types, T_init,[]); %class for robot kinematics


%% Load data
% measurements
% load 'P_m_3r'
% load 'Q_3r'

%Without noise
% load 'P_m_3r_2'
% load 'Q_3r_2'

%with noise
load 'P_m_3r_3' %joint values in R^7xm = 3D position;quaternion
load 'Q_3r_3' %joint values in R^njxm
m = length(P_m); %number of datapoints

%% Set DH bounds
% DH param limits for each link. 1 = min
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
Limits(1:n_joints,1:4,1) = ...
    [0 -5*pi/180 0.9 -5*pi/180;
    0 -5*pi/180 0.4 -5*pi/180;
    0 -5*pi/180 1.9 -5*pi/180];

Limits(1:n_joints,1:4,2) = ... 
    [0.1 5*pi/180 1.1 5*pi/180;
    0.1 5*pi/180 0.6 5*pi/180;
    0.1 5*pi/180 2.1 5*pi/180];

%Selection matrix (0 = do not optimize; 1 = optimize)
%in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
%d, theta, a alpha
w_p = [0 1 1 0;
    0 1 1 0;
    0 1 1 0];

% w_p = ones(n_joints,4);

% real DH parameters used for the simulated robot. here they are used as
% ground truth
DH_real = [0 0 1 0;
    0 0 0.5 0;
    0 0 2 0];

% initial DH estimates
%in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
%d,theta,a,alpha
DH = [0 0 0.9 0;
    0 0 0.4 0;
    0 0 1.9 0];

% DH = 1e-03*ones(3,4);

% motion components to consider: 1 = x,2 = y,3 = z,4 = qx,5 = qy, 6= qz, 7 = qw
dim = [1,2,6,7]; %%motion directions

%weights for each motion direction.
W = 1*ones(length(dim),m);

%data to consider
P_m = P_m(dim,:);


%% Solver selection
options.solver = "pinv";
options.damping = 1e-0; %initial damping factor for Levenberg-Marquardt algorithm
options.Visualize{1} = true; %enable plotting of the results
options.Visualize{2} = [0,0,0]; %joint values for visualization
[DH_params_pinv,P_pinv,W_pinv,Info_pinv] = Calibrate(Robot,dim,P_m,Q,DH,W,w_p,Limits,options);

options.solver = "qp";
options.damping = 1e-0;
options.Visualize{1} = true;
options.Visualize{2} = [];
[DH_params_qp,P_qp,W_qp,Info_qp] = Calibrate(Robot,dim,P_m,Q,DH,W,w_p,Limits,options);

DH_params_real = reshape(DH_real',4*n_joints,1);
DH_params_compare = [DH_params_real,DH_params_pinv,DH_params_qp];

err_DH_params = [DH_params_pinv-DH_params_real,DH_params_qp-DH_params_real];

% save('Info_qp_3r','Info_qp')
% save('Info_pinv_3r','Info_pinv')

end


function P = getEstimate(Robot,Q,DH_params)
[~,n_points] = size(Q);

P = zeros(7,n_points);

for i = 1:n_points
    
    q = Q(:,i);
    [Robot,~,P(:,i)] = Robot.getPoseNum(q,DH_params);
%     [Robot,~,P(:,i)] = Robot.getPose(q,DH_params);
    
end

end
