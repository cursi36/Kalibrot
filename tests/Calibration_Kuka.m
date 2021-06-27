function Calibration_Kuka
clear all
close all
clc

%% Load data
% measurements

%with noise
load 'P_m_Kuka_3' %pose in R^7xm = [3D position; quaternion]
load 'Q_Kuka_3' %joint values in R^njxm
m = length(P_m); %number data points

addpath("../")
%% Initiaizations
T_init = eye(4,4); %robot base transf
n_joints = 7; %number of joints
types = 'rrrrrrr'; % r = revolute; p = prismatic

Robot = RobotKinematics(n_joints, types, T_init,[]);

% real DH for simulation.
DH_real = ...
[0.34	 3.1415926	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.4	 0	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.4 	 0	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.1257	 0	 0	 0];

% DH param limits for each link. 1 = min
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
Limits(1:n_joints,1:4,1) = ...
[0.2	 pi-5*pi/180	 0	 -pi/2
0	 -pi	 0	 -pi/2
0 0	 0	 -pi/2
0	 -pi	 0	 -pi/2
0 	 0	 0	 -pi/2
0	 pi-5*pi/180	 0	 -pi/2
0	 0	 0	 0];

Limits(1:n_joints,1:4,2) = ...
[0.4	 pi	 0.1	 pi/2
0.1	 pi	 0	 pi/2
0.5	 0	 0	 pi/2
0	 pi	 0.1	 pi/2
0.5 	 0	 0	 pi/2
0	 pi	 0	 pi/2
0.3	 0	 0	 0];

% initial estimates
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
%d,theta,a,alpha
DH = ...
[0.2	 pi-2*pi/180	 0.05	 pi/2
0.05	 pi	 0	 pi/2
0.1	 0	 0	 pi/2
0	 pi	 0.1	 pi/2
0.1 	 0	 0	 pi/2
0	 pi-2*pi/180	 0	 pi/2
0.3	 0	 0	 0];

%parameters selection (0 = no optimize,1 = optimize)
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
w_p = [1 1 1 1;
    1 1 1 1;
    1 1 1 1;
    1 1 1 1;
    1 1 1 1;
    1 1 1 1
    1 1 1 1];

% motion components to consider: 1 = x,2 = y,3 = z,4 = qx,5 = qy, 6= qz, 7 = qw
dim = [1,2,3]; %%motion directions
%weights for each motion direction.
W = [1*ones(length(1),m);
    1*ones(length(1),m);
    1*ones(length(1),m);
%     1*ones(length(1),m);
%     1*ones(length(1),m);
%     1*ones(length(1),m);
%     1*ones(length(1),m)
];
P_m = P_m(dim,:); %data to consider

%% Solver selection
options.solver = "pinv";
options.damping = 1e-02; %initial damping factor for Levenberg-Marquardt algorithm
options.MaxIter = 1000; %Maximum number of iterations
options.Visualize{1} = true; %enable plotting of the results
options.Visualize{2} =[0,0,0,0,0,0,0]'; %joint values for visualization
[DH_params_pinv,P_pinv,W_pinv,Info_pinv] = Calibrate(Robot,dim,P_m,Q,DH,W,w_p,Limits,options);


options.solver = "qp";
options.damping = 1e-02;
options.MaxIter = 1000; %Maximum number of iterations
options.Visualize{1} = true;
options.Visualize{2} =[0,0,0,0,0,0,0]';
[DH_params_qp,P_qp,W_qp,Info_qp] = Calibrate(Robot,dim,P_m,Q,DH,W,w_p,Limits,options);


DH_params_real = reshape(DH_real',4*n_joints,1);
DH_params_compare = [DH_params_real,DH_params_pinv,DH_params_qp];
err_DH_params = [DH_params_pinv-DH_params_real,DH_params_qp-DH_params_real];

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





