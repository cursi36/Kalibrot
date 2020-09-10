function testKinematics_Stanford

addpath("../")
T_init = eye(4,4);
n_joints = 6;
types = 'rrprrr';

% T_tool = eye(4,4);

tic
Robot = RobotKinematics(n_joints, types, T_init,[]);
toc

q = zeros(n_joints,1);
DH = ...
    [1 -pi/2 0 -pi/2;
    1.5 pi 0 -pi/2;
    0 0 0 0;
    0.5 -pi/2 0 pi/2;
    0 -pi/2 0 -pi/2;
    0.1 0 0 0];

DH_params = reshape(DH',4*n_joints,1);
% 
tic
[Robot,T_val,P] = Robot.getPose(q,DH_params);
[Dp,Dor] = Robot.getDerivs(q,DH_params);
toc

tic
% [Robot,T_val_num,P_num] = Robot.getPoseNum(q,DH_params);
% [Dp_num,Dor_num] = Robot.getDerivsNum(q,DH_params);
% 
% [D,quat] = Robot.getQuatDerivNum(q,DH_params);

[Robot,T_num,P_num,Dp_num,Dor_num] = Robot.getKineDeriv_Ana(q,DH_params);
toc

q_max = [pi;pi/4;0.7;pi;pi/2;pi];
q_min = [-pi;-pi/4;0;-pi;-pi/2;-pi];

ns = 5; % number of samples

Q = zeros(n_joints,ns);

for i = 1:2

    Q(i,:) = linspace(q_min(i),q_max(i),ns);

end

C = Mycombvec(Q);
[~,nc] = size(C);

for i = 1:nc

q = C(:,i);
tic
% [Robot,T_val,P(:,i)] = Robot.getPose(q,DH_params);
[Robot,T_val,P(:,i)] = Robot.getPoseNum(q,DH_params);
toc

end

P_m = P;
save('P_m_stanford','P_m')
Q = C;
save('Q_stanford','Q')



end