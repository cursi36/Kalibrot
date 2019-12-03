function testKinematics_3r

T_init = eye(4,4);
n_joints = 3;
types = 'rrr';

tic
Robot = RobotKinematics(n_joints, types, T_init);

toc
% T = Robot.m_T_sym;
% Dp = Robot.m_Dp_sym;
% Dor = Robot.m_Dor1_sym;
% 
q = zeros(3,1);
q(1) = 45*pi/180;
q(2) = 90*pi/180;
q(3) = -45*pi/180;
DH = [0 0 1 0; 0 0 0.5 0; 0 0 2 0]';
% DH = rand(3,4)';
DH_params = reshape(DH,4*n_joints,1);
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

q_max = [pi/2;pi/4;pi/4];
q_min = -q_max;

ns = 11; % number of samples

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
save('P_m_3r','P_m')
Q = C;
save('Q_3r','Q')



end