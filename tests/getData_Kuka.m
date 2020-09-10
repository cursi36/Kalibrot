function getData_Kuka
clear all
close all

addpath("../")
T_init = eye(4,4);
n_joints = 7;
types = 'rrrrrrr';

% T_tool = eye(4,4);

tic
Robot = RobotKinematics(n_joints, types, T_init,[]);
toc

q = zeros(n_joints,1);
DH = ...
    [0.34	 3.1415926	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.4	 0	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.4 	 0	 0	 1.5707963
0	 3.1415926	 0	 1.5707963
0.1257	 0	 0	 0];

DH_params = reshape(DH',4*n_joints,1);
% 
tic
[Robot,T_val,P] = Robot.getPose(q,DH_params);

q_max = [pi;pi/2;pi/2;pi/2;pi/2;pi/2;pi/2];
q_min = [-pi;-pi/2;-pi/2;-pi/2;-pi/2;-pi/2;pi/2];

C(1:n_joints,1) = 1;
C(1:n_joints,2) = 0;
C = Mycombvec( C);
C(:,end) = [];

[~,nc] = size(C);
N_pnts = 50;

Q = [];
Poses = [];
for i = 1:nc
    n_pnts = 0;
    comb = C(:,i);
    while (n_pnts <= N_pnts)
      
        wave = (q_max-q_min)/2*sin(n_pnts/N_pnts*2*pi)+(q_max+q_min)/2;
        
        q = wave.*comb;
tic
[Robot,T_val,P] = Robot.getPoseNum(q,DH_params);
toc

n_pnts = n_pnts+1;

Q = [Q,q];
Poses = [Poses,P];
         
    end

end

figure()
for i = 1:n_joints
    subplot(4,2,i)
    plot(1:length(Q),Q(i,:))
end
figure()
for i = 1:7
    subplot(4,2,i)
    plot(1:length(Poses),Poses(i,:))
end

P_m = Poses;
save('P_m_Kuka','P_m')
save('Q_Kuka','Q')

%%%Add noise
noise = [10*1e-02*ones(3,1);5*1e-02*ones(4,1)].*randn(7,length(Poses));
P_m = Poses+noise;
figure()
for i = 1:7
    subplot(4,2,i)
    hold on
    plot(1:length(Poses),Poses(i,:),'-k')
    plot(1:length(P_m),P_m(i,:),'.','Color',[0.7 0.7 0.7])
end
save('P_m_Kuka_3','P_m')
save('Q_Kuka_3','Q')



end