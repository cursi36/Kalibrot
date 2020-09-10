function getData_3R
clear all
close all

addpath("../")
T_init = eye(4,4);
n_joints = 3;
types = 'rrr';

% T_tool = eye(4,4);
tic
Robot = RobotKinematics(n_joints, types, T_init,[]);
toc

q = zeros(n_joints,1);
DH = [0 0 1 0;
    0 0 0.5 0;
    0 0 2 0];
DH_params = reshape(DH',4*n_joints,1);
% 
[Robot,T_val,P] = Robot.getPose(q,DH_params);

q_max = [pi/2;pi/4;pi/4];
q_min = -q_max;

C = combvec( [1,0],[1,0],[1,0]);
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
    subplot(3,2,i)
    plot(1:length(Q),Q(i,:))
end
figure()
for i = 1:7
    subplot(4,2,i)
    plot(1:length(Poses),Poses(i,:))
end

P_m = Poses;
save('P_m_3r_2','P_m')
save('Q_3r_2','Q')

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
save('P_m_3r_3','P_m')
save('Q_3r_3','Q')



end