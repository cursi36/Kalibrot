function Calibration_3r


%% Initiaizations
T_init = eye(4,4);
n_joints = 3;
types = 'rrr';

Robot = RobotKinematics(n_joints, types, T_init);

% T = Robot.m_T_sym;
% Dp = Robot.m_Dp_sym;

% measurements
load 'P_m_3r'
load 'Q_3r'

% P_m = P_m(:,1:100);
% Q = Q(:,1:100);
m = length(P_m);

% DH param limits for each link. 1 = min
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
Limits(1:n_joints,1:4,1) = ...
    [0 -pi 0.5 -pi;
    0 -pi 0.2 -pi;
    0 -pi 1.9 -pi];

Limits(1:n_joints,1:4,2) = ... 
    [1 pi 1.5 pi;
    1 pi 0.6 pi;
    1 pi 2.1 pi];

w_p = [0 0 1 0;
    0 0 1 0;
    0 0 1 0];

% w_p = ones(n_joints,4);

DH_real = [0 0 1 0;
    0 0 0.5 0;
    0 0 2 0];

% initial estimates
%d,theta,a,alpha
% DH = [0 0 0.5 0;
%     0 0 0.3 0;
%     0 0 1.9 0];

DH = 1e-03*ones(3,4);


DH_params = reshape(DH',4*n_joints,1);
DH_params_real = reshape(DH_real',4*n_joints,1);

W = ones(7,m);
DH_param_lims(:,1) = reshape(Limits(:,:,1)',4*n_joints,1);
DH_param_lims(:,2) = reshape(Limits(:,:,2)',4*n_joints,1);


dim = [1,2,3,4,5,6,7]; %%motion directions
n_ref = 1;

[DH_params_pinv,W_pinv,iter_pinv] = getModel(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,n_ref,1);
[DH_params_qp,W_qp,iter_qp] = getModel(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,n_ref,2);
DH_params_compare = [DH_params_real,DH_params_pinv,DH_params_qp];

err_DH_params = [DH_params_pinv-DH_params_real;DH_params_qp-DH_params_real];

P = getEstimate(Robot,Q,DH_params);
R_2 = rsqrd(P_m,P);



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
% %%%%%%%%%%%%%%%%%%%%%%%
% % iterate also on weights
% function [DH_params,W] = getModel_iterWeight(Robot,dim,P_m,Q,x0,W,W_p,DH_param_lims)
% 
% DH_params = x0;
% n_var = length(DH_params); %4*nj
% 
% [n_dim,m] = size(P_m);
% 
% resid = zeros(n_dim,m);
% 
% err = 1;
% Err = 1e-06;
% dErr = 1e-06;
% iter = 1;
% Iter = 5000;
% 
% 
% Inliers_Idx = 1:m;
% 
% options = optimoptions('quadprog','display','notify');
% 
% while(err > Err && iter < Iter)
%     
%     n_points = length(Inliers_Idx);
%     
%     dP = zeros(n_dim*n_points,1);
%     D = zeros(n_dim*n_points,n_var);
%     
%     err_old = err;
%     
%     lb = DH_param_lims(:,1)-x0;
%     ub = DH_param_lims(:,2)-x0;
%     
%     for i = 1:n_points
%         
%         j = Inliers_Idx(i);
%         
%         q = Q(:,j);
%         [Robot,~,P_expect] = Robot.getPose(q,DH_params);
%         P_expect = P_expect(dim,1);
%         [Dp,Dor] = Robot.getDerivs(q,DH_params);
%         
%         v1 = n_dim*i-6;
%         v2 = v1+6;
%         
%         dP(v1:v2,1) = W(:,j).*(P_m(:,j)-P_expect);
%         D(v1:v2,:) = W(:,j).*[Dp;Dor];
%     end
%     
%     H = D'*D;
%     f = -D'*dP;
%     %compute model
%     [x,fval,exitflag,~] = quadprog(H,f,[],[],[],[],lb,ub,DH_params,options);
%     
%     if exitflag == -2
%         "*******INFEASIBLE******"
%         break;
%     end
%     
%     DH_params = DH_params+x;
%     
%     err = 0;
%     
%     for i = 1:m
%         q = Q(:,i);
%         [Robot,~,P_expect] = Robot.getPose(q,DH_params);
%         P_expect = P_expect(dim,1);
%         
%         resid(:,i) = P_m(:,i)-P_expect;
%         err_i = (norm(W(:,i).*resid(:,i)))^2;
%         err = err+err_i;
%         
%     end
%     err = err/m; %mse
%     
%     [W,Inliers_Idx] = getWeights(resid);
%     
%     derr = abs(err-err_old);
%     iter = iter+1;
%     if (derr < dErr)
%         
%         break;
%         
%     end
%     
% end
% 
% 
% end


