function Calibration_Stanford

% measurements
load 'P_m_stanford'
load 'Q_stanford'

addpath("../")
%% Initiaizations
T_init = eye(4,4);
n_joints = 6;
types = 'rrprrr';

Robot = RobotKinematics(n_joints, types, T_init);

% P_m = P_m(:,1:3000);
% Q = Q(:,1:3000);
m = length(P_m);


DH_real = ...
    [1 -pi/2 0 -pi/2;
    1.5 pi 0 -pi/2;
    0 0 0 0;
    0.5 -pi/2 0 pi/2;
    0 -pi/2 0 -pi/2;
    0.1 0 0 0];

% DH param limits for each link. 1 = min
% in R^njx4 [d1 theta1 a1 alpha1],;...[dn thetan an alphan]
Limits(1:n_joints,1:4,1) = ...
    [0.7 -pi/2 0 -pi/2;
    1.3 -pi 0 -pi/2;
    0 0 0 0;
    0.3 -pi/2 0 -pi/2;
    0 -pi/2 0 -pi/2;
    0.05 0 0 0];

Limits(1:n_joints,1:4,2) = ...
    [1.2 pi/2 0.1 pi/2;
    1.7 pi 0 pi/2;
    0 0 0 0;
    0.7 pi/2 0.1 pi/2;
    0 pi/2 0 pi/2;
    0.3 0 0 0];

% initial estimates
%d,theta,a,alpha
DH = ...
    [0.9 -pi/2 0 -pi/2;
    1.35 pi 0 -pi/2;
    0 0 0 0;
    0.3 -pi/2 0.05 pi/2;
    0 -pi/2 0 -pi/2;
    0.05 0 0 0];

% w_p = [0 0 1 0;
%     0 0 1 0;
%     0 0 1 0];

w_p = ones(n_joints,4);


DH_params = reshape(DH',4*n_joints,1);
DH_params_real = reshape(DH_real',4*n_joints,1);


W = 100*ones(7,m);
DH_param_lims(:,1) = reshape(Limits(:,:,1)',4*n_joints,1);
DH_param_lims(:,2) = reshape(Limits(:,:,2)',4*n_joints,1);


dim = [1,2,3,4,5,6,7]; %%motion directions
n_ref = 1;

[DH_params_pinv,W_pinv,iter_pinv] = getModel(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,n_ref,1);
[DH_params_qp,W_qp,iter_qp] = getModel(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,n_ref,2);
DH_params_compare = [DH_params_real,DH_params_pinv,DH_params_qp];

err_DH_params = [DH_params_pinv-DH_params_real,DH_params_qp-DH_params_real];

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

% %% dim = measurements dimension (generally 7)
% %% m = numb of measurements
% %% x0 = initlai DH_params estimate in R^4nj
% %% Q in R^njxm
% %% P_m = measured poses in R^dimxm [x,y,z,qx,qy,qz,qw]
% %% W = weighting matrix in R^dimxm
% %% W_p in R^4nj = weights for parameters
% %% method 1 = pinv; 2 = qp
% 
% function [DH_params,W] = getModel(Robot,dim,P_m,Q,x0,W,W_p,DH_param_lims,N_ref,method)
% 
% DH_params = x0;
% n_var = length(DH_params); %4*nj
% 
% [n_dim,m] = size(P_m);
% 
% resid = zeros(n_dim,m);
% 
% Err = 1e-08;
% dErr = 1e-10;
% Iter = 5000;
% 
% Inliers_Idx = 1:m;
% 
% % options = optimoptions('quadprog','Display','off','algorithm','trust-region-reflective');
% options = optimoptions('quadprog','Display','off');
% 
% % options_lsqlin = optimoptions('lsqlin','display','off');
% 
% 
% 
% for n_ref = 1:N_ref
%     
%     n_points = length(Inliers_Idx);
%     
%     dP = zeros(n_dim*n_points,1);
%     D = zeros(n_dim*n_points,n_var);
%     
%     iter = 1;
%     err = 1e10;
%     err_opt = err;
%     t = 0;
%     while(err > Err && iter < Iter)
%         tic
%         err_old = err;
%         
%         lb = W_p*(DH_param_lims(:,1)-x0)-1e-05;
%         ub = W_p*(DH_param_lims(:,2)-x0)+1e-05;
%         
%         iter
%         "LOADING DATA"
%         for i = 1:n_points
%             
%             j = Inliers_Idx(i);
%             
%             q = Q(:,j);
%             %             [Robot,~,P_expect] = Robot.getPose(q,DH_params);
%             
%             [Robot,~,P_expect,Dp,Dor] = Robot.getKineDeriv_Ana(q,DH_params);
%             P_expect = P_expect(dim,1);
%             
%             %             [Dp,Dor] = Robot.getDerivs(q,DH_params);
%             
%             v1 = n_dim*i-6;
%             v2 = v1+6;
%             
%             dP(v1:v2,1) = W(:,j).*(P_m(:,j)-P_expect);
%             D(v1:v2,:) = W(:,j).*[Dp;Dor]*W_p;
%         end
%         
%         lambda = 1e-03;
%         I = eye(n_var,n_var); %%needed to keep dx small
%         
%         H = D'*D+lambda*I;
%         f = -D'*dP;
%         
%                 "SOLVE QP"
%                 
%         if method == 1
% %         x = pinv(D)*dP;
%         x = -inv(H)*f;
%         
%         elseif method == 2
%        
%         %compute model
%                 [x,fval,exitflag,~] = quadprog(H,f,[],[],[],[],lb,ub,[],options);
%         
% %                  [x,fval,exitflag,] = lsqlin(D,dP,[],[],[],[],lb,ub,[],options_lsqlin);
%         
%         %         X = [x0,x,x_lsqlin]
%         
%         if exitflag == -2
%             "*******INFEASIBLE******"
%             break;
%         end
%         
%         end
%         
%         DH_params = DH_params+x
%         
%         err = 0;
%         
%         "MSE"
%         for i = 1:m
%             q = Q(:,i);
%             [Robot,~,P_expect] = Robot.getPoseNum(q,DH_params);
%             P_expect = P_expect(dim,1);
%             
%             resid(:,i) = P_m(:,i)-P_expect;
%             err_i = (norm(W(:,i).*resid(:,i)))^2;
%             err = err+err_i;
%             
%         end
%         err = err/m %mse on data with weights
%         
%         derr = abs(err-err_old)
%         err_old = err;
%         
%         iter = iter+1;
%         
%         %         t_err = toc;
%         
%         time_iter = toc
%         
%         t = t+time_iter;
%         if err < err_opt
%             err_opt = err
%             DH_params_opt = DH_params;
%         end
%         
%         if (derr < dErr || err > 10*err_opt)
%             
%             break;
%             
%         end
%         
%     end
%     t
%     
%     DH_params = DH_params_opt;
%     
%     % refine model
%     [W,Inliers_Idx] = getWeights(resid);
%     
% end
% 
% end
% 
% function [W,Inliers_Idx] = getWeights(R)
% [med,thresh] = getThreshold(R);
% 
% [n_dim,m] = size(R);
% 
% v = (R-med)./thresh;
% v = v.^2;
% 
% [~,Outliers_Idx] = find(v > 1);
% Inliers_Idx = 1:m;
% Inliers_Idx(Outliers_Idx) = [];
% 
% W = exp(-7*v); %%weights on all datapoints
% 
% end
% 
% function [m,thresh] = getThreshold(R)
% m = median(R,2);
% [n_dim,~] = size(R);
% 
% abs_dev = abs(R-m);
% 
% MAD = median(abs_dev,2);
% 
% thresh = 2*ones(n_dim,1);
% thresh = thresh.*MAD;
% 
% for i = 1:n_dim
%     if (thresh(i)) < 1e-10
%         thresh(i) = 1e10;
%     end
%     
% end
% 
% end
% 
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



