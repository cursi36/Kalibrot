%% dim = measurements dimension (generally 7)
%% m = numb of measurements
%% x0 = initlai DH_params estimate in R^4nj
%% Q in R^njxm
%% P_m = measured poses in R^dimxm [x,y,z,qx,qy,qz,qw]
%% W = weighting matrix in R^dimxm
%% W_p in R^4nj = weights for parameters
%% method 1 = pinv; 2 = qp

function [DH_params,W,iter] = getModel(Robot,dim,P_m,Q,x0,W,w_p,DH_param_lims,N_ref,method)

DH_params = x0;
n_var = length(DH_params); %4*nj
[n_joints,~] = size(Q);

W_p = reshape(w_p',4*n_joints,1);
f = find(W_p == 0);
W_p = diag(W_p);


W_bound = eye(n_var,n_var);
W_bound(f,f) = 0;

[n_dim,m] = size(P_m);

resid = zeros(n_dim,m);

Err = 1e-08;
dErr = 1e-10;
Err_rel = 1e-03;
Iter = 5000;

Inliers_Idx = 1:m;

% options = optimoptions('quadprog','display','off','algorithm','trust-region-reflective');
options = optimoptions('quadprog','display','off');
% options_lsqlin = optimoptions('lsqlin','display','off');



% for n_ref = 1:N_ref

n_points = length(Inliers_Idx);

dP = zeros(n_dim*n_points,1);
D = zeros(n_dim*n_points,n_var);

iter = 1;

%Initial error
for i = 1:m
    q = Q(:,i);
    [Robot,~,P_expect] = Robot.getPoseNum(q,DH_params);
    P_expect = P_expect(dim,1);
    
    resid(:,i) = P_m(:,i)-P_expect;
    err_i = (norm(W(:,i).*resid(:,i)))^2;
    err = err+err_i;
    
end
err = err/m %mse on data with weights
%     err = 1e10;
err_opt = err;
DH_params_opt = DH_params;

t = 0;
while(err > Err && iter < Iter)
    tic
    err_old = err;
    
    lb = W_bound*(DH_param_lims(:,1)-DH_params)-1e-06;
    ub = W_bound*(DH_param_lims(:,2)-DH_params)+1e-06;
    
    iter
    "LOADING DATA"
    for i = 1:n_points
        
        j = Inliers_Idx(i);
        
        q = Q(:,j);
        %             [Robot,~,P_expect] = Robot.getPose(q,DH_params);
        
        [Robot,~,P_expect,Dp,Dor] = Robot.getKineDeriv_Ana(q,DH_params);
        P_expect = P_expect(dim,1);
        
        %             [Dp,Dor] = Robot.getDerivs(q,DH_params);
        
        v1 = n_dim*i-6;
        v2 = v1+6;
        
        dP(v1:v2,1) = W(:,j).*(P_m(:,j)-P_expect);
        Der = [Dp;Dor];
        D(v1:v2,:) = W(:,j).*Der(dim,:)*W_p;
    end
    
    lambda = 1e-03;
    I = eye(n_var,n_var); %%needed to keep dx small
    
    H = D'*D+lambda*I;
    f = -D'*dP;
    
    "SOLVE QP"
    
    if method == 1
        x = -pinv(H)*f;
        
    elseif method == 2
        
        %compute model
        [x,fval,exitflag,~] = quadprog(H,f,[],[],[],[],lb,ub,[],options);
        
        %                  [x,fval,exitflag,] = lsqlin(D,dP,[],[],[],[],lb,ub,[],options_lsqlin);
        
        %         X = [x0,x,x_lsqlin]
        
        if exitflag == -2
            "*******INFEASIBLE******"
            break;
        end
        
    end
    
    DH_params = DH_params+x;
    
    err = 0;
    
    "MSE"
    for i = 1:m
        q = Q(:,i);
        [Robot,~,P_expect] = Robot.getPoseNum(q,DH_params);
        P_expect = P_expect(dim,1);
        
        resid(:,i) = P_m(:,i)-P_expect;
        err_i = (norm(W(:,i).*resid(:,i)))^2;
        err = err+err_i;
        
    end
    err = err/m %mse on data with weights
    
    derr = abs(err-err_old)
    err_old = err;
    
    err_rel = derr/err;
    
    iter = iter+1;
    
    %         t_err = toc;
    
    time_iter = toc;
    
    t = t+time_iter;
    if err < err_opt
        err_opt = err
        DH_params_opt = DH_params;
    end
    
    if (derr < dErr || err > 10*err_opt ||err_rel < Err_rel)
        
        break;
        
    end
    
end

DH_params = DH_params_opt;

%     % refine model
%     [W,Inliers_Idx] = getWeights(resid);
%
% end

end


%%%%%%%UNUSED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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