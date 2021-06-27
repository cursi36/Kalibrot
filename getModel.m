%% Iputs:
%% dim = measurements dimension (generally 7)
%% m = numb of measurements
%% x0 = initlai DH_params estimate in R^4nj
%% Q in R^njxm
%% P_m = measured poses in R^dimxm [x,y,z,qx,qy,qz,qw]
%% W = weighting matrix in R^dimxm
%% W_p in R^4nj = weights for parameters
%% method 1 = pinv; 2 = qp
%% returns:
%% DH_params = parameters_optimized
%% P = Estimated poses
%% W = weighting matrix
%% Info = optimization info:
%% iter =  number of iterations
%% mse = Measn squared pose error between measurements and estimates for each component
%% R2 = r squared for goodness of fit
%% err_iter = minimum cost function error for each iteration
%% dx = DH parameters change for each iteration;
%% dx_base = base DH parameters change for each iteration;
%% base_params = index of independent parameters for each iter;
%% base_params_values = base DH parameters values for each iter;
%% PermMat = identification matrix K at each iteration;
%% observability = 3 observability measures for each iteration.

%     n_pars = length(x(comp_base,:)); %paramters estimated
%     [obs,~] = AnalyzeData(dP,D(:,comp_base),x(comp_base,:),n_pars,m);
%     Info.Observability(:,iter) = obs;

function [DH_params,P,W,Info] = getModel(Robot,dim,P_m,Q,x0,W,w_p,DH_param_lims,options)

damping = options.damping;
solver = options.solver;

DH_params = x0;
n_var = length(DH_params); %4*nj
[n_joints,~] = size(Q);

W_p = reshape(w_p',4*n_joints,1);
f = find(W_p == 0);
f_optimize = find(W_p == 1);
W_p = diag(W_p);
f_notOptimize = f;

% W_bound = eye(n_var,n_var);
% W_bound(f,f) = [];

[n_dim,m] = size(P_m);

resid = zeros(n_dim,m);

Err = 1e-08;
dErr = 1e-10;
Err_rel = 1e-03;
Iter = options.MaxIter;

Inliers_Idx = 1:m;

% options = optimoptions('quadprog','display','off','algorithm','trust-region-reflective');
options = optimoptions('quadprog','display','off');
% options_lsqlin = optimoptions('lsqlin','display','off');



% for n_ref = 1:N_ref

n_points = length(Inliers_Idx);

iter = 1;

%Initial error
err = 0;
for i = 1:m
    q = Q(:,i);
    [Robot,~,P_expect] = Robot.getPoseNum(q,DH_params);
    P_expect = P_expect(dim,1);
    
    resid(:,i) = P_m(:,i)-P_expect;
    err_i = (norm(W(:,i).*resid(:,i)))^2;
    err = err+err_i;
    
end
err = err/m; %mse on data with weights

disp("Initial err = "+num2str(err));
%     err = 1e10;
err_opt = err;
DH_params_opt = DH_params;

t = 0;

lambda = damping;
err_old = err;
err_prev = err_old;

while(err > Err && iter < Iter)
    tic
    
    dP = zeros(n_dim*n_points,1);
    D = zeros(n_dim*n_points,n_var);
    x_tot = zeros(n_var,1);
    x_base_tot = zeros(n_var,1);
    
    DH_params(2:4:end,1) = atan2(sin(DH_params(2:4:end,1)),cos(DH_params(2:4:end,1)));
    DH_params(4:4:end,1) = atan2(sin(DH_params(4:4:end,1)),cos(DH_params(4:4:end,1)));
    
    lb = (DH_param_lims(:,1)-DH_params)-1e-06;
    ub = (DH_param_lims(:,2)-DH_params)+1e-06;
    
    lb(f_notOptimize) = [];
    ub(f_notOptimize) = [];
    
    %     iter
    %     "LOADING DATA"
    for i = 1:n_points
        
        j = Inliers_Idx(i);
        
        q = Q(:,j);
        %             [Robot,~,P_expect] = Robot.getPose(q,DH_params);
        
        [Robot,~,P_expect,Dp,Dor] = Robot.getKineDeriv_Ana(q,DH_params);
        P_expect = P_expect(dim,1);
        
        %             [Dp,Dor] = Robot.getDerivs(q,DH_params);
        
        v1 = n_dim*i-n_dim+1;
        v2 = v1+n_dim-1;
        
        dP(v1:v2,1) = W(:,j).*(P_m(:,j)-P_expect);
        Der = [Dp;Dor];
        D(v1:v2,:) = W(:,j).*Der(dim,:)*W_p;
    end
    
    D(:,f_notOptimize) = [];
    
    H = D'*D;
    [U,S,V] = svd(H);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%reduce singularities
    %     eig = diag(S);
    %     flag_eig = find(eig/max(eig) < 1e-04);
    %     S(flag_eig,flag_eig) = 0;
    %     H = U*S*V';
    %     [~,S,V] = svd(H);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ran = rank(H);
    
    if ran < length(f_optimize)
        
        [V21,V22,comb] = getPermutationMat(H,ran); %%comb = parameters that are in combination with others
        
        comp = 1:length(f_optimize);
        comp_base = comp;
        comp_base(comb) = []; %% parameters that can be optimized independetly
        
    else
        
        comp = 1:length(f_optimize);
        comp_base = comp;
    end
    
    I = eye(length(f_optimize),length(f_optimize)); %%needed to keep dx small
    
    H = D'*D+lambda*I; %add penalization matrix for each parameter
    f = -D'*dP;
    
    disp("*****SOLVE OPTIMIZATION")
    disp("iter = "+num2str(iter))
    
    if solver == "pinv"
        x = -pinv(H)*f;
        
    elseif solver == "qp"
        
        x0_qp = initializeSol(H,f,lb,ub);
        %compute model
        [x,fval,exitflag,~] = quadprog(H,f,[],[],[],[],lb,ub,x0_qp,options);
        
        %                  [x,fval,exitflag,] = lsqlin(D,dP,[],[],[],[],lb,ub,[],options_lsqlin);
        
        %         X = [x0,x,x_lsqlin]
        
        if exitflag == -2
            "*******INFEASIBLE******"
            break;
        end
    elseif solver == "gd"
        
        
    end
    
    
    if ran < length(f_optimize)
        Perm_comp_base = zeros(length(comp_base),length(x));
        for n_cb = 1:length(comp_base)
            Perm_comp_base(n_cb,comp_base(n_cb)) = 1;
        end
        Perm_comb = zeros(length(comb),length(x));
        
        for n_c = 1:length(comb)
            Perm_comb(n_c,comb(n_c)) = 1;
        end
        
        x_base = [x(comp_base,:)-V21*inv(V22)*x(comb,:);zeros(length(comb),1)];
        
        K = [Perm_comp_base-V21*inv(V22)*Perm_comb;zeros(length(comb),length(x))];
    else
        x_base = x;
        K = eye(length(f_optimize));
        
    end
    K_tot = zeros(n_var,n_var);
    K_tot(f_optimize,f_optimize) = K;
    
    
    x_tot(f_optimize,1) = x;
    x_base_tot(f_optimize,1) = x_base;
    DH_params_old = DH_params;
    %     DH_params_base = [DH_params(comp_base,:);DH_params(comb,:)]+x_base;
    DH_params_base = DH_params+x_base_tot;
    DH_params = DH_params+x_tot;
    
    err = 0;
    
    %     "MSE"
    for i = 1:m
        q = Q(:,i);
        [Robot,~,P_expect] = Robot.getPoseNum(q,DH_params);
        P_expect = P_expect(dim,1);
        
        resid(:,i) = P_m(:,i)-P_expect;
        err_i = (norm(W(:,i).*resid(:,i)))^2;
        err = err+err_i;
        
    end
    err = err/m; %mse on data with weights
    
    
    disp("err = "+num2str(err));
    disp("lambda = "+num2str(lambda));
    
    if solver ~= "gd"
        
        if err > err_old
            
            lambda = lambda*10;
            
            lambda = min(lambda, 1e03);
            
            DH_params = DH_params_old;
            
            derr = abs(err-err_prev);
            err_prev = err;
            
            % STops solvers if function not decreasing anymore
            if (lambda >= 1e06 && derr < dErr )
                disp("ERROR NOT DECREASING. Error change: "+num2str(derr))
                break;
            end
            
        elseif err <= err_old
            
            lambda = lambda/5;
            lambda = max(lambda, 1e-06);
            
            derr = abs(err-err_old);
            err_old = err;
            err_rel = derr/err;
            
            if lambda <= 1e-06
                % STops solvers if function not decreasing anymore
                %                 if (derr < dErr ||err_rel < Err_rel)
                if (derr < dErr )
                    disp("ERROR NOT DECREASING. Error change: "+num2str(derr)+" < "+num2str(dErr))
                    break;
                    
                end
                
            end
            
        end
        
    end
    
    if err < err_opt
        err_opt =  err;
        disp("err opt = "+num2str(err_opt));
        DH_params_opt = DH_params;
        x_opt = x;
        dP_opt= dP;
        D_opt= D;
    end
    
    Info.err_iter(iter) = err_opt;
    
    Info.dx(:,iter) = x_tot;
    Info.dx_base(:,iter) = x_base_tot;
    Info.base_params{iter} = f_optimize(comp_base);
    Info.base_params_values(:,iter) = DH_params_base;
    Info.PermMat(:,:,iter) = K_tot;
    
    n_pars = length(x(comp_base,:)); %paramters estimated
    [obs,~] = AnalyzeData(dP,D(:,comp_base),x(comp_base,:),n_pars,m);
    Info.Observability(:,iter) = obs;
    %     Info.Confidence = conf;
    
    iter = iter+1;
    
    time_iter = toc;
    
    t = t+time_iter;
    
    if norm(x) <= 1e-08
        disp("PARAMETERS NOT CHANGING. Paramteres change: "+num2str(norm(x))+" < "+num2str(1e-08))
        break;
    end
    
end
disp("MAX ITER REACHED")
DH_params = DH_params_opt;


P = getEstimate(Robot,Q,DH_params);
R2 = rsqrd(P_m,P(dim,:));

err = zeros(length(dim),1);
for i = 1:length(dim)
    err(i) = mse(P_m(i,:)-P(dim(i),:));
    
end

Info.iter = iter-1;
Info.R2 = R2;
Info.mse = err;


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


function [V21,V22,comb] = getPermutationMat(H,ran)

[U,S,V] = svd(H);
V2 = V(:,ran+1:end);

n = size(V2,1); %%number of variables
dim = n-ran;

% Comb = nchoosek(n-1:-1:1,dim-1);

Comb = nchoosek(n:-1:1,dim);

for i = 1:size(Comb,1)
    comb = Comb(i,:);
    
    %     V22 = [V2(comb,:);V2(end,:)];
    V22 = [V2(comb,:)];
    
    H_base = H;
    %     comb = [comb,n];
    H_base(:,comb) = [];
    
    det_V22 = det(V22);
    if(abs(det_V22) > 1e-06)
        ran_V22 = rank(V22);
    else
        ran_V22 = rank(V22)-1;
    end
    
    
    if ran_V22 == dim && rank(H_base) == ran
        
        break;
        
    end
end

V21 = V2;
V21(comb,:) = [];


end

function x0_qp = initializeSol(H,f,lb,ub)

x0_qp = -pinv(H)*f;

for i = 1:length(x0_qp)
    if (x0_qp(i) > lb(i) && x0_qp(i) < ub(i))
        
    else
        
        x0_qp(i) = (ub(i)+lb(i))/2;
        
    end
    
end

end