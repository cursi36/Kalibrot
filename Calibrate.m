function [DH_params,P,W,Info] = Calibrate(Robot,dim,P_m,Q,DH,W,w_p,Limits,options)
clc
disp("*********CALIBRATING********")

[n_joints,~] = size(Q);
DH_params = reshape(DH',4*n_joints,1);

DH_params_init = DH_params;

DH_param_lims(:,1) = reshape(Limits(:,:,1)',4*n_joints,1);
DH_param_lims(:,2) = reshape(Limits(:,:,2)',4*n_joints,1);

if isempty(options.damping) == 1
    
    options.damping = 1e-03;
end

if options.solver == "qp"
     f = find(DH_params < DH_param_lims(:,1));
     if (isempty(f) == 0)
         disp("VALUES OUT OF BOUNDS "+num2str(f))
         disp(num2str(DH_params(f)')+" < "+num2str(DH_param_lims(f,1)'))
     end
     f = find(DH_params > DH_param_lims(:,2));
     if (isempty(f) == 0)
          disp("VALUES OUT OF BOUNDS "+num2str(f))
         disp(num2str(DH_params(f)')+" > "+num2str(DH_param_lims(f,1))')
     end
    
end

% if options.OptimizeBaseParams == false
[DH_params,P,W,Info] = getModel(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,options);
% else
%    [DH_params,P,W,Info] = getModelBaseParams(Robot,dim,P_m,Q,DH_params,W,w_p,DH_param_lims,options); 
% end

if options.Visualize{1} == true

    VisualizeResults(Robot,DH_params_init,DH_params,DH_param_lims,Info,options);
    
end


end