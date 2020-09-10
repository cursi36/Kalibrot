%% Analyzes data and returns observability and Confidence Measure
% dP = P_measured-P_approximated
% D = jacobian 
% x = parameter estimates
% n_pars = number parameters estimated
% m = number of measurements
% returns:
% conf = confidence measure of parameters estimates
% obs = observability measures
function [obs,conf] = AnalyzeData(dP,D,x,n_pars,m)
conf = [];

% chi_2 = (dP-D*x)'*(dP-D*x);
% v = m-n_pars;
% sigma_2 = chi_2/v;
% H = D'*D;
% conf = sigma_2*inv(H);

%observability measures
[~,S,~] = svd(D'*D);
sigma = sqrt(diag(S));
obs(1) = prod(sigma);

m = min(sigma);
if ( m > 1e-10)
    
else
    
    m = 1e-10;
end
obs(2) = max(sigma)/m;
obs(3) = m;

end