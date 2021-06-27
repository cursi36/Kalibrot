function VisualizeResults(Robot,DH_params_init,DH_params,DH_params_lims,Info,options)

joint_types = Robot.m_joint_types;
n_jnts = length(joint_types);

%%%%%ROBOT PLOTTING
sigma = zeros(length(joint_types),1);
f_prism = find(joint_types == 'p');
f_rev = find(joint_types == 'r');
sigma(f_prism) = 1;

%d, theta, a, alpha
DH_tab = reshape(DH_params,4,n_jnts)';
eps = 1e-15;
%convert to Robotic Toolbox convention
theta = eps*ones(n_jnts,1);
d = eps*ones(n_jnts,1);
offset = zeros(n_jnts,1);

offset(f_rev) = DH_tab (f_rev,2);
offset(f_prism) = DH_tab (f_prism,1);
theta(f_prism) = DH_tab (f_prism,2);
d = DH_tab (:,1);
a = DH_tab (:,3);
alpha = DH_tab (:,4);

% %theta,d,a,alpha,sigma (0 = rev) offset
dh = [theta,d,a,alpha,sigma,offset];

robot_RT = SerialLink(dh,'name',"Final Robot");
if isempty(f_prism) == 0
    robot_RT.qlim(f_prism,:) = [zeros(length(f_prism),1), 0.11*ones(length(f_prism),1)];
end
q = options.Visualize{2};
[n_r,~] = size(q);
if (n_r > 1)
    q = q';
end
if isempty(q) == 1
    q = zeros(1,n_jnts);
end
disp("PLOTTING ROBOT STRUCTURE")
figure()
robot_RT.plot(q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%DH PARAMETERS
disp("PLOTTING DH PARMAS")
n_vars = length(DH_params);
figure()
sgtitle("DH Parameters")
if options.solver == "qp"
    
    hold on
    plot(1:n_vars,DH_params,'ob','MarkerSize',10,'LineWidth',3)
    center = (DH_params_lims(:,1)+DH_params_lims(:,2))/2;
    delta = (DH_params_lims(:,2)-DH_params_lims(:,1))/2;
    errorbar(1:n_vars,center,delta,delta,...
        '.r','LineWidth',1,'MarkerSize',0.1)
    plot(1:n_vars,DH_params_init,'.k','MarkerSize',10)
    % plotBounds(DH_params,DH_params_init,DH_params_lims)
    legend("Final DH params","limits","Intial DH Params")
    
else
    hold on
    plot(1:n_vars,DH_params,'ob','MarkerSize',8,'LineWidth',3)
    plot(1:n_vars,DH_params_init,'.k','MarkerSize',8)
    legend("Final DH params","Intial DH Params")
    
end
box on
grid on
labels = cell(n_vars,1);
labels(1:4:n_vars) = {"d"};
labels(2:4:n_vars) = {"theta"};
labels(3:4:n_vars) = {"a"};
labels(4:4:n_vars) = {"alpha"};
xticks(1:n_vars)
xticklabels(labels)
xlim([0, n_vars+1])
xlabel("DH Params")
ylabel("DH Params Values")
set(gca,'fontname','Times New Romans')


%%%%%%%%%%%%%%%%%%%%%
if isempty(Info.PermMat) == 0
    disp("PLOTTING MATRIX OF BASE PARAMETERS")
    K = abs(Info.PermMat(:,:,end));
    f = find(K(:) < 1e-06);
    K(f) = 0;
    f = find(K(:) >= 1e-06);
    K(f) = 1;
    Color = flipud(gray);
    
    figure()
    sgtitle("Combination Matrix")
    colormap(Color)
    % pcolor(K)
    imagesc(K)
    colorbar
    ylabel("Base Params")
    xlabel("DH Params")
    grid on
    xticks(1:n_vars)
    yticks(1:n_vars)
    xticklabels(labels)
    hold on
    plot([1 n_vars],[1 n_vars],'--k','LineWidth',2)
    set(gca,'fontname','Times New Romans')
    axis equal
    
    pause(1)
    
end

end

function plotBounds(DH_params,DH_params_init,DH_params_lims)

n_vars = length(DH_params);

hold on
plot(1:n_vars,DH_params,'ob','MarkerSize',10,'LineWidth',3)

for i = 1:length(DH_params)
    plot([i i],[DH_params_lims(i,1) DH_params_lims(i,2)],'-r','LineWidth',1)
    plot([i-0.1 i+0.1],[DH_params_lims(i,1) DH_params_lims(i,1)],'-r','LineWidth',1)
    plot([i-0.1 i+0.1],[DH_params_lims(i,2) DH_params_lims(i,2)],'-r','LineWidth',1)
end

plot(1:n_vars,DH_params_init,'.k','MarkerSize',10)


end
