
classdef RobotKinematics
    
    properties
        
        m_n_joints = [];
        m_joint_types = [];
        m_T_init = [];
        m_DH_params = [];
        m_q = [];
        
        m_T_sym = [];
        m_Dp_sym = [];  %derivative of position wrt DH params
        m_Dor1_sym = []; %derivative of oritnt (condition 1) wrt DH params
        m_Dor2_sym = []; %derivative of oritnt (condition 2) wrt DH params
        m_Dor3_sym = []; %derivative of oritnt (condition 1) wrt DH params
        m_Dor4_sym = []; %derivative of oritnt (condition 2) wrt DH params
        
        m_orient_cond = []
        
    end
    
    methods
        
        function obj = RobotKinematics(n_joints, types, T_init)
            
            obj.m_n_joints = n_joints;
            obj.m_joint_types = types;
            obj.m_T_init = T_init;
            
            obj.m_DH_params = sym('m_DH_params',[4*n_joints,1]);
            obj.m_q = sym('m_q',[n_joints,1]);
            
            [obj.m_T_sym,obj.m_Dp_sym, ...
                obj.m_Dor1_sym,obj.m_Dor2_sym,obj.m_Dor3_sym,obj.m_Dor4_sym] = getKineSym(obj);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% substitute values in symbolic functions
        function [obj,T_val,P] = getPose(obj,q,DH_params)
            
            T_val =  subs(obj.m_T_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            T_val = double(T_val);
            
            R = T_val(1:3,1:3);
            
            [quat,obj.m_orient_cond] = Rot2Quat(R);
                        
            P = [T_val(1:3,4);quat];
            
        end
        
        %% substitue in symbolic function
        function [Dp,Dor] = getDerivs(obj,q,DH_params)
            
             Dp = subs(obj.m_Dp_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            
            if obj.m_orient_cond == 1

                Dor = subs(obj.m_Dor1_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            elseif obj.m_orient_cond == 2

                Dor = subs(obj.m_Dor2_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            
            elseif obj.m_orient_cond == 3

                Dor = subs(obj.m_Dor3_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            
            elseif obj.m_orient_cond == 4

                Dor = subs(obj.m_Dor4_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
            end
            
            Dp = double(Dp);
            Dor = double(Dor);
        end
        
        %% Symbolic FWD Kine
        function [T,Dp,Dor1,Dor2,Dor3,Dor4] = getKineSym(obj)
            
            T = obj.m_T_init;
            T = sym(T);
            
            for i = 1:obj.m_n_joints
                
                v1 = 4*i-3;
                xi = obj.m_DH_params(v1:v1+3);
                
                if obj.m_joint_types(i) == 'r'
                    
                    xi(2) = xi(2)+obj.m_q(i);
                else
                    xi(1) = xi(1)+obj.m_q(i);
                end
                Ti = DH_mat(xi);
                
                T = T*Ti;
                
            end
            
            %             T = simplify(T);
            
            P = T(1:3,4);
            R = T(1:3,1:3);
            quat_1 = Rot2Quat_1(R);
            quat_2 = Rot2Quat_2(R);
            quat_3 = Rot2Quat_3(R);
            quat_4= Rot2Quat_4(R);
            
            Dp = jacobian(P,obj.m_DH_params);

            Dor1 = jacobian(quat_1,obj.m_DH_params);
            Dor2 = jacobian(quat_2,obj.m_DH_params);
            Dor3 = jacobian(quat_3,obj.m_DH_params);
            Dor4 = jacobian(quat_4,obj.m_DH_params);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%% NUMERICAL COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get pose numerical. Do not use symbolic functions
        function [obj,T_val,P] = getPoseNum(obj,q,DH_params)
            
            T_val = obj.m_T_init;
            
            for i = 1:obj.m_n_joints
                
                v1 = 4*i-3;
                xi = DH_params(v1:v1+3);
                
                if obj.m_joint_types(i) == 'r'
                    
                    xi(2) = xi(2)+q(i);
                else
                    xi(1) = xi(1)+q(i);
                end
                Ti = DH_matNum(xi);
                
                T_val = T_val*Ti;
                
            end
            R = T_val(1:3,1:3);
            
            [quat,~] = Rot2Quat(R);           
            P = [T_val(1:3,4);quat];
            
        end
        
        %% Position Derivative numerical
        function [D] = getPDerivNum(obj,q,DH_params,P_e)
            
            n_var = length(DH_params);
            D = zeros(3,n_var);
            
            if isempty(P_e) == 0
                
            else
                [~,~,P_e] = getPoseNum(obj,q,DH_params);
                P_e = P_e(1:3,1); %%ee position in 0
                
            end
            
            T = obj.m_T_init;
            
            for i = 1:obj.m_n_joints
                
                R_i_1_0 = T(1:3,1:3); %R_i-1_0
                
                v1 = 4*i-3;
                xi = DH_params(v1:v1+3);
                
                type = obj.m_joint_types(i);
                [DP_i,dR_i] = Pose_deriv(q(i),xi,type); % of i wrt i-1
                
                if type == 'r'
                    
                    xi(2) = xi(2)+q(i);
                else
                    xi(1) = xi(1)+q(i);
                end
                
                Ti = DH_matNum(xi);
                
                T = T*Ti; %T_i_0
                R_i_0 = T(1:3,1:3); %R_i_0
                
                P_i = T(1:3,4);
                P_e_i = R_i_0'*(P_e-P_i);
                DR_i = zeros(3,4);
                
                for dim = 1:3
                    
                    DR_i(dim,:) = P_e_i'*dR_i(:,:,dim);

                end
                Di = R_i_1_0*(DP_i+DR_i);
                
                D(:, v1:v1+3) = Di;
                
            end
            
        end
        
        %% Quaternion Derivative numerical
        function [Dquat,quat] = getQuatDerivNum(obj,q,DH_params)
            
            n_var = length(DH_params);
            
            T = obj.m_T_init;
            
            R_i = T(1:3,1:3);
            
            [quat,~] = Rot2Quat(R_i); 
                        
            Dquat = zeros(4,n_var);
            
            for i = 1:obj.m_n_joints
                
                v1 = 4*i-3;
                xi = DH_params(v1:v1+3);
                
                type = obj.m_joint_types(i);
                dR_i = RotMatDeriv(q(i),xi,type); %%in R^9x4
                
                dR = zeros(9,n_var);
                dR(:,v1:v1+3) = dR_i;
                
                if type == 'r'
                    
                    xi(2) = xi(2)+q(i);
                else
                    xi(1) = xi(1)+q(i);
                end
                
                Ti = DH_matNum(xi);
                
                R_i = Ti(1:3,1:3);
                
                [quat_i,cond] = Rot2Quat(R_i); 
                
                if cond == 1 

                    dq_dR_i = quatDeriv_Rot_1(R_i); %% in R^4x9
                    
                elseif cond == 2

                    dq_dR_i = quatDeriv_Rot_2(R_i);
                    
                elseif cond == 3

                    dq_dR_i = quatDeriv_Rot_3(R_i);
                    
                elseif cond == 4

                    dq_dR_i = quatDeriv_Rot_4(R_i);
                    
                end
                
                n = quat(4,1);
                e = quat(1:3,1);
                
                n_i = quat_i(4,1);
                e_i = quat_i(1:3,1);
                
                quat(1:3,1) = n*e_i+n_i*e+cross(e,e_i);
                quat(4,1) = n*n_i-e'*e_i;
                
                Dn = Dquat(4,:);
                De = Dquat(1:3,:);
                
                Dn_i = dq_dR_i(4,:)*dR;
                De_i =  dq_dR_i(1:3,:)*dR;
                
                %derivative of qw
                Dquat(4,:) = Dn*n_i+n*Dn_i;
                der = zeros(1,n_var);
                for dim = 1:3
                    der = der+De(dim,:)*e_i(dim)+e(dim)*De_i(dim,:);
                end
                Dquat(4,:) = Dquat(4,:)-der;
                
                %derivative of qx,y,z
                Dquat(1:3,:) = n*De_i+n_i*De;
                der = zeros(3,n_var);
                for dim = 1:3
                    der(dim,:) = Dn*e_i(dim)+Dn_i*e(dim);
                end
                Dquat(1:3,:) = Dquat(1:3,:)+der;
                
                %%cross prodcut deriv
                Dcp(1,:) = -De(3,:)*e_i(2)+De(2,:)*e_i(3);
                Dcp(2,:) = De(3,:)*e_i(1)-De(1,:)*e_i(3);
                Dcp(3,:) = -De(2,:)*e_i(1)+De(1,:)*e_i(2);
                
                Dquat(1:3,:) = Dquat(1:3,:)+Dcp;
                
                Dcp(1,:) = -De_i(3,:)*e(2)+De_i(2,:)*e(3);
                Dcp(2,:) = De_i(3,:)*e(1)-De_i(1,:)*e(3);
                Dcp(3,:) = -De_i(2,:)*e(1)+De_i(1,:)*e(2);
                
                Dquat(1:3,:) = Dquat(1:3,:)-Dcp;
            end
            
        end
        
        %%GetPoseDerivNum: no symbolic used
        %%Get Pose, DP, Dquat derivs together
        
        function [obj,T,P,Dp,Dquat] = getKineDeriv_Ana(obj,q,DH_params)
            
            n_var = length(DH_params);
            
            T = obj.m_T_init;
            
            R_i = T(1:3,1:3);
                        
            [quat,~] = Rot2Quat(R_i); 
            
            Dquat = zeros(4,n_var);
            
            for i = 1:obj.m_n_joints
                
                v1 = 4*i-3;
                xi = DH_params(v1:v1+3);
                
                type = obj.m_joint_types(i);
                dR_i = RotMatDeriv(q(i),xi,type); %%in R^9x4
                
                dR = zeros(9,n_var);
                dR(:,v1:v1+3) = dR_i;
                
                if type == 'r'
                    
                    xi(2) = xi(2)+q(i);
                else
                    xi(1) = xi(1)+q(i);
                end
                
                Ti = DH_matNum(xi);
                
                T = T*Ti; %%FKine
                
                %%%%%%Quat Derivatives
                R_i = Ti(1:3,1:3);
                
                [quat_i,cond] = Rot2Quat(R_i); 
                
                if cond == 1 

                    dq_dR_i = quatDeriv_Rot_1(R_i); %% in R^4x9
                    
                elseif cond == 2

                    dq_dR_i = quatDeriv_Rot_2(R_i);
                    
                elseif cond == 3

                    dq_dR_i = quatDeriv_Rot_3(R_i);
                    
                elseif cond == 4

                    dq_dR_i = quatDeriv_Rot_4(R_i);
                    
                end
                
                n = quat(4,1);
                e = quat(1:3,1);
                
                n_i = quat_i(4,1);
                e_i = quat_i(1:3,1);
                
                quat(1:3,1) = n*e_i+n_i*e+cross(e,e_i);
                quat(4,1) = n*n_i-e'*e_i;
                
                Dn = Dquat(4,:);
                De = Dquat(1:3,:);
                
                Dn_i = dq_dR_i(4,:)*dR;
                De_i =  dq_dR_i(1:3,:)*dR;
                
                %derivative of qw
                Dquat(4,:) = Dn*n_i+n*Dn_i;
                der = zeros(1,n_var);
                for dim = 1:3
                    der = der+De(dim,:)*e_i(dim)+e(dim)*De_i(dim,:);
                end
                Dquat(4,:) = Dquat(4,:)-der;
                
                %derivative of qx,y,z
                Dquat(1:3,:) = n*De_i+n_i*De;
                der = zeros(3,n_var);
                for dim = 1:3
                    der(dim,:) = Dn*e_i(dim)+Dn_i*e(dim);
                end
                Dquat(1:3,:) = Dquat(1:3,:)+der;
                
                %%cross prodcut deriv
                Dcp(1,:) = -De(3,:)*e_i(2)+De(2,:)*e_i(3);
                Dcp(2,:) = De(3,:)*e_i(1)-De(1,:)*e_i(3);
                Dcp(3,:) = -De(2,:)*e_i(1)+De(1,:)*e_i(2);
                
                Dquat(1:3,:) = Dquat(1:3,:)+Dcp;
                
                Dcp(1,:) = -De_i(3,:)*e(2)+De_i(2,:)*e(3);
                Dcp(2,:) = De_i(3,:)*e(1)-De_i(1,:)*e(3);
                Dcp(3,:) = -De_i(2,:)*e(1)+De_i(1,:)*e(2);
                
                Dquat(1:3,:) = Dquat(1:3,:)-Dcp;
            end
            
            P(1:3,1) = T(1:3,4);
            P(4:7,:) = quat;
            
            P_e = P(1:3,1);
            [Dp] = obj.getPDerivNum(q,DH_params,P_e);
            
            
        end
        
        
        
%         %% get derivatives with Dp computed numerically without symbolic
%         function [Dp,Dor] = getDerivsNum(obj,q,DH_params)
%             
%             Dp =  getPDerivNum(obj,q,DH_params);
%             
%             if obj.m_orient_cond == 1
%                 Dor = subs(obj.m_Dor1_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
%             else
%                 Dor = subs(obj.m_Dor2_sym,[obj.m_DH_params;obj.m_q],[DH_params;q]);
%             end
%             Dor = double(Dor);
%         end
        
        
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = DH_matNum(x)
d = x(1);
theta = x(2);
a = x(3);
alpha = x(4);

T = zeros(4,4);
T(4,4) = 1;

T(1,1) = cos(theta);    T(1,2) = -sin(theta)*cos(alpha);    T(1,3) = sin(theta)*sin(alpha); T(1,4) = a*cos(theta);
T(2,1) = sin(theta);    T(2,2) = cos(theta)*cos(alpha); T(2,3) = -cos(theta)*sin(alpha);    T(2,4) = a*sin(theta);
T(3,2) = sin(alpha);    T(3,3) = cos(alpha);    T(3,4) = d;

end

function T = DH_mat(x)
d = x(1);
theta = x(2);
a = x(3);
alpha = x(4);

T = zeros(4,4);
T(4,4) = 1;

T = sym(T);
T(1,1) = cos(theta);    T(1,2) = -sin(theta)*cos(alpha);    T(1,3) = sin(theta)*sin(alpha); T(1,4) = a*cos(theta);
T(2,1) = sin(theta);    T(2,2) = cos(theta)*cos(alpha); T(2,3) = -cos(theta)*sin(alpha);    T(2,4) = a*sin(theta);
T(3,2) = sin(alpha);    T(3,3) = cos(alpha);    T(3,4) = d;

end

% quat = [x,y,z,w]
% if 1+R(1,1)-R(2,2)-R(3,3) > 0
function quat = Rot2Quat_1(R)
%
%https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Quaternions
t = R(1,1)-R(2,2)-R(3,3);

quat(1,1) = 0.5*sqrt(1+t); %%qx
s = 1/(4*quat(1));

quat(2,1) = s*(R(1,2)+R(2,1)); %%qy
quat(3,1) = s*(R(1,3)+R(3,1)); %%qz
quat(4,1) = s*(R(3,2)-R(2,3)); %%qw

end

% quat = [x,y,z,w]
% if 1+R(1,1)+R(2,2)+R(3,3) > 0
function quat = Rot2Quat_2(R)
%
%https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Quaternions
t = R(1,1)+R(2,2)+R(3,3);

quat(4,1) = 0.5*sqrt(1+t); %%qw
s = 1/(4*quat(4));

quat(1,1) = s*(R(3,2)-R(2,3)); %%qx
quat(2,1) = s*(R(1,3)-R(3,1)); %%qy
quat(3,1) = s*(R(2,1)-R(1,2)); %%qz

end

function quat = Rot2Quat_3(R)
%
%http://docs.ros.org/indigo/api/orocos_kdl/html/frames_8cpp_source.html#l00205t = R(1,1)+R(2,2)+R(3,3);

t = R(2,2)-R(1,1)-R(3,3);

quat(2,1) = 0.5*sqrt(1+t); %%qw
s = 1/(4*quat(2));

quat(1,1) = s*(R(1,2)+R(2,1)); %%qx
quat(3,1) = s*(R(2,3)+R(3,2)); %%qz
quat(4,1) = s*(R(1,3)-R(3,1)); %%qw

end

function quat = Rot2Quat_4(R)
%
%http://docs.ros.org/indigo/api/orocos_kdl/html/frames_8cpp_source.html#l00205t = R(1,1)+R(2,2)+R(3,3);

t = R(3,3)-R(2,2)-R(1,1);

quat(3,1) = 0.5*sqrt(1+t); %%qw
s = 1/(4*quat(3));

quat(1,1) = s*(R(1,3)+R(3,1)); %%qx
quat(2,1) = s*(R(2,3)+R(3,2)); %%qz
quat(4,1) = s*(R(2,1)-R(1,2)); %%qw

end

function [quat,condition] = Rot2Quat(R)

t = R(1,1)+R(2,2)+R(3,3);

if (t > 1e-12)
    
    quat = Rot2Quat_2(R);
    condition = 2;
    
else
   
    if (R(1,1) > R(2,2) && R(1,1) > R(3,3))
        
        quat = Rot2Quat_1(R);
        condition = 1;
        
    elseif (R(2,2) > R(3,3))
        
        quat = Rot2Quat_3(R);
        condition = 3;
        
    else 
        quat = Rot2Quat_4(R);
        condition = 4;
    end
    
end
end

function R = Quat2Rot(quat)
x = quat(1);    y = quat(2);    z = quat(3); w = quat(4);

R(1,1) = 1-2*y^2-2*z^2;
R(1,2) = 2*(x*y-z*w);
R(1,3) = 2*(x*z+y*w);

R(2,1) = 2*(x*y+z*w);
R(2,2) = 1-2*x^2-z^2;
R(2,3) = 2*(y*z-x*w);

R(3,1) = 2*(x*z-y*w);
R(3,2) = 2*(y*z+x*w);
R(3,3) = 1-x^2-y^2;
end


%%%%DERIVATIVES NUmerical

% derivative of P_i+1_i and R_i+1_i wrt DH_param of link i
function [dP,dR] = Pose_deriv(q_i,DH_params_i,type)

dP = zeros(3,4);
dR = zeros(3,4,3); %% for each row of R

if type == 'r'
    
    DH_params_i(2) = DH_params_i(2)+q_i;
    
else
    
    DH_params_i(1) = DH_params_i(1)+q_i;
    
    
end
d = DH_params_i(1);   theta = DH_params_i(2); a = DH_params_i(3);   alpha = DH_params_i(4);
% dd,dtheta,da,dalpha
dP(1,:) =  [0 -a*sin(theta) cos(theta) 0];
dP(2,:) =  [0 a*cos(theta) sin(theta) 0];
dP(3,:) =  [1 0 0 0];

%derivative of first row
% [dR11/dx; dR12/dx; dR13/dx]
dR(:,:,1) = ...
    [0 -sin(theta) 0 0;
    0 -cos(theta)*cos(alpha) 0 sin(theta)*sin(alpha);
    0 cos(theta)*sin(alpha) 0 sin(theta)*cos(alpha)];

% [dR21/dx; dR22/dx; dR23/dx]
dR(:,:,2) = ...
    [0 cos(theta) 0 0;
    0 -sin(theta)*cos(alpha) 0 -cos(theta)*sin(alpha);
    0 sin(theta)*sin(alpha) 0 -cos(theta)*cos(alpha)];

% [dR31/dx; dR32/dx; dR33/dx]
dR(:,:,3) = ...
    [0 0 0 0;
    0 0 0 cos(alpha);
    0 0 0 -sin(alpha)];

end

% quaternion derivatives wrt to rotation matrix
% d in R^4x9
function dq_dR = quatDeriv_Rot_1(R)

dq_dR = zeros(4,9);

t = 1+R(1,1)-R(2,2)-R(3,3);

quat(1,1) = 0.5*sqrt(t); %%qx

if quat(1,1) < 1e-10
    quat(1,1) =  1e-10;
end

dq_dR(1,:) = 1/4*t^-(1/2)*[1 0 0 0 -1 0 0 0 -1];

dq_dR(2,:) = -1/(4*quat(1,1)^2)*(R(1,2)+R(2,1))*dq_dR(1,:)+...
    1/(4*quat(1,1))*[0 1 0 1 0 0 0 0 0];

dq_dR(3,:) = -1/(4*quat(1,1)^2)*(R(1,3)+R(3,1))*dq_dR(1,:)+...
    1/(4*quat(1,1))*[0 0 1 0 0 0 1 0 0];

dq_dR(4,:) = -1/(4*quat(1,1)^2)*(R(3,2)-R(2,3))*dq_dR(1,:)+...
    1/(4*quat(1,1))*[0 0 0 0 0 -1 0 1 0];

end


function dq_dR = quatDeriv_Rot_2(R)

dq_dR = zeros(4,9);

t = 1+R(1,1)+R(2,2)+R(3,3);

quat(4,1) = 0.5*sqrt(t); %%qw

if quat(4,1) < 1e-10
    quat(4,1) =  1e-10;
end


dq_dR(4,:) = 1/4*t^-(1/2)*[1 0 0 0 1 0 0 0 1];

dq_dR(1,:) = -1/(4*quat(4,1)^2)*(R(3,2)-R(2,3))*dq_dR(4,:)+...
    1/(4*quat(4,1))*[0 0 0 0 0 -1 0 1 0];

dq_dR(2,:) = -1/(4*quat(4,1)^2)*(R(1,3)-R(3,1))*dq_dR(4,:)+...
    1/(4*quat(4,1))*[0 0 1 0 0 0 -1 0 0];

dq_dR(3,:) = -1/(4*quat(4,1)^2)*(R(2,1)-R(1,2))*dq_dR(4,:)+...
    1/(4*quat(4,1))*[0 -1 0 1 0 0 0 0 0];

end

function dq_dR = quatDeriv_Rot_3(R)

dq_dR = zeros(4,9);

t = 1+R(2,2)-R(1,1)-R(3,3);

quat(2,1) = 0.5*sqrt(t); %%qw

if quat(2,1) < 1e-10
    quat(2,1) =  1e-10;
end


dq_dR(2,:) = 1/4*t^-(1/2)*[-1 0 0 0 1 0 0 0 -1];

dq_dR(1,:) = -1/(4*quat(2,1)^2)*(R(1,2)+R(2,1))*dq_dR(2,:)+...
    1/(4*quat(2,1))*[0 1 0 1 0 0 0 0 0];

dq_dR(3,:) = -1/(4*quat(2,1)^2)*(R(2,3)+R(3,2))*dq_dR(2,:)+...
    1/(4*quat(2,1))*[0 0 0 0 0 1 0 1 0];

dq_dR(4,:) = -1/(4*quat(2,1)^2)*(R(1,3)-R(3,1))*dq_dR(2,:)+...
    1/(4*quat(2,1))*[0 0 1 0 0 0 -1 0 0];

end

function dq_dR = quatDeriv_Rot_4(R)

dq_dR = zeros(4,9);

t = 1+R(3,3)-R(1,1)-R(2,2);

quat(3,1) = 0.5*sqrt(t); %%qw

if quat(3,1) < 1e-10
    quat(3,1) =  1e-10;
end


dq_dR(3,:) = 1/4*t^-(1/2)*[-1 0 0 0 -1 0 0 0 1];

dq_dR(1,:) = -1/(4*quat(3,1)^2)*(R(1,3)+R(3,1))*dq_dR(3,:)+...
    1/(4*quat(3,1))*[0 0 1 0 0 0 1 0 0];

dq_dR(2,:) = -1/(4*quat(3,1)^2)*(R(2,3)+R(3,2))*dq_dR(3,:)+...
    1/(4*quat(3,1))*[0 0 0 0 0 1 0 1 0];

dq_dR(4,:) = -1/(4*quat(3,1)^2)*(R(2,1)-R(1,2))*dq_dR(3,:)+...
    1/(4*quat(3,1))*[0 -1 0 1 0 0 0 0 0];

end




% Rot mat deriv wrt DH params in vector form
function  dR = RotMatDeriv(q_i,DH_params_i,type)

dR = zeros(9,4); %% for each row of R

if type == 'r'
    
    DH_params_i(2) = DH_params_i(2)+q_i;
    
else
    
    DH_params_i(1) = DH_params_i(1)+q_i;
    
    
end
d = DH_params_i(1);   theta = DH_params_i(2); a = DH_params_i(3);   alpha = DH_params_i(4);

%derivative of first row
% [dR11/dx; dR12/dx; dR13/dx]
dR_mat(:,:,1) = ...
    [0 -sin(theta) 0 0;
    0 -cos(theta)*cos(alpha) 0 sin(theta)*sin(alpha);
    0 cos(theta)*sin(alpha) 0 sin(theta)*cos(alpha)];

% [dR21/dx; dR22/dx; dR23/dx]
dR_mat(:,:,2) = ...
    [0 cos(theta) 0 0;
    0 -sin(theta)*cos(alpha) 0 -cos(theta)*sin(alpha);
    0 sin(theta)*sin(alpha) 0 -cos(theta)*cos(alpha)];

% [dR31/dx; dR32/dx; dR33/dx]
dR_mat(:,:,3) = ...
    [0 0 0 0;
    0 0 0 cos(alpha);
    0 0 0 -sin(alpha)];

dR = [dR_mat(:,:,1);dR_mat(:,:,2);dR_mat(:,:,3)];

end



