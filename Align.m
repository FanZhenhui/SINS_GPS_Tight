function [ psi, theta, gamma ] = Align(IMUdata)

%% Init
glvs;
Weie = [0;0;7.292115e-5];
g_unit = 9.7803;
deltaT = 10e-3;
g0 = 9.7803;
Len_IMUdata = length(IMUdata);
lambda = IMUdata(1,20);
L = IMUdata(1,19);
H = IMUdata(1,21);
g = g0*(1 + 0.0052884*(sin(L))^2 - 0.0000059*(sin(2*L))^2) - 0.0003086*H;
RM = glv.Re*(1 - 2*glv.e + 3*glv.e*(sin(L))^2);
RN = glv.Re*(1 + glv.e*(sin(L))^2);
Cne = [-sin(lambda)        cos(lambda)         0;...
       -sin(L)*cos(lambda) -sin(L)*sin(lambda) cos(L);...
       cos(L)*cos(lambda)  cos(L)*sin(lambda)  sin(L)];

Wbib = [];
Fb = [];
for i = 1:2500
    Wbib = [Wbib [IMUdata(i,1);IMUdata(i,2);IMUdata(i,3)]]; 
    Fb = [Fb [IMUdata(i,4);IMUdata(i,5);IMUdata(i,6)]];
end
Wbib = sum(Wbib,2)./size(Wbib,2);
Fb = sum(Fb,2)./size(Fb,2);
Fb = -Fb;

%% 粗对准
Vn_T_inv = [0 0 1/(g*Weie(3)*cos(L));tan(L)/g 1/(Weie(3)*cos(L)) 0;-1/g 0 0];
Vb_T = [Fb(1) Fb(2) Fb(3);Wbib(1) Wbib(2) Wbib(3);...
    Fb(2)*Wbib(3) - Fb(3)*Wbib(2) Fb(3)*Wbib(1) - Fb(1)*Wbib(3) Fb(1)*Wbib(2) - Fb(2)*Wbib(1)];
Tnb = Vn_T_inv*Vb_T;
[psi,theta,gamma] = Trans_attm2att_ptg(Tnb);
Tnb = Trans_att2attm_ptg(psi,theta,gamma);

Wbib = [];
Fb = [];
for i = 1:2500
    Wbib = [Wbib [IMUdata(i,1);IMUdata(i,2);IMUdata(i,3)]]; 
    Fb = [Fb [IMUdata(i,4);IMUdata(i,5);IMUdata(i,6)]];
end
Wbib = sum(Wbib,2)./size(Wbib,2);
Fb = sum(Fb,2)./size(Fb,2);
Fb = -Fb;
Fn = Tnb*Fb;
Wnib = Tnb*Wbib;
phi = [-Fn(2)/g;Fn(1)/g;Wnib(1)/(Weie(3)*cos(L)) + Fn(1)*tan(L)/g];
Mphi = [0 -phi(3) phi(2);phi(3) 0 -phi(1);-phi(2) phi(1) 0];
Tnb = (eye(3,3) + Mphi)*Tnb;
[psi,theta,gamma] = Trans_attm2att_ptg(Tnb);
anttitude_res = [psi theta gamma];

%% 精对准
V_delta = 0.01;
W_epsilon = glv.D2R*(0.5)/3600;
W_d = glv.D2R*(0.5)/3600;
F_delta = 1e-5*g_unit;
F_d = 1e-5*g_unit;
PHI = glv.D2R*(1);

%数值初始化
Xk = [0 0 0 0 0 0 0 0 0 0]';
PHIx_0 = Xk(3);
f_PHIx_0 = PHIx_0;
Pk = diag([V_delta V_delta PHI PHI PHI F_delta F_delta W_epsilon W_epsilon W_epsilon].^2);
Q = diag([F_d F_d W_d W_d W_d 0 0 0 0 0].^2);
R = diag([V_delta V_delta].^2);

%系统矩阵
F = zeros(10,10);
F(1,2) = 2*Weie(3)*sin(L);
F(2,1) = -2*Weie(3)*sin(L);
F(3,4) = Weie(3)*sin(L);
F(3,5) = -Weie(3)*cos(L);
F(4,3) = -Weie(3)*sin(L);
F(5,3) = Weie(3)*cos(L);
F(1,4) = -g;
F(2,3) = g;
F(3,2) = -1/(RM + H);
F(4,1) = 1/(RN + H);
F(5,1) = tan(L)/(RN + H);
F(1:2,6:7) = Tnb(1:2,1:2);
F(3:5,8:10) = Tnb;
G = eye(10);
Hk = zeros(2,10);
Hk(1,1) = 1;
Hk(2,2) = 1;
PHIk_k_1 = eye(10) + F.*deltaT + F^2.*(deltaT^2/2);
GAMMAk_1 = deltaT.*(eye(10) + F.*(deltaT/2) + F^2.*(deltaT^2/6))*G;

Vn = [0;0;0];
Fb = -Fb;

RM = glv.Re*(1 - 2*glv.e + 3*glv.e*(sin(L))^2);
RN = glv.Re*(1 + glv.e*(sin(L))^2);
Wnen = [-Vn(2)/(RM + H);Vn(1)/(RN + H);Vn(1)/(RN + H)*tan(L)];

d_Vn_0 = d_V_N(Tnb,Fb,Wnen,Cne,Vn,g);

prog = 1;
PHI_res = [];
Weie = [0;0;7.292115e-5];
while((prog <= Len_IMUdata))

    Fb = [IMUdata(prog,4);IMUdata(prog,5);IMUdata(prog,6)];
    d_Vn = d_V_N(Tnb,Fb,Wnen,Cne,Vn,g);
    

    d_Vn = Tnb*Fb - (2.*antisym_mat(Cne*Weie) + antisym_mat(Wnen))*Vn - [0;0;g];
    Vn(1) = R_K_2(deltaT,Vn(1),d_Vn_0(1),d_Vn(1));
    Vn(2) = R_K_2(deltaT,Vn(2),d_Vn_0(2),d_Vn(2));
    Vn = [0;0;0];
    d_Vn_0 = d_Vn;
    
    RM = glv.Re*(1 - 2*glv.e + 3*glv.e*(sin(L))^2);
    RN = glv.Re*(1 + glv.e*(sin(L))^2);
    Wnen = [-Vn(2)/(RM + H);Vn(1)/(RN + H);Vn(1)/(RN + H)*tan(L)];
    
    Pk_k_1 = PHIk_k_1*Pk*PHIk_k_1' + GAMMAk_1*Q*GAMMAk_1';
    Kk = Pk_k_1*Hk'/(Hk*Pk_k_1*Hk' + R);
    Pk = (eye(10) - Kk*Hk)*Pk_k_1*(eye(10) - Kk*Hk)' + Kk*R*Kk';
    Xk_k_1 = PHIk_k_1*Xk;
    Xk = Xk_k_1 + Kk*(Vn(1:2,:) - Hk*Xk_k_1);
    
    f_PHIx = Filter_PHIx(PHIx_0,Xk(3),f_PHIx_0);
    d_phi_x = (f_PHIx - f_PHIx_0)./deltaT;
    PHIx_0 = Xk(3);
    f_PHIx_0 = f_PHIx;
    
    %计算姿态角
    phi_z = (-d_phi_x + Xk(4)*Weie(3)*sin(L) - Xk(2)/RM)/(Weie(3)*cos(L));
    phi = [Xk(3:4);phi_z];
    Tnb_k = (eye(3,3) + antisym_mat(phi))*Tnb;
    [psi,theta,gamma] = Trans_attm2att_ptg(Tnb_k);
    anttitude_res = [anttitude_res;[psi theta gamma]];
    
    PHI_res = [PHI_res [Xk(3:5);phi_z]];
    prog = prog + 1;
end

figure;
anttitude_res = glv.R2D*(anttitude_res);
index = (0:size(anttitude_res,1) - 1).*deltaT;
PHI_res = glv.R2D*(PHI_res); 
end

%% F3
function [ Yk ] = Filter_PHIx( Xk_1, Xk, Yk_1 )
%对水平误差角PHIx进行一阶低通滤波（Ref：房建成.一种新的惯导系统静基座快速初始对准方法[J].北京航空航天大学学报,1999(06):728-731.）
T = 10e-3;
Wdc = 0.001256;
Wac = tan(Wdc*T/2)*2/T;

Yk = Wac*T/2/(1 + Wac*T/2)*(Xk + Xk_1) + (1 - Wac*T/2)/(1 + Wac*T/2)*Yk_1;
end