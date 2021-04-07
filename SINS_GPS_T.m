%% ENU
clear all;
clc;
glvs;
%% 导航原始数据
addpath('Data');
%24星位置
load('matG1_G24.mat'); 
%24星速度
load('matV1_V24.mat');
%24星伪距误差
load('matTu.mat');
%24星伪距率误差
load('matTru.mat');
%24星伪距
load('matPI.mat');
%24星伪距率
load('matPIs.mat');
%IMU数据
load('IMUdata.mat');
N = 4;                       
Ts = 0.01; 
An = 300/Ts;
TS = N*Ts;                  
Num=745/Ts;

IMU = imu(1:Num,1:6);                      
Wdata = [IMU(:,1) IMU(:,2) IMU(:,3)];    
Fdata = [IMU(:,4) IMU(:,5) IMU(:,6)];      
GPS = [imu(:,19) imu(:,20) imu(:,21)];     
PosR = [imu(:,13) imu(:,14) imu(:,15)];    
disp('Step.1:--- 数据导入结束 ---');

%% 初始对准
AIMU = imu(1:An,:);
[Apsi,Atheta,Agamma] = Align(AIMU);
disp('Step.2:--- 初始对准结束 ---');

%% 组合导航初始化 
L = length(IMU); 
Cnb = Trans_att2attm([Apsi,Atheta,Agamma]*glv.D2R);  
[Vn,Pos, vn_1,vn_2, pos_1,pos_2,posIn_1,posIn_2,Qt,P_K,X_K,xc] = SINS_initT(matTu,matTru);
countNa=1;
countFuse=1;
NUM = 5;
%% 组合导航
countNa = 1;
for count=1:N:L
    FmIn=Ts*(Fdata(count:count+N-1,:))';   
    WmIn=Ts*(Wdata(count:count+N-1,:))';    
    %导航解算
    [Cnb,Vn, Pos, posIn] = SINS(Cnb,vn_1, vn_2, pos_1, pos_2, posIn_1, posIn_2, WmIn, FmIn, Ts); 
    if rem(countNa,NUM)==0
        Fn=Cnb*Fdata(countNa*4,:)'; 
        Rm=glv.Re*(1-2*glv.e+3*glv.e*sin(Pos(1))^2);
        Rn=glv.Re*(1+glv.e*sin(Pos(1))^2);  
        
        Ft = KalmanPhi_T(Vn,Cnb,Pos,Fn,Rm,Rn);
        F_K = eye(17)+Ft*0.01*20;
        Q_K = diag([(0.005*glv.D2R/3600)^2;(0.005*glv.D2R/3600)^2;(0.005*glv.D2R/3600)^2;(5e-5)^2;(5e-5)^2;(5e-5)^2;0;0;0;0;0;0;0;0;0;(1)^2;(0.1)^2]); 
        X_K = F_K*X_K;
        P_K = F_K*P_K*F_K'+Q_K*Ts;
        
        [resZK,Cnb,Vn,Pos,P_K,X_K,xc] = Zk_update(Cnb,Vn,Pos,Rm,Rn,matPI,matPIs,matTu,matTru,matG1_G24,matV1_V24,count,xc,P_K,X_K,PosR(countNa*4,:));
        Na_res_Tulti(:,countFuse)=[Trans_attm2att(Cnb);Vn;Pos;X_K(10:15)];
        PosTTTulti(countFuse,:)=PosR(count,:);
        countFuse=countFuse+1;
    end    
    vn_2 = vn_1;                
    vn_1 = Vn;                   
    pos_2 = pos_1;                   pos_1 = Pos;                 
    posIn_2 = posIn_1;          
    posIn_1 = posIn;            
    countNa=countNa+1;          
end
disp('Step.3:--- 组合导航结束 ---');
%% 组合结果
plot3(PosTTTulti(:,2)*glv.R2D,PosTTTulti(:,1)*glv.R2D,PosTTTulti(:,3));
hold on;
plot3(Na_res_Tulti(8,:)'*glv.R2D,Na_res_Tulti(7,:)*glv.R2D,Na_res_Tulti(9,:));
grid on;
xlabel('经度/°');ylabel('纬度/°');

lenRes = [1:length(Na_res_Tulti)];
figure name '位置误差'
subplot(311);plot(lenRes,Na_res_Tulti(7,:)*glv.R2D-PosTTTulti(:,1)'*glv.R2D,'r'); 
xlabel('t/s');ylabel('纬度误差/°');  title('纬度误差');  
grid on; 
subplot(312);plot(lenRes,Na_res_Tulti(8,:)*glv.R2D-PosTTTulti(:,2)'*glv.R2D,'r');
xlabel('t/s');ylabel('经度误差/°');  title('经度误差');   
grid on;
subplot(313);plot(lenRes,Na_res_Tulti(9,:)-PosTTTulti(:,3)','r'); 
xlabel('t/s');ylabel('高程误差/m');  title('高程误差'); 
grid on;   