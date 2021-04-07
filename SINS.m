function [cnb,vn,pos,dRnm]=SINS(cnb_1,vn_1,vn_2,pos_1,pos_2,dRnm_1,dRnm_2,wm,vm,ts)
glvs;
Tm=4*ts;

%% 姿态
Lm_half = 1.5*pos_1(1)-0.5*pos_2(1); 
hm_half = 1.5*pos_1(3)-0.5*pos_2(3); 
wniem_half = [0;glv.Wie*cos(Lm_half);glv.Wie*sin(Lm_half)];
wnen = [-vn_1(2)/(glv.Re+hm_half);vn_1(1)/(glv.Re+hm_half);vn_1(1)*tan(Lm_half)/(glv.Re+hm_half)];
Fcm_half = [0,-1/(glv.Re+hm_half),0;1/(glv.Re+hm_half),0,0;tan(Lm_half)/(glv.Re+hm_half),0,0];
cge = [-sin(pos_1(2)),cos(pos_1(2)),0;-sin(pos_1(1))*cos(pos_1(2)),-sin(pos_1(1))*sin(pos_1(2)),cos(pos_1(1));cos(pos_1(1))*cos(pos_1(2)),cos(pos_1(1))*sin(pos_1(2)),sin(pos_1(1))];
Sm = wniem_half*Tm+Fcm_half*(2*dRnm_1-dRnm_2);
qnb_1 = Trans_att2quat(Trans_attm2att(cnb_1));
phi=wm(:,1)+wm(:,2)+wm(:,3)+wm(:,4)+736/945*(cross(wm(:,1),wm(:,2))+cross(wm(:,3),wm(:,4)))+334/945*(cross(wm(:,1),wm(:,3))+cross(wm(:,2),wm(:,4)))+526/945*cross(wm(:,1),wm(:,4))+654/945*cross(wm(:,2),wm(:,3));
psi_m=wniem_half*Tm+Fcm_half*(2*dRnm_1-dRnm_2);
Cn=eye(3)+sin(norm(psi_m))/norm(psi_m)*antisym_mat(psi_m)+(1-cos(norm(psi_m)))/(norm(psi_m))^2*antisym_mat(psi_m)*antisym_mat(psi_m);
Cb=eye(3)+sin(norm(phi))/norm(phi)*antisym_mat(phi)+(1-cos(norm(phi)))/(norm(phi))^2*antisym_mat(phi)*antisym_mat(phi);
Cnb=Cn'*Trans_quat2attm(qnb_1)*Cb;
cnb = Cn'*cnb_1*Cb;

%% 速度更新
vnm_half=1.5*vn_1-0.5*vn_2;
gm_half=glv.G;
gnm_half=[0;0;-gm_half];
dvngcorm=gnm_half*Tm-cross((Fcm_half*vnm_half+2*wniem_half),(2*dRnm_1)-dRnm_2);
dvbm = sum(vm,2)+0.5*cross(sum(wm,2),sum(vm,2))+cross((54/105*wm(:,1)+92/105*wm(:,2)+214/105*wm(:,3)),vm(:,4))+cross((54/105*vm(:,1)+92/105*vm(:,2)+214/105*vm(:,3)),wm(:,4));
dvnsfm=quat_mulVec(qnb_1,dvbm)-0.5*antisym_mat(Sm)*Trans_quat2attm(qnb_1)*sum(vm,2);
vn=vn_1+dvnsfm+dvngcorm;

%% 位置更新
dthet1=wm(:,1)+wm(:,2);
dthet2=wm(:,3)+wm(:,4);
dv1=vm(:,1)+vm(:,2);
dv2=vm(:,3)+vm(:,4);
A=1/Tm*(3*dv1-dv2);
B=4/Tm^2*(dv2-dv1);
Sdvm=Tm^2*A/2+Tm^3*B/6;
dRrotm=Tm*(cross(dthet1,(5/18*dv1+1/6*dv2))+cross(dthet2,(1/6*dv1+1/18*dv2)));
dRsculm=Tm*(cross(dthet1,(11/90*dv1+1/10*dv2))+cross(dthet2,(1/90*dv2-7/30*dv1)));
dRbsfmm_1=Sdvm+dRrotm+dRsculm;
dRnsfm=-1/3*antisym_mat(Sm)*quat_mulVec(qnb_1,dvbm)*Tm+quat_mulVec(qnb_1,dRbsfmm_1);
dRnm=(vn_1+0.5*dvngcorm)*Tm+dRnsfm;
zetam=Fcm_half*dRnm;
cne=eye(3)-antisym_mat(zetam);
cpos=cne*cge;
cpos=cpos*(cpos'*cpos)^(-1/2);
Latitude=asin(cpos(3,3))*180/pi;
longmain=atan(cpos(3,2)/cpos(3,1))*180/pi;
if longmain<0 && cpos(3,1)<0
    Longtitude=longmain+180;
end
if longmain>0 && cpos(3,1)<0
    Longtitude=longmain-180;
end
if cpos(3,1)>0
    Longtitude=longmain;
end
h=pos_1(3)+dRnm(3);
pos=[Latitude*pi/180;Longtitude*pi/180;h];