function [resZK,Cnb,Vn,Pos,P_K,X_K,xc] = Zk_update(Cnb,Vn,Pos,Rm,Rn,matPI,matPIs,matTu,matTru,matG1_G24,matV1_V24,count,xc,P_K,X_K,PosR)
    glvs;
    e = glv.e;
    sL  = sin(Pos(1)); cL = cos(Pos(1));  sL2 = sL^2;   
    RMh = Rm + Pos(3); RNh = Rn + Pos(3);
    RN = Rn;
    sJ  = sin(Pos(2)); cJ = cos(Pos(2)); H = Pos(3);
    xI  = [RNh*cL*cJ;RNh*cL*sJ;(Rn*(1-e)^2+H)*sL];
    
    %-- Dr.LU
    D   = [-RNh*sL*cJ,       -RNh*cL*sJ, cL*cJ;
           -RNh*sL*sJ,       RNh*cL*cJ,  cL*sJ; 
           (RN*(1-e)^2+H)*cL,0,          sL    ];
    CneT = [-sL*cJ,cL*cJ,-sJ;-sL*sJ,cL*sJ,cJ;cL,sL,0]; 
    Cne  = [CneT(3,3),CneT(3,1),CneT(3,2);
            CneT(1,3),CneT(1,1),CneT(1,2);
            CneT(2,3),CneT(2,1),CneT(2,2)];
    Vnt = [Vn(1);Vn(2);Vn(3)];
    dxI = Cne*Vn;
    F  = [-Vnt(1)*cL*cJ-Vnt(2)*sL*cJ,-Vnt(3)*cJ+Vnt(1)*sL*sJ-Vnt(2)*cL*sJ,0;
            -Vnt(1)*cL*sJ-Vnt(2)*sL*sJ,-Vnt(3)*sJ-Vnt(1)*sL*cJ+Vnt(2)*cL*cJ,0;
            -Vnt(1)*sL+Vnt(2)*cL,0,0];
     for i = 1 : 8
         [xI(1),xI(2),xI(3)] = TransN2Ecef(Pos(1),Pos(2),Pos(3));
         rI = sqrt((xI(1)-matG1_G24(count,i*3-2))^2+(xI(2)-matG1_G24(count,i*3-1))^2+(xI(3)-matG1_G24(count,i*3-0))^2);
         E  = [xI(1)-matG1_G24(count,i*3-2),xI(2)-matG1_G24(count,i*3-1),xI(3)-matG1_G24(count,i*3-0)]/rI;
         H1 = [zeros(1,6) E*D zeros(1,6) -1 0];
         drI = E(1)*(dxI(1)-matV1_V24(count,i*3-2))+E(2)*(dxI(2)-matV1_V24(count,i*3-1))+E(3)*(dxI(3)-matV1_V24(count,i*3-0));
         M   = [dxI(1)-matV1_V24(count,i*3-2)-drI*E(1),dxI(2)-matV1_V24(count,i*3-1)-drI*E(2),dxI(3)-matV1_V24(count,i*3-0)-drI*E(3)]/rI;
         H2 = [zeros(1,3) E*Cne E*F+M*D zeros(1,6) 0 -1];
         H2 = [E*Cne zeros(1,3) E*F+M*D zeros(1,6) 0 -1];
         if(1)
             Z_K = rI-matPI(count,i)- matTu(count,1)-xc(1);
             resZK((count-1)/4+1,1) = Z_K;
             H_K  = H1; 
             R_K  = 9;
             K_K = P_K*H_K'*inv(H_K*P_K*H_K'+R_K);
             P_K = (eye(17)-K_K*H_K)*P_K;
             Zn = Z_K-H_K*X_K;
             X_K = X_K + K_K*Zn;
         end
         if(1)
             Z_K = drI-matPIs(count,i)- matTru(count,1);
             resZK((count-1)/4+1,2) = Z_K;
             H_K  = H2; 
             R_K  = 0.09;
             K_K = P_K*H_K'*inv(H_K*P_K*H_K'+R_K);
             P_K = (eye(17)-K_K*H_K)*P_K;
             Zn = Z_K-H_K*X_K;
             X_K = X_K + K_K*Zn;
         end 
     end
     %--×ËÌ¬ÐÞÕý-- Dr.Lu
%      phi_e   = X_K(4:6,1);
%      qnne = Rv2Quat(phi_e);
%      Cnb = quat_mul(qnne,Trans_att2quat(Trans_attm2att(Cnb)));
%      Cnb = Trans_quat2attm(Cnb);
     
     %--×ËÌ¬ÐÞÕý--Hui
     Cnb=(eye(3)+antisym_mat(X_K(4:6)))*(Cnb);
%      attm=(eye(3)+antisym_mat(X_K(4:6)))*Trans_quat2attm(Qnb);
%      att=Trans_attm2att(attm);
%      Qnb=Trans_att2quat(att); 

     %--ËÙ¶ÈÐÞÕý
     Vn = Vn - X_K(1:3,1);
     
     %--Î»ÖÃÐÞÕý
     Pos = Pos - X_K(7:9,1);
     
     %--Î±¾àÐÞÕý
     xc = xc -X_K(16:17,1);
     
     %--×´Ì¬ÇåÁã
     X_K(1:3,1) = X_K(1:3,1) - X_K(1:3,1);
     X_K(4:9,1) = X_K(4:9,1) - X_K(4:9,1);
     X_K(16,1) = zeros(1,1);
     
end
%--×ËÌ¬ÐÞÕý-- Dr.Mei
% function qof = Rv2Quat(Rv)
% %%
%     NormRv = norm(Rv);
%     if NormRv==0
%         f = 0.5;
%     else
%         f = sin(NormRv/2)/NormRv;
%     end
%     qof = [cos(NormRv/2);f*Rv(1);f*Rv(2);f*Rv(3)];
% end