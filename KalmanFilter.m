function [Xk, Pk, Kk] = KalmanFilter(Fhikk_1, Qk, Xk_1, Pk_1, Hk, Rk, Zk)
    Xkk_1=Fhikk_1*Xk_1;    
    Pkk_1 = Fhikk_1*Pk_1*Fhikk_1' + Qk; 
    Pxz = Pkk_1*Hk';                  
    Pzz = Hk*Pxz + Rk;                 
    Kk = Pxz*Pzz^-1;                  
    Xk = Xkk_1 + Kk*(Zk-Hk*Xkk_1);
    Pk = Pkk_1-Kk*Hk*Pkk_1;   
end
 
