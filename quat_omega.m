function Q = quat_omega(wnb, Qini)
Q = 0.5*[0 -wnb(1) -wnb(2) -wnb(3);
         wnb(1) 0 wnb(3) -wnb(2);
         wnb(2) -wnb(3) 0 wnb(1);
         wnb(3) wnb(2) -wnb(1) 0]*Qini;
     return;