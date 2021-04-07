function attm = Trans_att2attm(att)
    attm(1,1)=cos(att(2))*cos(att(3))+sin(att(2))*sin(att(3))*sin(att(1));
    attm(1,2)=sin(att(3))*cos(att(1));
    attm(1,3)=sin(att(2))*cos(att(3))-cos(att(2))*sin(att(3))*sin(att(1));
    attm(2,1)=-cos(att(2))*sin(att(3))+sin(att(2))*cos(att(3))*sin(att(1));
    attm(2,2)=cos(att(3))*cos(att(1));
    attm(2,3)=-sin(att(2))*sin(att(3))-cos(att(2))*cos(att(3))*sin(att(1));
    attm(3,1)=-sin(att(2))*cos(att(1));
    attm(3,2)=sin(att(1));
    attm(3,3)=cos(att(2))*cos(att(1));
end