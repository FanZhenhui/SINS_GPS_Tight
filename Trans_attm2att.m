function att = Trans_attm2att(attm)
    glv.R2D=180/pi;         
    att = [ asin(attm(3,2)); 
            atan2(-attm(3,1),attm(3,3)); 
            atan2(attm(1,2),attm(2,2))];
    if att(3)<0
        att(3)=att(3)+2*pi;
    end

    if att(3)>2*pi
        att(3)=att(3)-2*pi;
    end 
end
