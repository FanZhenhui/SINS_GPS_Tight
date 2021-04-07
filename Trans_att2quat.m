function quat = Trans_att2quat(att)
    attm=Trans_att2attm(att);    

    U=[1+attm(1,1)+attm(2,2)+attm(3,3);
        1+attm(1,1)-attm(2,2)-attm(3,3);
        1-attm(1,1)+attm(2,2)-attm(3,3);
        1-attm(1,1)-attm(2,2)+attm(3,3)];
    if U(1)==max(U)
        Q=[U(1);attm(3,2)-attm(2,3);attm(1,3)-attm(3,1);attm(2,1)-attm(1,2)];
    end
    if U(2)==max(U)
            Q=[attm(3,2)-attm(2,3);U(2);attm(1,2)+attm(2,1);attm(1,3)+attm(3,1)];
    end
        if U(3)==max(U)
                Q=[attm(1,3)-attm(3,1);attm(2,1)+attm(1,2);U(3);attm(3,2)+attm(2,3)];
        end
            if U(4)==max(U)
                    Q=[attm(2,1)-attm(1,2);attm(1,3)+attm(3,1);attm(3,2)+attm(2,3);U(4)];
            end
    quat=Q/norm(Q);
end