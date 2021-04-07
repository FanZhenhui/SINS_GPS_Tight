function att = Trans_quat2att(quat)
    attm = Trans_quat2attm(quat);
    att = Trans_attm2att(attm);
end