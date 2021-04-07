function qo = quat_conj(qi) %求取共轭四元数
    qo = [qi(1); -qi(2:4)];