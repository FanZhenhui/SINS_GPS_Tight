function vo = quat_mulVec(q, vi) %×ø±ê±ä»»
    qi = [0;vi];
    qo = quat_mul(quat_mul(q,qi),quat_conj(q));
    vo = qo(2:4,1);