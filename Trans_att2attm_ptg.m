function [ Tnb ] = Trans_att2attm_ptg( psi, theta, gamma )
%��ƫ���Ǧף������ǦȺͺ���Ǧü��㷽�����Ҿ���Tnb

S_psi = sin(psi);
S_theta = sin(theta);
S_gamma = sin(gamma);
C_psi = cos(psi);
C_theta = cos(theta);
C_gamma = cos(gamma);

Tnb = [C_gamma*C_psi - S_gamma*S_theta*S_psi C_gamma*S_psi + S_gamma*S_theta*C_psi -S_gamma*C_theta;...
    -C_theta*S_psi C_theta*C_psi S_theta;...
    S_gamma*C_psi + C_gamma*S_theta*S_psi S_gamma*S_psi - C_gamma*S_theta*C_psi C_gamma*C_theta]';
end