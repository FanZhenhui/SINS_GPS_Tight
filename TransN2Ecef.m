function [x_u,y_u,z_u] =  TransN2Ecef(psi,lambda,h)
a = 6378137.0;
e = 1/298.257;
x_u = a*cos(lambda)/sqrt(1+(1-e^2)*tan(psi)*tan(psi)) + h*cos(lambda)*cos(psi);
y_u = a*sin(lambda)/sqrt(1+(1-e^2)*tan(psi)*tan(psi)) + h*sin(lambda)*cos(psi);
z_u = a*(1-e^2)*sin(psi)/sqrt(1-e^2*sin(psi)*sin(psi)) + h*sin(psi);
end