function F=Force(X,Ns,ds)
% stiff Hooke's spring between markers

kappa = 1;
kp = [2:Ns,1];
km = [Ns,1:Ns-1];

F = kappa * (X(kp,:)+X(km,:)-2*X)/(ds*ds);