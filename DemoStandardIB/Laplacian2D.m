function f = Laplacian2D(phi,N,h)

Nx=N(1);
Ny=N(2);

xp=[2:Nx,1];
xm=[Nx,1:Nx-1];

yp=[2:Ny,1];
ym=[Ny,1:Ny-1];


f = (phi(xp,:)+phi(xm,:)+phi(:,yp)+phi(:,ym)-4*phi)/(h*h);