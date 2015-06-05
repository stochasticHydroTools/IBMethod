function adv = advection2D(u,N,h)

Nx=N(1);
Ny=N(2);

xp=[2:Nx,1];
xm=[Nx,1:Nx-1];

yp=[2:Ny,1];
ym=[Ny,1:Ny-1];

u1=u(:,:,1);
u2=u(:,:,2);

% 1st component of advection
Dxu1=(u1(xp,:)-u1(xm,:))/(2*h);
Dyu1=(u1(:,yp)-u1(:,ym))/(2*h);
u2i = (u2(xm,:)+u2)/2;
u2i = (u2i+u2i(:,yp))/2;
u1sq=u1.^2; Dxu1sq=(u1sq(xp,:)-u1sq(xm,:))/(2*h);
u1u2=u1.*u2i; Dyu1u2=(u1u2(:,yp)-u1u2(:,ym))/(2*h);
adv(:,:,1)=0.5*(u1.*Dxu1+u2i.*Dyu1)+0.5*(Dxu1sq+Dyu1u2);

% 2nd component of advection
u1i = (u1(:,ym)+u1)/2;
u1i = (u1i(xp,:)+u1i)/2;
Dxu2= (u2(xp,:)-u2(xm,:))/(2*h);
Dyu2= (u2(:,yp)-u2(:,ym))/(2*h);
u1u2= u1i.*u2; Dxu1u2=(u1u2(xp,:)-u1u2(xm,:))/(2*h);
u2sq= u2 .*u2; Dyu2sq=(u2sq(:,yp)-u2sq(:,ym))/(2*h);
adv(:,:,2)=0.5*(u1i.*Dxu2 + u2.*Dyu2) + 0.5*(Dxu1u2 + Dyu2sq);