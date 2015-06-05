% ----------------------------------------------------------------------- %
% NS2D_IBMAC.m
% Yuanxun Bill Bao
% Nov 09, 2013
% 
% Description: solve the fluid-structure interaction of an elastic
% membrane by solving a 2D NS flow. The standard immersed 
% boundary on a staggered grid is used (IB-MAC). 
% A second-order centered scheme is used for temporal integration.
%
% Inputs:
% L      -- domain length
% N      -- grid size N=[Nx,Ny], Nx=Ny for now
% nu     -- viscosity
% Kernel -- [KernelID, DerivativeID, K]
% X0     -- initial configuration of the memberane
% u0     -- intial velocity
% tend   -- final time
% dt     -- time step
% ----------------------------------------------------------------------- %

function NS2D_IBMAC(L,N,nu,Kernel,X,u,tend,dt,Nf,MpC,showplot,mexFlag)

Nx=N(1); Ny=N(2);
h=L/Nx;
rho=1;
Ns=length(X); 
ds=2*pi/Ns;

x=(0:Nx-1)*h;
y=(0:Ny-1)*h;
[yy,xx]=meshgrid(y,x);
xx2=xx+h/2; yy2=yy+h/2;
sk=1:1:Nx;

kp = [2:Nx,    1];
km = [Nx  ,1:Nx-1];
xp=[2:Nx,1];
xm=[Nx,1:Nx-1];
yp=[2:Ny,1];
ym=[Ny,1:Ny-1];

% set Kernel
kID=Kernel{1}; K=Kernel{3};

% 2D Laplacian in Fourier space
[m2,m1] = meshgrid(fftshift(-Ny/2:Ny/2-1), fftshift(-Nx/2:Nx/2-1));
L_hat1  = -4/(h*h) * (sin(pi/N(1)*m1).^2 + sin(pi/N(2)*m2).^2);
L_hat2  = L_hat1; L_hat2(1,1)=1;

% % add passive tracers
% Ntracer = 2048*4;
% s = (0:Ntracer-1)*(2*pi/Ntracer);
% alpha=1/4; beta=1/4;
% %alpha = 5/28; beta = 7/20;
% % Xtracer = [alpha*cos(s'), beta*sin(s')]*L+L/2;
% % Xtracer0 = Xtracer;
% % %areaTracer(1) = polyarea([Xtracer(:,1); Xtracer(1,1)], [Xtracer(:,2);Xtracer(1,2)]);
% % %[ampl,ppint,xf,yf] = makeIBcurve(Xtracer(:,1)-L/2,Xtracer(:,2)-L/2,100);
% % areaTracer(1) = ppint;

% time step
Nt = floor(tend/dt);

uIBMAC=zeros(Nx,Ny,2,Nf+1); XIBMAC=zeros(Ns,2,Nf+1);
areaIBMAC=zeros(1,Nt+1);
uIBMAC(:,:,:,1)=u;
XIBMAC(:,:,1)=X;
areaIBMAC(1) = polyarea([X(:,1); X(1,1)], [X(:,2);X(1,2)]);

% use 2nd-order RK to get one extra initial value
uold=u; Xold=X; % Xtracer_old = Xtracer;
% compute X(n+1/2)
XX=X+(dt/2)*interpMAC2Dvector(u,X,N,h,kID,K,mexFlag);
%XXtracer = Xtracer + (dt/2) * interpMAC2Dvector(u,Xtracer,N,h,kID,K,mexFlag);
% force-spreading f(n+1/2)
FF=Force(XX,Ns,ds)*ds;%*(1+2*0.45*sin(9*0.5*dt))*4;
ff=spreadMAC2Dvector(FF,XX,N,h,kID,K,mexFlag);
% solve for uu=u(n+1/2)
ADV=advection2D(u,N,h);
w1=(-2*rho/dt)*u(:,:,1)+rho*ADV(:,:,1)-ff(:,:,1);
w2=(-2*rho/dt)*u(:,:,2)+rho*ADV(:,:,2)-ff(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);
c1=dt/(2*rho); c2=nu*dt/(2*rho);
r1=u(:,:,1)-(dt/2)*ADV(:,:,1)+c1*ff(:,:,1)-c1*(pp-pp(xm,:))/h;
r2=u(:,:,2)-(dt/2)*ADV(:,:,2)+c1*ff(:,:,2)-c1*(pp-pp(:,ym))/h;
uu1_hat=fft2(r1)./(1-c2*L_hat1);
uu2_hat=fft2(r2)./(1-c2*L_hat1);
uu(:,:,1)=real(ifft2(uu1_hat));
uu(:,:,2)=real(ifft2(uu2_hat));

% solve for u(n+1)
ADVold = ADV;
ADV=advection2D(uu,N,h);
w1=(-rho/dt)*u(:,:,1)+rho*ADV(:,:,1)-nu/2*Laplacian2D(u(:,:,1),N,h)-ff(:,:,1);
w2=(-rho/dt)*u(:,:,2)+rho*ADV(:,:,2)-nu/2*Laplacian2D(u(:,:,2),N,h)-ff(:,:,2);
Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
pp=PoissonSolver2D(Divw,L_hat2);

c3=dt/rho;
r1=u(:,:,1)-dt*ADV(:,:,1)+c2*Laplacian2D(u(:,:,1),N,h)+c3*ff(:,:,1)-c3*(pp-pp(xm,:))/h;
r2=u(:,:,2)-dt*ADV(:,:,2)+c2*Laplacian2D(u(:,:,2),N,h)+c3*ff(:,:,2)-c3*(pp-pp(:,ym))/h;
u1_hat=fft2(r1)./(1-c2*L_hat1);
u2_hat=fft2(r2)./(1-c2*L_hat1);
u(:,:,1)=real(ifft2(u1_hat));
u(:,:,2)=real(ifft2(u2_hat));
ADV=advection2D(u,N,h);

% compute X(n+1)
X=X+dt*interpMAC2Dvector(uu,XX,N,h,kID,K,mexFlag);
% areaIBMAC(2) = polyarea([X(:,1); X(1,1)],[X(:,2); X(1,2)]);
% 
% Xtracer = Xtracer + dt * interpMAC2Dvector(uu,XXtracer,N,h,kID,K,mexFlag);
%[ampl,ppint,xf,yf] = makeIBcurve(Xtracer(:,1)-L/2,Xtracer(:,2)-L/2,100);
%areaTracer(2) = ppint;
%areaTracer(2) = polyarea([Xtracer(:,1); Xtracer(1,1)], [Xtracer(:,2);Xtracer(1,2)]);

% vol(1)=abs(areaIBMAC(2)-areaIBMAC(1))/areaIBMAC(1);
% voltracer(1) = abs(areaTracer(2)-areaTracer(1))/areaTracer(1);


for n=2:Nt
    
    % compute X(n+1) using 2nd order AB
    Uold=interpMAC2Dvector(uold,Xold,N,h,kID,K,mexFlag);
    U   =interpMAC2Dvector(u,X,N,h,kID,K,mexFlag);
    UU  = 3/2*U - 1/2*Uold;
    Xold = X;
    X=X+dt*UU;
    
%     Utracer_old = interpMAC2Dvector(uold,Xtracer_old,N,h,kID,K);
%     Utracer = interpMAC2Dvector(u,Xtracer,N,h,kID,K);
%     UUtracer = 3/2*Utracer - 1/2*Utracer_old;
%     Xtracer_old = Xtracer;
%     Xtracer = Xtracer + dt * UUtracer;
    %[ampl,ppint,xf,yf] = makeIBcurve(Xtracer(:,1)-L/2,Xtracer(:,2)-L/2,100);
    % areaTracer(n+1) = ppint;
    
    % compute f(n+1/2)
    XX=(X+Xold)/2;
    FF=Force(XX,Ns,ds)*ds;
    ff=spreadMAC2Dvector(FF,XX,N,h,kID,K,mexFlag);
    
    % fluid solve
    S=3/2*ADV-1/2*ADVold;
    w1=(-rho/dt)*u(:,:,1)+rho*S(:,:,1)-nu/2*Laplacian2D(u(:,:,1),N,h)-ff(:,:,1);
    w2=(-rho/dt)*u(:,:,2)+rho*S(:,:,2)-nu/2*Laplacian2D(u(:,:,2),N,h)-ff(:,:,2);
    Divw=(w1(xp,:)-w1)/h + (w2(:,yp)-w2)/h;
    pp=PoissonSolver2D(Divw,L_hat2);
    c3=dt/rho;
    r1=u(:,:,1)-dt*S(:,:,1)+c2*Laplacian2D(u(:,:,1),N,h)+c3*ff(:,:,1)-c3*(pp-pp(xm,:))/h;
    r2=u(:,:,2)-dt*S(:,:,2)+c2*Laplacian2D(u(:,:,2),N,h)+c3*ff(:,:,2)-c3*(pp-pp(:,ym))/h;
    u1_hat=fft2(r1)./(1-c2*L_hat1);
    u2_hat=fft2(r2)./(1-c2*L_hat1);
    uold=u;
    u(:,:,1)=real(ifft2(u1_hat));
    u(:,:,2)=real(ifft2(u2_hat));
    
    ADVold = ADV;
    ADV=advection2D(u,N,h);
        
%     areaIBMAC(n+1) = polyarea([X(:,1); X(1,1)], [X(:,2);X(1,2)]);
%     
%     %%% caculate fore on the marker
%     Fmarker=Force(X,Ns,ds)*ds;%*(1+2*0.45*sin(9*(n-0.5)*dt))*4;
%     for k = 1 : Ns
%         [th,rr] = cart2pol(X(k,1)-L/2,X(k,2)-L/2);
%         P = [cos(th), sin(th); -sin(th), cos(th)];
%         Fpol(k,:) = (P * Fmarker(k,:)')'./norm(Fmarker(k,:));
%         TH(k) = th;
%     end
    

    % plotting
    if mod(n,floor(Nt/Nf)) == 0
        disp([num2str(n/Nt*100),'%'])
        
        if strcmp(showplot, 'on')
            
%             % ADD passive tracers
%             Utracer=interpMAC2Dvector(u,Xtracer,N,h,kID,K);
%             Xtracer=Xtracer+dt*Utracer;
%             areaTracer(n+1) = polyarea([Xtracer(:,1); Xtracer(1,1)], [Xtracer(:,2);Xtracer(1,2)]);
        
            j=n/(Nt/Nf);
            XIBMAC(:,:,j+1)=X;
            uIBMAC(:,:,:,j+1)=u;
            uu = (u(kp,:,1) + u(:,:,1))/2;
            vv = (u(:,kp,2) + u(:,:,2))/2;
            uMax=max(max(sqrt(uu.^2+vv.^2)));
            Re=rho*uMax*(1/4)/nu;
            %vort = (u(:,:,2)-u(km,:,2) - u(:,:,1) + u(:,km,1) )/h;
            %contourf(xx,yy,vort ,100); hold on
            %shading flat;  
            %colormap(jet); colorbar; %caxis([-10,10])
            quiver(xx2(sk,sk),yy2(sk,sk),uu(sk,sk),vv(sk,sk),1); hold on;
%             plot(mod(Xtracer(:,1),L),mod(Xtracer(:,2),L),'r-','linewidth',2)
%             plot(mod(Xtracer0(:,1),L),mod(Xtracer0(:,2),L),'g-','linewidth',2)
            plot(mod(X(:,1),L),mod(X(:,2),L),'k-','linewidth',2);    
            %plot([X(1,1),X(Ns/2,1)],[X(1,2),X(Ns/2,2)],'ro','markerfacecolor','r');
            xlabel(['t = ',num2str(n*dt), ', Re = ',num2str(Re),', N = ',num2str(Nx),...
                    ', # marker per cell = ',num2str(MpC)])
            title(['NS, IB-MAC, ',kID])
            axis equal
            axis([0,L,0,L])
            %axis([.2,.5,.5,.8]);
            %axis([.25,.35,.6,.7])
            drawnow
            hold off


%             subplot(6,2,[9,11]);
%             magU= sqrt(uu.^2 + vv.^2);
%             contourf(xx2,yy2,magU);
%             shading flat;
%             colorbar;
%             axis equal;
%             axis([0,L,0,L])
%             title({'magnitude of velocity: |u|',['max |U|: ',num2str(uMax)]})
%             drawnow
% 
%     %         subplot(2,2,3); hold on;
%     %         plot(n*dt,norm(X(1,:)-X(Ns/2,:),2),'k.');
%     %         axis([0,tend,0,5])
%     %         ylabel('distance between red dots')
%     %         hold off;
% 
% 
%             subplot(6,2,[10,12]); 
%             vol(n)=abs(areaIBMAC(n+1)-areaIBMAC(1))/areaIBMAC(1);
%             voltracer(n) = abs(areaTracer(n+1)-areaTracer(1))/areaTracer(1);
%             semilogy((1:n)*dt,vol, 'k-','linewidth',1); hold on;
%             semilogy((1:n)*dt,voltracer, 'r-','linewidth',1); 
%             ylabel('area loss');
%             set(gca,'Ytick',10.^[-13:2:-1]);
%             axis([0,tend,1e-13,1e-1])
%             legend('IB points', 'Tracers')

            %print('-dpng',['./movies/IBMAC1/',num2str(j+10000),'.png'])
        
        end
    end    
end


