% 
%   This code constructs a signed distance function from
%   the interface defined by {phi0 = 0}.
%
clear all
close all
N=100;         % number of grid points in one direction
numc=2;        % number of circles for the initial condition
R1=.5;         % initial radius of circles 
h=2/(N-1);     % grid spacing
dt=.2*h;       % time step
tfin=.6;       % total simulation time
nit=tfin/dt;   % number of time steps
vn=1;           % normal speed this must be positive
x=-1:h:1;
y=x;
[X,Y]=meshgrid(x);
%
%
%   Initialize the level set function
a=.2+2.*abs(atan2(X,Y))/(2.*pi);
phi0=a.*((1.5*(X.*X)+(Y.*Y)).^(1/2)-.75);
phi=phi0;
figure(1)
c=-1.5:.05:1.5;
contour(X,Y,phi,c);
hold on
contour(X,Y,phi,[0,0],'k');
hold off
axis([-1 1 -1 1])
axis('square')
pause
%
%      arrays for the periodic boundary conditions
       for i=1:N
         ip(i)=i+1;
         im(i)=i-1;
       end
       im(1)=N;
       ip(N)=1;
%
%      begin simulation loop
       for iter=1:nit
           for i=1:N
             for j=1:N
              dmx=(phi(i,j)-phi(im(i),j))/h;                     % x backward difference
              dpx=(phi(ip(i),j)-phi(i,j))/h;                     % x forward difference
              dmy=(phi(i,j)-phi(i,im(j)))/h;                     % y backward difference
              dpy=(phi(i,ip(j))-phi(i,j))/h;                     % y forward difference
              if(phi0(i,j) > 0)
              fluxx=max(abs(max(dmx,0)),abs(min(dpx,0)));        % Godunov Flux x direction
              fluxy=max(abs(max(dmy,0)),abs(min(dpy,0)));        % Godunov Flux y direction
              else
              fluxx=max(abs(max(dpx,0)),abs(min(dmx,0)));        % Godunov Flux x direction
              fluxy=max(abs(max(dpy,0)),abs(min(dmy,0)));        % Godunov Flux y direction
              end
              sgn=phi0(i,j)/(phi0(i,j)^2+h^2)^(1/2);
              grp=(fluxx^2+fluxy^2)^(1/2);         
              phin(i,j)=phi(i,j)-(grp-1)*sgn*dt;                   % advance by dt
            end
           end 
           phi=phin;                                             % update
%
%         Plotting
          c=-1.5:.05:1.5;
          contour(X,Y,phi,c);
          hold on
          contour(X,Y,phi,[0,0],'k');
          hold off
          axis([-1 1 -1 1])
          axis('square')
          pause(.001)
       end