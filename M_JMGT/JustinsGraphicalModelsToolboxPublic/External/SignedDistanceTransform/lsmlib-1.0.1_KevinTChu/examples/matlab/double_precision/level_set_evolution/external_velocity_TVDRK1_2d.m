%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        external_velocity_TVDRK1_2d.m
% Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
%                  Regents of the University of Texas.  All rights reserved.
%              (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:    $Revision: 149 $
% Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
% Description: MATLAB demo program for level set evolution functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script demos level set evolution functions when an external 
% velocity field drives the motion of the zero level set.
%
% The initial condition is the signed distance from the ellipse:
%
%    (x-0.25).^2 + 0.1*(y-0.25).^2 = 0.01
%
% The velocity is a constant vector (V_x, V_y).  The boundary conditions 
% are periodic in both coordinate directions.
%
% In this code, time advection is done using TVD RK1 (forward euler).
%
% NOTES:
% - All data arrays are in the ordered generated by the MATLAB meshgrid()
%   function.  That is, the data corresponding to the point (x_i,y_j) 
%   is stored at index (j,i).  This affects the dimensions of the velocity
%   arrays and how the ghostcell values are set.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup environment
clear
format long

% set up physical parameters for velocity and width of initial conditions
V_x = 0; V_y = 0.33;
V_x = 0.33; V_y = 0;
V_x = 0.33; V_y = -0.48;

% set up spatial grid parameters
Nx = 100;
Ny = 100;
spatial_derivative_order = 5;
ghostcell_width = 3;
Nx_with_ghostcells = Nx+2*ghostcell_width;
Ny_with_ghostcells = Ny+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
y_lo = -1;
y_hi = 1;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;
dX = [dx dy];
X = (x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx)';
Y = (y_lo-(ghostcell_width-0.5)*dy:dy:y_hi+ghostcell_width*dy)';
[x,y] = meshgrid(X,Y); 

% set advection velocity function
v_x = V_x*ones(Ny,Nx);
v_y = V_y*ones(Ny,Nx);
velocity{1} = v_x;
velocity{2} = v_y;

% set up time integration parameters
cfl_number = 0.4;
dt = cfl_number/(abs(V_x)/dx+abs(V_y)/dy);
t_start = 0;
t_end = 1;

% initialize phi
phi_init = (x-0.25).^2 + 0.1*(y-0.25).^2 - 0.01;
phi = phi_init;

% initialize time
t_cur = t_start;

% initialize flag to indicate if color axis has been set
color_axis_set = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK1 (forward euler)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (t_cur < t_end)

  % fill boundary cells 
  phi(1:ghostcell_width,:) = ...
    phi(Ny+1:ghostcell_width+Ny,:);
  phi(Ny+ghostcell_width+1:end,:) = ...
    phi(ghostcell_width+1:2*ghostcell_width,:);
  phi(:,1:ghostcell_width) = ...
    phi(:,Nx+1:ghostcell_width+Nx);
  phi(:,Nx+ghostcell_width+1:end) = ...
    phi(:,ghostcell_width+1:2*ghostcell_width);

  % advance level set function
  phi = advancePhiTVDRK1(phi, velocity, ghostcell_width, ...
                         dX, dt, cfl_number, spatial_derivative_order);

  % update current time
  t_cur = t_cur + dt
  
  % plot the current level set function and zero level set
  figure(1); clf;
  pcolor(x,y,phi);
  shading interp
  hold on 
  contour(x,y,phi,[0 0],'c','linewidth',2);
  xlabel('x');
  ylabel('y');
  axis([-1 1 -1 1]);
  axis square;
  if (color_axis_set == 0)
    color_axis = caxis;  % save color axis for initial data
  else
    caxis(color_axis);   % set color axis equal to the one for initial data
  end
  drawnow;

end

% plot results
figure(1); clf;
pcolor(x,y,phi);
shading interp
hold on
contour(x,y,phi,[0 0],'c','linewidth',2);
contour(x,y,phi_init,[0 0],'r','linewidth',2);
xlabel('x');
ylabel('y');
axis([-1 1 -1 1]);
axis square;
