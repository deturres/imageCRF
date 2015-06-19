%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File:        const_adv_1d_hat.m
% Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
%                  Regents of the University of Texas.  All rights reserved.
%              (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:    $Revision: 149 $
% Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
% Description: MATLAB demo program for 1D ENO/WENO spatial derivatives 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script demos the 1D ENO/WENO derivative calculations using 
% the advection equation:
%
%    u_t + v u_x = 0
%
% where v is a constant advection velocity.  The initial condition
% is a "hat"-function:
%
%   u(t=0) = 0      for -1 <= x < -0.5
%          = 0.5+x  for -0.5 <= x < 0
%          = 0.5-x  for 0 <= x < 0.5
%          = 0      for 0.5 <= x < 1
%
% and the boundary conditions are periodic.
%
% In this code, time advection is done using forward euler (TVD RK1), 
% TVD RK2, and TVD RK3.
%
% Kevin Chu
% Dept of Mathematics, MIT
% March 2005
%

% setup environment
clear
format long

% set physical parameters
v = 0.33;

% set up spatial grid parameters
N = 500;
max_ENO_order = 3;
ghostcell_width = max_ENO_order;
N_with_ghostcells = N+2*ghostcell_width;
x_lo = -1;
x_hi = 1;
dx = (x_hi-x_lo)/N;
x = (x_lo-(ghostcell_width-0.5)*dx:dx:x_hi+ghostcell_width*dx)';

% set advection velocity function
vel = v*ones(N_with_ghostcells,1);

% set up time integration parameters
cfl_RK1 = 0.1;
cfl_RK2 = 0.6;
cfl_RK3 = 1;
dt_RK1 = cfl_RK1*dx/abs(v);
dt_RK2 = cfl_RK2*dx/abs(v);
dt_RK3 = cfl_RK3*dx/abs(v);
t_i = 0;
t_f = 5;
t_RK1 = t_i:dt_RK1:t_f;
if (t_RK1(end) ~= t_f)
  t_RK1 = [t_RK1 t_f];
end
t_RK2 = t_i:dt_RK2:t_f;
if (t_RK2(end) ~= t_f)
  t_RK2 = [t_RK2 t_f];
end
t_RK3 = t_i:dt_RK3:t_f;
if (t_RK3(end) ~= t_f)
  t_RK3 = [t_RK3 t_f];
end

% initialize u
u_ENO1_1d_RK1 = zeros(N_with_ghostcells,1);
idx = find(x < -0.5);
u_ENO1_1d_RK1(idx) = 0;
idx = find(x >= -0.5 & x < 0);
u_ENO1_1d_RK1(idx) = 0.5+x(idx);
idx = find(x >= 0 & x < 0.5);
u_ENO1_1d_RK1(idx) = 0.5-x(idx);
idx = find(x >= 0.5);
u_ENO1_1d_RK1(idx) = 0;
u_ENO2_1d_RK1 = u_ENO1_1d_RK1;
u_ENO3_1d_RK1 = u_ENO1_1d_RK1;
u_WENO5_1d_RK1 = u_ENO1_1d_RK1;
u_ENO1_1d_RK2 = u_ENO1_1d_RK1;
u_ENO2_1d_RK2 = u_ENO1_1d_RK1;
u_ENO3_1d_RK2 = u_ENO1_1d_RK1;
u_WENO5_1d_RK2 = u_ENO1_1d_RK1;
u_ENO1_1d_RK3 = u_ENO1_1d_RK1;
u_ENO2_1d_RK3 = u_ENO1_1d_RK1;
u_ENO3_1d_RK3 = u_ENO1_1d_RK1;
u_WENO5_1d_RK3 = u_ENO1_1d_RK1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for forward euler (TVD RK1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK1)

  % fill boundary cells
  u_ENO1_1d_RK1(1:ghostcell_width) = u_ENO1_1d_RK1(N+1:ghostcell_width+N);
  u_ENO1_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK1(1:ghostcell_width) = u_ENO2_1d_RK1(N+1:ghostcell_width+N);
  u_ENO2_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK1(1:ghostcell_width) = u_ENO3_1d_RK1(N+1:ghostcell_width+N);
  u_ENO3_1d_RK1(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK1(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK1(1:ghostcell_width) = u_WENO5_1d_RK1(N+1:ghostcell_width+N);
  u_WENO5_1d_RK1(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK1(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK1,vel,ghostcell_width,dx);

  % advance solution
  u_ENO1_1d_RK1 = u_ENO1_1d_RK1 - dt_RK1*vel.*u_x_ENO1_1d;
  u_ENO2_1d_RK1 = u_ENO2_1d_RK1 - dt_RK1*vel.*u_x_ENO2_1d;
  u_ENO3_1d_RK1 = u_ENO3_1d_RK1 - dt_RK1*vel.*u_x_ENO3_1d;
  u_WENO5_1d_RK1 = u_WENO5_1d_RK1 - dt_RK1*vel.*u_x_WENO5_1d;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK2)

  % fill boundary cells
  u_ENO1_1d_RK2(1:ghostcell_width) = u_ENO1_1d_RK2(N+1:ghostcell_width+N);
  u_ENO1_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK2(1:ghostcell_width) = u_ENO2_1d_RK2(N+1:ghostcell_width+N);
  u_ENO2_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK2(1:ghostcell_width) = u_ENO3_1d_RK2(N+1:ghostcell_width+N);
  u_ENO3_1d_RK2(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK2(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK2(1:ghostcell_width) = u_WENO5_1d_RK2(N+1:ghostcell_width+N);
  u_WENO5_1d_RK2(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK2(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x for first stage of TVD RK2
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK2,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK2,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK2,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK2,vel,ghostcell_width,dx);

  % advance first stage of TVD RK2
  u_ENO1_1d_tmp1 = u_ENO1_1d_RK2 - dt_RK2*vel.*u_x_ENO1_1d;
  u_ENO2_1d_tmp1 = u_ENO2_1d_RK2 - dt_RK2*vel.*u_x_ENO2_1d;
  u_ENO3_1d_tmp1 = u_ENO3_1d_RK2 - dt_RK2*vel.*u_x_ENO3_1d;
  u_WENO5_1d_tmp1 = u_WENO5_1d_RK2 - dt_RK2*vel.*u_x_WENO5_1d;

  % compute approximations to u_x for second stage of TVD RK2
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK2
  u_ENO1_1d_RK2 = 0.5*(u_ENO1_1d_RK2 + u_ENO1_1d_tmp1 - dt_RK2*vel.*u_x_ENO1_1d);
  u_ENO2_1d_RK2 = 0.5*(u_ENO2_1d_RK2 + u_ENO2_1d_tmp1 - dt_RK2*vel.*u_x_ENO2_1d);
  u_ENO3_1d_RK2 = 0.5*(u_ENO3_1d_RK2 + u_ENO3_1d_tmp1 - dt_RK2*vel.*u_x_ENO3_1d);
  u_WENO5_1d_RK2 = 0.5*(u_WENO5_1d_RK2 + u_WENO5_1d_tmp1 - dt_RK2*vel.*u_x_WENO5_1d);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main time integration loop for TVD RK3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(t_RK3)

  % fill boundary cells
  u_ENO1_1d_RK3(1:ghostcell_width) = u_ENO1_1d_RK3(N+1:ghostcell_width+N);
  u_ENO1_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO1_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_ENO2_1d_RK3(1:ghostcell_width) = u_ENO2_1d_RK3(N+1:ghostcell_width+N);
  u_ENO2_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO2_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_ENO3_1d_RK3(1:ghostcell_width) = u_ENO3_1d_RK3(N+1:ghostcell_width+N);
  u_ENO3_1d_RK3(N+ghostcell_width+1:end) = ...
    u_ENO3_1d_RK3(ghostcell_width+1:2*ghostcell_width);
  u_WENO5_1d_RK3(1:ghostcell_width) = u_WENO5_1d_RK3(N+1:ghostcell_width+N);
  u_WENO5_1d_RK3(N+ghostcell_width+1:end) = ...
    u_WENO5_1d_RK3(ghostcell_width+1:2*ghostcell_width);

  % compute approximations to u_x
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_RK3,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_RK3,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_RK3,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_RK3,vel,ghostcell_width,dx);

  % advance first stage of TVD RK3
  u_ENO1_1d_tmp1 = u_ENO1_1d_RK3 - dt_RK3*vel.*u_x_ENO1_1d;
  u_ENO2_1d_tmp1 = u_ENO2_1d_RK3 - dt_RK3*vel.*u_x_ENO2_1d;
  u_ENO3_1d_tmp1 = u_ENO3_1d_RK3 - dt_RK3*vel.*u_x_ENO3_1d;
  u_WENO5_1d_tmp1 = u_WENO5_1d_RK3 - dt_RK3*vel.*u_x_WENO5_1d;

  % compute approximations to u_x for second stage of TVD RK3
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp1,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp1,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp1,vel,ghostcell_width,dx);

  % advance second stage of TVD RK3
  u_ENO1_1d_tmp2 = 0.75*u_ENO1_1d_RK3 + 0.25*(u_ENO1_1d_tmp1 - dt_RK3*vel.*u_x_ENO1_1d);
  u_ENO2_1d_tmp2 = 0.75*u_ENO2_1d_RK3 + 0.25*(u_ENO2_1d_tmp1 - dt_RK3*vel.*u_x_ENO2_1d);
  u_ENO3_1d_tmp2 = 0.75*u_ENO3_1d_RK3 + 0.25*(u_ENO3_1d_tmp1 - dt_RK3*vel.*u_x_ENO3_1d);
  u_WENO5_1d_tmp2 = 0.75*u_WENO5_1d_RK3 + 0.25*(u_WENO5_1d_tmp1 - dt_RK3*vel.*u_x_WENO5_1d);

  % compute approximations to u_x for third stage of TVD RK3
  u_x_ENO1_1d = UPWIND_HJ_ENO1_1D(u_ENO1_1d_tmp2,vel,ghostcell_width,dx);
  u_x_ENO2_1d = UPWIND_HJ_ENO2_1D(u_ENO2_1d_tmp2,vel,ghostcell_width,dx);
  u_x_ENO3_1d = UPWIND_HJ_ENO3_1D(u_ENO3_1d_tmp2,vel,ghostcell_width,dx);
  u_x_WENO5_1d = UPWIND_HJ_WENO5_1D(u_WENO5_1d_tmp2,vel,ghostcell_width,dx);

  % advance third stage of TVD RK3
  u_ENO1_1d_RK3 = 1/3*u_ENO1_1d_RK3 + 2/3*(u_ENO1_1d_tmp2 - dt_RK3*vel.*u_x_ENO1_1d);
  u_ENO2_1d_RK3 = 1/3*u_ENO2_1d_RK3 + 2/3*(u_ENO2_1d_tmp2 - dt_RK3*vel.*u_x_ENO2_1d);
  u_ENO3_1d_RK3 = 1/3*u_ENO3_1d_RK3 + 2/3*(u_ENO3_1d_tmp2 - dt_RK3*vel.*u_x_ENO3_1d);
  u_WENO5_1d_RK3 = 1/3*u_WENO5_1d_RK3 + 2/3*(u_WENO5_1d_tmp2 - dt_RK3*vel.*u_x_WENO5_1d);

end

% plot results
figure(1); clf;
plot(x,u_ENO1_1d_RK1,'b-.');
hold on;
plot(x,u_ENO2_1d_RK1,'g-.');
plot(x,u_ENO3_1d_RK1,'r-.');
plot(x,u_WENO5_1d_RK1,'m-.');

plot(x,u_ENO1_1d_RK2,'b--');
plot(x,u_ENO2_1d_RK2,'g--');
plot(x,u_ENO3_1d_RK2,'r--');
plot(x,u_WENO5_1d_RK2,'m--');

plot(x,u_ENO1_1d_RK3,'b');
plot(x,u_ENO2_1d_RK3,'g');
plot(x,u_ENO3_1d_RK3,'r');
plot(x,u_WENO5_1d_RK3,'m');
axis([-1 1 -0.1 0.6]);

