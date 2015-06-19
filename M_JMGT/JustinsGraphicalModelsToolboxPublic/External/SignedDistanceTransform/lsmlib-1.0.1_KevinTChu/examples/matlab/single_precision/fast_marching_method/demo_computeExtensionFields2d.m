%
% File:        demo_computeExtensionFields2d.m
% Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
%                  Regents of the University of Texas.  All rights reserved.
%              (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:    $Revision: 149 $
% Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
% Description: MATLAB demo code for computeExtensionFields2d MEX file
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script demos the computeExtensionFields2d MATLAB MEX-function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% setup environment
clear
format long

% plotting parameters
num_plots_per_example = 5;

% grid parameters
Nx = 200;
Ny = 200;
x_lo = -1;
x_hi = 1;
y_lo = -1;
y_hi = 1;
dx = (x_hi-x_lo)/Nx;
dy = (y_hi-y_lo)/Ny;
x = (x_lo:dx:x_hi)';
y = (y_lo:dy:y_hi)';
dX = [dx dy];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example Problem 1 
% -----------------
%  * interface: two circles centered at (0.25,0.25) and (-0.25,-0.25)
%               both with radius 0.2
%  * velocity:  uniform velocity.  U = (1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Example Problem 1...');
demo_problem_num = 1;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.25]; radius1 = 0.2;
center2 = [-0.25 -0.25]; radius2 = 0.2;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = single(min(dist1, dist2));
vel_x = single(ones(size(X)));
vel_y = single(ones(size(X)));
source_fields = cell(2,1);
source_fields{1} = vel_x;
source_fields{2} = vel_y;
[distance_function, ext_fields] = ...
  computeExtensionFields2d(phi, source_fields, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((demo_problem_num-1)*num_plots_per_example+1); clf;
contourf(X,Y,distance_function,10);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+2); clf;
contourf(X,Y,phi,10);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example Problem 2
% -----------------
%  * interface: two circles centered at (0.25,0.25) and (-0.25,-0.25)
%               both with radius 0.2
%  * velocity:  "expanding velocity".  U = (x,y)/r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Example Problem 2...');
demo_problem_num = 2;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.25]; radius1 = 0.2;
center2 = [-0.25 -0.25]; radius2 = 0.2;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = single(min(dist1, dist2));
r = sqrt(X.^2 + Y.^2);
idx = find(r ~= 0);
vel_x = single(zeros(size(X))); vel_y = single(zeros(size(Y)));
vel_x(idx) = X(idx)./r(idx);
vel_y(idx) = Y(idx)./r(idx);
idx = find(r == 0);
vel_x(idx) = 0;
vel_y(idx) = 0;
source_fields = cell(2,1);
source_fields{1} = vel_x;
source_fields{2} = vel_y;
[distance_function, ext_fields] = ...
  computeExtensionFields2d(phi, source_fields, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((demo_problem_num-1)*num_plots_per_example+1); clf;
contourf(X,Y,distance_function,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+2); clf;
contourf(X,Y,phi,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+3); clf;
contourf(X,Y,ext_fields{1},50);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+4); clf;
contourf(X,Y,source_fields{1},20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example Problem 3
% -----------------
%  * interface: two circles centered at (0.25,0.1) and (-0.35,0.25)
%               both with radii of 0.2 and 0.3 respectively
%  * velocity:  "expanding velocity".  U = (x,y)/r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Example Problem 3...');
demo_problem_num = 3;
[X,Y] = meshgrid(x,y);  
center1 = [0.25 0.1]; radius1 = 0.2;
center2 = [-0.35 0.25]; radius2 = 0.3;
dist1 = sqrt( (X-center1(1)).^2 + (Y-center1(2)).^2 ) - radius1;
dist2 = sqrt( (X-center2(1)).^2 + (Y-center2(2)).^2 ) - radius2;
phi = single(min(dist1, dist2));
r = sqrt(X.^2 + Y.^2);
idx = find(r ~= 0);
vel_x = single(zeros(size(X))); vel_y = single(zeros(size(Y)));
vel_x(idx) = X(idx)./r(idx);
vel_y(idx) = Y(idx)./r(idx);
idx = find(r == 0);
vel_x(idx) = 0;
vel_y(idx) = 0;
source_fields = cell(2,1);
source_fields{1} = vel_x;
source_fields{2} = vel_y;
[distance_function, ext_fields] = ...
  computeExtensionFields2d(phi, source_fields, dX);

% compare FMM distance function with exact distance function
err = max(max(abs(distance_function-phi)))

% plot results
figure((demo_problem_num-1)*num_plots_per_example+1); clf;
contourf(X,Y,distance_function,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+2); clf;
contourf(X,Y,phi,20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+3); clf;
contourf(X,Y,ext_fields{1},50);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

figure((demo_problem_num-1)*num_plots_per_example+4); clf;
contourf(X,Y,source_fields{1},20);
axis([x_lo x_hi y_lo y_hi]);
pbaspect([1 1 1]);
colorbar;

