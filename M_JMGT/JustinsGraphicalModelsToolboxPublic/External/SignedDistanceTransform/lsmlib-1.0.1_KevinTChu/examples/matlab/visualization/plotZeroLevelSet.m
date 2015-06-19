%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function plotZeroLevelSet(data_filename, grid_filename)
%
% Arguments:
% data_filename - file name for data array
% grid_filename - file name for binary Grid array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision: 149 $
% Modified:   $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotZeroLevelSet(data_filename, grid_filename)

[data, nx, ny, nz] = readDataArray(data_filename);
grid = readGridFromBinaryFile(grid_filename);

if( nz == 1) % 2D visualization
    
  x = grid.x_lo_ghostbox(1)*ones(1,nx) + grid.dx(1)*(0:nx-1);
  y = grid.x_lo_ghostbox(2)*ones(1,ny) + grid.dx(2)*(0:ny-1);
  [X Y]= meshgrid(x,y);
  
  data  = data';
  
  figure;
  level = [ 0 0 ];
  [ garbage, hI ] = contour(X, Y, data, level,'b-'); xlabel('x'); ylabel('y');
  
else % 3D visualization
    
  x = grid.x_lo_ghostbox(1)*ones(1,nx) + grid.dx(1)*(0:nx-1);
  y = grid.x_lo_ghostbox(2)*ones(1,ny) + grid.dx(2)*(0:ny-1);
  z = grid.x_lo_ghostbox(3)*ones(1,nz) + grid.dx(3)*(0:nz-1);
  [Y X Z] = meshgrid(y,x,z);
  
  figure;
  level = 0;
  figure, h = patch(isosurface(Y,X,Z,data,level));
  set(h,'FaceColor',[0.9 0.9 0.9],'EdgeColor', 'none');
  daspect([1 1 1]); camlight;
  xlabel('x'); ylabel('y'); zlabel('z');
  view(3)
  
end    
