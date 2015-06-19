%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function grid = readGridFromBinaryFile(file_name)
%
% Arguments
% file_name - file name 
% 
% Returns
% grid      - Grid structure
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

function grid = readGridFromBinaryFile(file_name)

fid = fopen(file_name,'r');

grid.num_dims = fread(fid,1,'int');
grid.x_lo = fread(fid,3,'double');
grid.x_hi = fread(fid,3,'double');
grid.x_lo_ghostbox = fread(fid,3,'double');
grid.x_hi_ghostbox = fread(fid,3,'double');
grid.grid_dims = fread(fid,3,'int');
grid.grid_dims_ghostbox = fread(fid,3,'int');
grid.dx = fread(fid,3,'double');
grid.num_gridpts = fread(fid,1,'int');

fclose(fid);
