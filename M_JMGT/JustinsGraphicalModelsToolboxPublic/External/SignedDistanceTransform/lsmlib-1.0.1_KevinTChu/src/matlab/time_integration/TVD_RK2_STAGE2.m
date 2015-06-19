%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TVD_RK2_STAGE2 advances the solution through the final stage of a 
%                second-order TVD Runge-Kutta step.
%
% Usage: u_next = TVD_RK2_STAGE2(u_stage1, u_cur, rhs, dt)
%
% Arguments:
% - u_stage1:    u_approx(t_cur+dt)
% - u_cur:       u(t_cur)
% - rhs:         right-hand side of time evolution equation at (t_cur+dt)
% - dt:          step size
%
% Return value:
% - u_next:      u(t_cur + dt)
%
% NOTES:
% - u_stage1, u_cur and rhs are assumed to have the same dimensions.
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

function u_next = TVD_RK2_STAGE2(u_stage1, u_cur, rhs, dt)

u_next = 0.5*(u_cur + u_stage1 + dt*rhs);
