c***********************************************************************
c
c  File:        lsm_field_extension1d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision: 149 $
c  Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
c  Description: 2D F77 routines for extending fields off of the 
c               zero level set 
c
c***********************************************************************

c***********************************************************************
c
c  lsm1dComputeFieldExtensionEqnRHS() computes right-hand side of the 
c  field extension equation when it is written in the form:
c
c  S_t = -sgn(phi) N dot grad(S)
c
c  Arguments:
c    rhs (out):             right-hand side of field extension equation
c    S (in):                field to be extended off of the zero level set
c    phi (in):              level set function used to compute normal vector
c    S_*_upwind (in):       upwind spatial derivatives for grad(S)
c    signed_normal_* (in):  signed normal 
c    *_gb (in):             index range for ghostbox
c    *_fb (in):             index range for fillbox
c
c***********************************************************************
      subroutine lsm1dComputeFieldExtensionEqnRHS(
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  S,
     &  ilo_S_gb, ihi_S_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  S_x_upwind, 
     &  ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb,
     &  signed_normal_x,
     &  ilo_signed_normal_gb, ihi_signed_normal_gb,
     &  ilo_fb, ihi_fb,
     &  dx)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_rhs_gb, ihi_rhs_gb
      integer ilo_S_gb, ihi_S_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb
      integer ilo_signed_normal_gb, ihi_signed_normal_gb
      integer ilo_fb, ihi_fb
      real rhs(ilo_rhs_gb:ihi_rhs_gb)
      real S(ilo_S_gb:ihi_S_gb)
      real phi(ilo_phi_gb:ihi_phi_gb)
      real S_x_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb)
      real signed_normal_x(ilo_signed_normal_gb:ihi_signed_normal_gb)
      real dx
      integer i
      real zero
      parameter (zero=0.0d0)
      real zero_level_set_cutoff

c     set zero_level_set_cutoff to 3*dx
      zero_level_set_cutoff = 3.0d0*dx

c     compute RHS
c     { begin loop over grid
      do i=ilo_fb,ihi_fb

        if ( abs(phi(i)) .gt. zero_level_set_cutoff ) then

          rhs(i) = -signed_normal_x(i)*S_x_upwind(i)

        else

          if ( (phi(i+1)*phi(i) .le. zero) .and. 
     &         (phi(i-1)*phi(i) .le. zero) ) then

            rhs(i) = -signed_normal_x(i)*S_x_upwind(i)

          else

            rhs(i) = zero

          endif
        endif

      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
