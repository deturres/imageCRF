c***********************************************************************
c
c  File:        lsm_field_extension2d.f
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
c  lsm2dComputeFieldExtensionEqnRHS() computes right-hand side of the 
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
      subroutine lsm2dComputeFieldExtensionEqnRHS(
     &  rhs,
     &  ilo_rhs_gb, ihi_rhs_gb,
     &  jlo_rhs_gb, jhi_rhs_gb,
     &  S,
     &  ilo_S_gb, ihi_S_gb,
     &  jlo_S_gb, jhi_S_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  S_x_upwind, S_y_upwind,
     &  ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb,
     &  jlo_grad_S_upwind_gb, jhi_grad_S_upwind_gb,
     &  signed_normal_x, signed_normal_y,
     &  ilo_signed_normal_gb, ihi_signed_normal_gb,
     &  jlo_signed_normal_gb, jhi_signed_normal_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_rhs_gb, ihi_rhs_gb
      integer jlo_rhs_gb, jhi_rhs_gb
      integer ilo_S_gb, ihi_S_gb
      integer jlo_S_gb, jhi_S_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_S_upwind_gb, ihi_grad_S_upwind_gb
      integer jlo_grad_S_upwind_gb, jhi_grad_S_upwind_gb
      integer ilo_signed_normal_gb, ihi_signed_normal_gb
      integer jlo_signed_normal_gb, jhi_signed_normal_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real rhs(ilo_rhs_gb:ihi_rhs_gb,
     &         jlo_rhs_gb:jhi_rhs_gb)
      real S(ilo_S_gb:ihi_S_gb,
     &       jlo_S_gb:jhi_S_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real S_x_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb,
     &                jlo_grad_S_upwind_gb:jhi_grad_S_upwind_gb)
      real S_y_upwind(ilo_grad_S_upwind_gb:ihi_grad_S_upwind_gb,
     &                jlo_grad_S_upwind_gb:jhi_grad_S_upwind_gb)
      real signed_normal_x(ilo_signed_normal_gb:ihi_signed_normal_gb,
     &                     jlo_signed_normal_gb:jhi_signed_normal_gb)
      real signed_normal_y(ilo_signed_normal_gb:ihi_signed_normal_gb,
     &                     jlo_signed_normal_gb:jhi_signed_normal_gb)
      real dx, dy
      integer i,j
      real zero
      parameter (zero=0.0d0)
      real zero_level_set_cutoff

c     set zero_level_set_cutoff to 3*max(dx,dy)
      zero_level_set_cutoff = 3.0d0*max(dx,dy)

c     compute RHS
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          if ( abs(phi(i,j)) .gt. zero_level_set_cutoff ) then

            rhs(i,j) = -( signed_normal_x(i,j)*S_x_upwind(i,j)
     &                  + signed_normal_y(i,j)*S_y_upwind(i,j) )

          else

            if ( (phi(i,j)*phi(i-1,j) .gt. zero) .and. 
     &           (phi(i,j)*phi(i+1,j) .gt. zero) .and.
     &           (phi(i,j)*phi(i,j-1) .gt. zero) .and.
     &           (phi(i,j)*phi(i,j+1) .gt. zero) ) then

              rhs(i,j) = -( signed_normal_x(i,j)*S_x_upwind(i,j)
     &                    + signed_normal_y(i,j)*S_y_upwind(i,j) )

            else

              rhs(i,j) = zero

            endif
          endif

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
