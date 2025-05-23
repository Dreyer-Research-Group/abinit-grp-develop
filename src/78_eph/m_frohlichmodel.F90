!!****m* ABINIT/m_frohlichmodel
!! NAME
!!  m_frohlichmodel
!!
!! FUNCTION
!!  Compute ZPR, temperature-dependent electronic structure, and other properties
!!  using the Frohlich model
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2022 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_frohlichmodel

 use defs_basis
 use m_abicore
 use m_errors
 use m_crystal
 use m_ebands
 use m_efmas_defs
 use m_ifc
 use m_dtset

 use m_fstrings,            only : sjoin, itoa
 use m_gaussian_quadrature, only : cgqf

 implicit none

 private 

 public :: hamiltonian, frohlichmodel, polaronmass

contains
!!***

!!****f* m_frohlichmodel/frohlichmodel
!! NAME
!!  frohlichmodel
!!
!! FUNCTION
!! Main routine to compute properties based on the Frohlich model
!!
!! INPUTS
!! cryst<crystal_t>=Structure defining the unit cell
!! dtset<dataset_type>=All input variables for this dataset.
!! efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band degeneracy at each k point
!! efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!   efmasval(:,:)%eig2_diag band curvature double tensor
!! ifc<ifc_type>=contains the dynamical matrix and the IFCs.
!!
!! PARENTS
!!      m_eph_driver
!!
!! CHILDREN
!!      cgqf,ifc%calcnwrite_nana_terms,zheev
!!
!! SOURCE

subroutine frohlichmodel(cryst, dtset, efmasdeg, efmasval, ifc)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ifc_type),intent(in) :: ifc
!arrays
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables ------------------------------
!scalars
 logical :: sign_warn
 integer :: deg_dim,iband,ideg,idir,ikpt,imode,info,ipar,iphi,iqdir,itheta
 integer :: jband,lwork,nphi,nqdir,ntheta
 real(dp) :: angle_phi,cosph,costh,sinph,sinth,weight,weight_phi
 real(dp) :: zpr_frohlich,zpr_q0_avg,zpr_q0_fact
 !character(len=500) :: msg
!arrays
 logical, allocatable :: saddle_warn(:), start_eigf3d_pos(:)
 logical :: lutt_found(3), lutt_warn(3)
 real(dp) :: kpt(3), lutt_params(3), lutt_unit_kdir(3,3)
 real(dp), allocatable :: eigenval(:), rwork(:), unit_qdir(:,:)
 real(dp), allocatable :: lutt_dij(:,:), lutt_eigenval(:,:)
 real(dp), allocatable :: m_avg(:), m_avg_frohlich(:)
 real(dp), allocatable :: gq_points_th(:),gq_weights_th(:)
 real(dp), allocatable :: gq_points_cosph(:),gq_points_sinph(:)
 real(dp), allocatable :: weight_qdir(:)
 real(dp), allocatable :: polarity_qdir(:,:,:)
 real(dp), allocatable :: proj_polarity_qdir(:,:)
 real(dp), allocatable :: zpr_q0_phononfactor_qdir(:)
 real(dp), allocatable :: frohlich_phononfactor_qdir(:)
 real(dp), allocatable :: phfrq_qdir(:,:)
 real(dp), allocatable :: dielt_qdir(:)
 real(dp), allocatable :: zpr_frohlich_avg(:)
 complex(dpc), allocatable :: eigenvec(:,:), work(:)
 complex(dpc), allocatable :: eig2_diag_cart(:,:,:,:)
 complex(dpc), allocatable :: f3d(:,:)

!************************************************************************

 !!! Initialization of integrals
 ntheta   = dtset%efmas_ntheta
 nphi     = 2*ntheta
 nqdir     = nphi*ntheta

 ABI_MALLOC(gq_points_th,(ntheta))
 ABI_MALLOC(gq_weights_th,(ntheta))
 ABI_MALLOC(gq_points_cosph,(nphi))
 ABI_MALLOC(gq_points_sinph,(nphi))

 ABI_MALLOC(unit_qdir,(3,nqdir))
 ABI_MALLOC(weight_qdir,(nqdir))

 call cgqf(ntheta,1,zero,zero,zero,pi,gq_points_th,gq_weights_th)
 weight_phi=two*pi/real(nphi,dp)
 do iphi=1,nphi
   angle_phi=weight_phi*(iphi-1)
   gq_points_cosph(iphi)=cos(angle_phi)
   gq_points_sinph(iphi)=sin(angle_phi)
 enddo
 nqdir=0
 do itheta=1,ntheta
   costh=cos(gq_points_th(itheta))
   sinth=sin(gq_points_th(itheta))
   weight=gq_weights_th(itheta)*weight_phi*sinth
   do iphi=1,nphi
     cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)
     nqdir=nqdir+1

     unit_qdir(1,nqdir)=sinth*cosph
     unit_qdir(2,nqdir)=sinth*sinph
     unit_qdir(3,nqdir)=costh
     weight_qdir(nqdir)=weight

   enddo
 enddo

 ABI_FREE(gq_points_th)
 ABI_FREE(gq_weights_th)
 ABI_FREE(gq_points_cosph)
 ABI_FREE(gq_points_sinph)

 ABI_MALLOC(polarity_qdir,(3,3*cryst%natom,nqdir))
 ABI_MALLOC(proj_polarity_qdir,(3*cryst%natom,nqdir))
 ABI_MALLOC(zpr_q0_phononfactor_qdir,(nqdir))
 ABI_MALLOC(frohlich_phononfactor_qdir,(nqdir))
 ABI_MALLOC(phfrq_qdir,(3*cryst%natom,nqdir))
 ABI_MALLOC(dielt_qdir,(nqdir))

 !Compute phonon frequencies and mode-polarity for each qdir
 call ifc%calcnwrite_nana_terms(cryst, nqdir, unit_qdir, phfrq2l=phfrq_qdir, polarity2l=polarity_qdir)

 !Compute dielectric tensor for each qdir
 do iqdir=1,nqdir
   dielt_qdir(iqdir)=DOT_PRODUCT(unit_qdir(:,iqdir),MATMUL(ifc%dielt(:,:),unit_qdir(:,iqdir)))
 enddo

 !Compute projection of mode-polarity on qdir, and other derived quantities summed over phonon branches for each iqdir.
 !Note that acoustic modes are discarded (imode sum starts only from 4)
 zpr_q0_phononfactor_qdir=zero
 zpr_q0_avg=zero
 frohlich_phononfactor_qdir=zero
 do iqdir=1,nqdir
   do imode=4,3*cryst%natom
     proj_polarity_qdir(imode,iqdir)=DOT_PRODUCT(unit_qdir(:,iqdir),polarity_qdir(:,imode,iqdir))
     zpr_q0_phononfactor_qdir(iqdir)=zpr_q0_phononfactor_qdir(iqdir)+&
&      proj_polarity_qdir(imode,iqdir)**2 / phfrq_qdir(imode,iqdir) **2
     frohlich_phononfactor_qdir(iqdir)=frohlich_phononfactor_qdir(iqdir)+&
&      proj_polarity_qdir(imode,iqdir)**2 / phfrq_qdir(imode,iqdir) **(three*half)
   enddo
   zpr_q0_avg=zpr_q0_avg+&
&    weight_qdir(iqdir)*zpr_q0_phononfactor_qdir(iqdir)/dielt_qdir(iqdir)**2
 enddo
 zpr_q0_avg=zpr_q0_avg*quarter*piinv
 zpr_q0_fact=zpr_q0_avg*eight*pi*(three*quarter*piinv)**third*cryst%ucvol**(-four*third)

!DEBUG
! do iqdir=1,nqdir,513
!   write(std_out,'(a,3f8.4,3es12.4)')' unit_qdir,dielt_qdir,zpr_q0_phononfactor_qdir,frohlich_phononfactor=',&
!&    unit_qdir(:,iqdir),dielt_qdir(iqdir),zpr_q0_phononfactor_qdir(iqdir),frohlich_phononfactor_qdir(iqdir)
!   do imode=1,3*cryst%natom
!     write(std_out,'(a,i5,6es12.4)')'   imode,phfrq_qdir,phfrq(cmm1),polarity_qdir=',&
!&     imode,phfrq_qdir(imode,iqdir),phfrq_qdir(imode,iqdir)*Ha_cmm1,polarity_qdir(:,imode,iqdir),proj_polarity_qdir(imode,iqdir)
!   enddo
! enddo
! write(std_out,'(2a,3es12.4)')ch10,&
!& ' zpr_q0_avg, zpr_q0_fact, zpr_q0_fact (eV) =',zpr_q0_avg, zpr_q0_fact, zpr_q0_fact*Ha_eV
!ENDDEBUG

 write(ab_out,'(6a,f14.6,a,f14.6,a)') ch10,&
&  ' Rough correction to the ZPR, to take into account the missing q=0 piece using Frohlich model:',ch10,&
&  ' (+ for occupied states, - for unoccupied states) * zpr_q0_fact / (Nqpt_full_bz)**(1/3) ',ch10,&
&  ' where Nqpt_full_bz=number of q wavevectors in full BZ, and zpr_q0_fact=',zpr_q0_fact,' Ha=',zpr_q0_fact*Ha_eV,' eV'

 !Compute effective masses, and integrate the Frohlich model
 do ikpt=1,dtset%nkpt

   kpt(:)=dtset%kptns(:,ikpt)
   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)

     deg_dim    = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1

     ABI_MALLOC(eig2_diag_cart,(3,3,deg_dim,deg_dim))

     !Convert eig2_diag to cartesian coordinates
     do iband=1,deg_dim
        do jband=1,deg_dim
          eig2_diag_cart(:,:,iband,jband)=efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)
          eig2_diag_cart(:,:,iband,jband)=&
&           matmul(matmul(cryst%rprimd,eig2_diag_cart(:,:,iband,jband)),transpose(cryst%rprimd))/two_pi**2
        enddo
     enddo

     ABI_MALLOC(f3d,(deg_dim,deg_dim))
     ABI_MALLOC(m_avg,(deg_dim))
     ABI_MALLOC(m_avg_frohlich,(deg_dim))
     ABI_MALLOC(zpr_frohlich_avg,(deg_dim))
     ABI_MALLOC(eigenval,(deg_dim))
     ABI_MALLOC(saddle_warn,(deg_dim))
     ABI_MALLOC(start_eigf3d_pos,(deg_dim))

     m_avg=zero
     m_avg_frohlich=zero
     saddle_warn=.false.

     !Initializations for the diagonalization routine
     if(deg_dim>1)then

       ABI_MALLOC(eigenvec,(deg_dim,deg_dim))
       lwork=-1
       ABI_MALLOC(rwork,(3*deg_dim-2))
       ABI_MALLOC(work,(1))
       call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
       lwork=int(work(1))
       ABI_FREE(work)
       ABI_MALLOC(work,(lwork))

     endif

     !Compute the Luttinger parameters for the cubic case (deg_dim=3)
     if(deg_dim==3) then

       ABI_MALLOC(lutt_eigenval, (3,deg_dim))
       ABI_MALLOC(lutt_dij, (deg_dim, deg_dim))

       !Define unit_kdir for Luttinger parameters
       lutt_unit_kdir(:,1) = (/1,0,0/)
       lutt_unit_kdir(:,2) = 1/sqrt(2.0)*(/1,1,0/)
       lutt_unit_kdir(:,3) = 1/sqrt(3.0)*(/1,1,1/)

       !Degeneracy problems warning
       lutt_warn=(/.false.,.false.,.false./)

       !Inverse effective mass tensor eigenvalues in lutt_unit_kdir directions
       do idir=1,3
         do iband=1,deg_dim
           do jband=1,deg_dim
             lutt_dij(iband,jband)=&
&             DOT_PRODUCT(lutt_unit_kdir(:,idir),MATMUL(eig2_diag_cart(:,:,iband,jband),lutt_unit_kdir(:,idir)))
           enddo
         enddo

         eigenvec=lutt_dij ; lutt_eigenval(idir,:)=zero
         work=zero     ; rwork=zero
         call zheev('V','U',deg_dim,eigenvec,deg_dim,lutt_eigenval(idir,:),work,lwork,rwork,info)
         ABI_CHECK(info == 0, sjoin("zheev returned info:", itoa(info)))
       enddo

       !Check degeneracies in (100) direction, and evaluate A and B.
       !Eigenvalues are 2*A (d=1), 2*B (d=2)
       if(abs(lutt_eigenval(1,2)-lutt_eigenval(1,3))<tol5) then
         lutt_params(2)=0.5*((lutt_eigenval(1,2)+lutt_eigenval(1,3))/2)
         lutt_params(1)=0.5*lutt_eigenval(1,1)
       else if(abs(lutt_eigenval(1,2)-lutt_eigenval(1,1))<tol5) then
         lutt_params(2)=0.5*((lutt_eigenval(1,2)+lutt_eigenval(1,1))/2)
         lutt_params(1)=0.5*lutt_eigenval(1,3)
       else
         lutt_warn(1)=.true.
       endif

       !Check degeneracies in (111) direction and evaluate C
       !Eigenvalues are 2/3*(A+2B-C) (d=2), 2/3*(A+2B+2C) (d=1)
       if(abs(lutt_eigenval(3,2)-lutt_eigenval(3,3))<tol5) then
         lutt_params(3)=lutt_params(1)+2*lutt_params(2)-1.5*(0.5*(lutt_eigenval(3,2)+lutt_eigenval(3,3)))
       else if(abs(lutt_eigenval(3,2)-lutt_eigenval(3,1))<tol5) then
         lutt_params(3)=lutt_params(1)+2*lutt_params(2)-1.5*(0.5*(lutt_eigenval(3,2)+lutt_eigenval(3,1)))
       else
         lutt_warn(2)=.true.
       endif

       !Verify that the (110) direction eigenvalues are coherent with Luttinger parameters
       !Eigenvalues are 2B, A+B-C, A+B+C
       lutt_found=(/.false.,.false.,.false./)
       do ipar=1,deg_dim
         if(abs(lutt_eigenval(2,ipar)-2*lutt_params(2))<tol4) then
           lutt_found(1)=.true.
         else if(abs(lutt_eigenval(2,ipar)-(lutt_params(1)+lutt_params(2)-lutt_params(3)))<tol4) then
           lutt_found(2)=.true.
         else if(abs(lutt_eigenval(2,ipar)-(lutt_params(1)+lutt_params(2)+lutt_params(3)))<tol4) then
           lutt_found(3)=.true.
         endif
       enddo

       if(.not. (all(lutt_found))) then
         lutt_warn(3)=.true.
       endif

       ABI_FREE(lutt_eigenval)
       ABI_FREE(lutt_dij)

     endif !Luttinger parameters

     !Perform the integral over the sphere
     zpr_frohlich_avg=zero
     do iqdir=1,nqdir
       do iband=1,deg_dim
         do jband=1,deg_dim
           f3d(iband,jband)=DOT_PRODUCT(unit_qdir(:,iqdir),MATMUL(eig2_diag_cart(:,:,iband,jband),unit_qdir(:,iqdir)))
         enddo
       enddo

       if(deg_dim==1)then
         eigenval(1)=f3d(1,1)
       else
         eigenvec = f3d ; eigenval = zero
         work=zero      ; rwork=zero
         call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
         ABI_CHECK(info == 0, sjoin("zheev returned info:", itoa(info)))
       endif

       m_avg = m_avg + weight_qdir(iqdir)*eigenval
       m_avg_frohlich = m_avg_frohlich + weight_qdir(iqdir)/(abs(eigenval))**half
       zpr_frohlich_avg = zpr_frohlich_avg + &
&        weight_qdir(iqdir) * frohlich_phononfactor_qdir(iqdir)/((abs(eigenval))**half *dielt_qdir(iqdir)**2)

       if(iqdir==1) start_eigf3d_pos = (eigenval > 0)

       do iband=1,deg_dim
         if(start_eigf3d_pos(iband) .neqv. (eigenval(iband)>0)) then
           saddle_warn(iband)=.true.
         end if
       end do

     enddo

     if(deg_dim>1)then
       ABI_FREE(eigenvec)
       ABI_FREE(rwork)
       ABI_FREE(work)
     endif

     m_avg = quarter*piinv*m_avg
     m_avg = one/m_avg

     m_avg_frohlich = quarter*piinv * m_avg_frohlich
     m_avg_frohlich = m_avg_frohlich**2

     zpr_frohlich_avg = quarter*piinv * zpr_frohlich_avg

     if(deg_dim==1)then
       write(ab_out,'(2a,3(f6.3,a),i5)')ch10,&
&        ' - At k-point (',kpt(1),',',kpt(2),',',kpt(3),'), band ',&
&        efmasdeg(ikpt)%degs_bounds(1,ideg)
     else
       write(ab_out,'(2a,3(f6.3,a),i5,a,i5)')ch10,&
&        ' - At k-point (',kpt(1),',',kpt(2),',',kpt(3),'), bands ',&
&        efmasdeg(ikpt)%degs_bounds(1,ideg),' through ',efmasdeg(ikpt)%degs_bounds(2,ideg)
     endif

     !Print the Luttinger for the cubic case (deg_dim=3)
     if(deg_dim==3) then
       if (.not. (any(saddle_warn))) then
         if(any(lutt_warn)) then
           ! Warn for degeneracy breaking in inverse effective mass tensor eigenvalues
           write(ab_out, '(2a)') ch10, ' Luttinger parameters could not be determined:'
           if (lutt_warn(1)) then
             write(ab_out, '(a)') '     Predicted degeneracies for deg_dim = 3 are not met for (100) direction.'
           endif
           if (lutt_warn(2)) then
             write(ab_out, '(a)') '     Predicted degeneracies for deg_dim = 3 are not met for (111) direction.'
           endif
           if (lutt_warn(3)) then
             write(ab_out, '(a)') '     Predicted inverse effective mass tensor eigenvalues for direction (110) are not met.'
           endif
           write(ab_out, '(a)') ch10
         else
           write(ab_out, '(a,3f14.6)') ' Luttinger parameters (A, B, C) [at. units]: ',lutt_params(:)
         endif
       endif
     endif

     sign_warn=.false.
     do iband=1,deg_dim
       if(saddle_warn(iband)) then
         write(ab_out,'(a,i5,a)') ' Band ',efmasdeg(ikpt)%degs_bounds(1,ideg)+iband-1,&
&          ' SADDLE POINT - Frohlich effective mass and ZPR cannot be defined. '
         sign_warn=.true.
       else
         m_avg_frohlich(iband) = DSIGN(m_avg_frohlich(iband),m_avg(iband))
         zpr_frohlich_avg(iband) = -DSIGN(zpr_frohlich_avg(iband),m_avg(iband))
         write(ab_out,'(a,i5,a,f14.10)') &
&          ' Band ',efmasdeg(ikpt)%degs_bounds(1,ideg)+iband-1,&
&          ' Angular average effective mass for Frohlich model (<m**0.5>)**2= ',m_avg_frohlich(iband)
       endif
       if(start_eigf3d_pos(iband) .neqv. start_eigf3d_pos(1))then
         sign_warn=.true.
       endif
     enddo

     if(sign_warn .eqv. .false.)then
       zpr_frohlich = four*pi* two**(-half) * (sum(zpr_frohlich_avg(1:deg_dim))/deg_dim) / cryst%ucvol
       write(ab_out,'(2a)')&
&       ' Angular and band average effective mass and ZPR for Frohlich model.'
       write(ab_out,'(a,es16.6)') &
&       ' Value of     (<<m**0.5>>)**2 = ',(sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim)**2
       write(ab_out,'(a,es16.6)') &
&       ' Absolute Value of <<m**0.5>> = ', sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim
       write(ab_out,'(a,es16.6,a,es16.6,a)') &
&       ' ZPR from Frohlich model      = ',zpr_frohlich,' Ha=',zpr_frohlich*Ha_eV,' eV'
     else
       write(ab_out,'(a)')&
&        ' Angular and band average effective mass for Frohlich model cannot be defined because of a sign problem.'
     endif

     ABI_FREE(eig2_diag_cart)
     ABI_FREE(f3d)
     ABI_FREE(m_avg)
     ABI_FREE(m_avg_frohlich)
     ABI_FREE(zpr_frohlich_avg)
     ABI_FREE(eigenval)
     ABI_FREE(saddle_warn)
     ABI_FREE(start_eigf3d_pos)

   enddo ! ideg
 enddo ! ikpt

 ABI_FREE(unit_qdir)
 ABI_FREE(weight_qdir)
 ABI_FREE(polarity_qdir)
 ABI_FREE(proj_polarity_qdir)
 ABI_FREE(phfrq_qdir)
 ABI_FREE(dielt_qdir)
 ABI_FREE(zpr_q0_phononfactor_qdir)
 ABI_FREE(frohlich_phononfactor_qdir)

 end subroutine frohlichmodel
 
!!***

!!****f* m_frohlichmodel/polaronmass
!! NAME
!!  polaronmass
!!
!! FUNCTION
!! Improved routine to compute properties based on the Frohlich model, including effective masses in the cubic case.
!!
!! INPUTS
!! cryst<crystal_t>=Structure defining the unit cell
!! dtset<dataset_type>=All input variables for this dataset.
!! efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band degeneracy at each k point
!! efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!   efmasval(:,:)%eig2_diag band curvature double tensor
!! ifc<ifc_type>=contains the dynamical matrix and the IFCs.
!!
!! PARENTS
!!      m_eph_driver
!!
!! CHILDREN
!!      cgqf,ifc%calcnwrite_nana_terms,zheev
!!
!! SOURCE

subroutine polaronmass(cryst, dtset, efmasdeg, efmasval, ifc)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ifc_type),intent(in) :: ifc
!arrays
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables ------------------------------
!scalars
 integer  :: deg_dim, signpm
 integer  :: i, iband, jband, ideg, idir, iqdir, ieig
 integer  :: ikpt, ixi, ipar, iphi, iphon, imode, itheta, ik
 integer  :: nkdir, nqdir, ntheta, nphi, nxi, nkgrid
 integer  :: info, lwork
 real(dp) :: angle_phi,cosph,costh,sinph,sinth,weight,weight_phi
 real(dp) :: qpt, krange, nq_factor
 !character(len=500) :: msg
!arrays
!Electronic
 real(dp), allocatable :: eigenvec(:,:,:,:), eigenval(:,:,:)
 real(dp), allocatable :: rwork(:), unit_qdir(:,:)
 complex(dpc), allocatable :: leigenvec(:,:), work(:)
 complex(dpc), allocatable :: eig2_diag_cart(:,:,:,:)
!Luttinger
 logical  :: lutt_found(3), lutt_warn(3)
 real(dp) :: lutt_params(3)
 real(dp) :: kpoint(3)
 real(dp), allocatable :: lutt_dij(:,:), lutt_eigenval(:,:), leigenval(:)
!Dielectric
 real(dp), allocatable :: dielt_qdir(:)
!Phonons
 real(dp), allocatable :: gq_points_th(:),gq_weights_th(:)
 real(dp), allocatable :: gq_points_cosph(:),gq_points_sinph(:)
 real(dp), allocatable :: weight_qdir(:)
 real(dp), allocatable :: polarity_qdir(:,:,:)
 real(dp), allocatable :: phfrq_qdir(:,:)
!Self-energy and polaron mass
 real(dp) :: xi
 real(dp) :: temporary1(3,3), temporary2(3), temporary3
 real(dp) :: unitary_33(3,3)
 real(dp) :: minelecmass,eham(3,3)
 real(dp) :: kpt(3), k_plus_q(3)
 real(dp), allocatable :: lutt_unit_kdir(:,:)
 real(dp), allocatable :: omega_zero(:)
 real(dp), allocatable :: intsum(:,:,:,:)
 real(dp), allocatable :: sigma(:,:,:), d2sigmadk2(:,:)
 real(dp), allocatable :: invepsilonstar(:)
 real(dp), allocatable :: invemass(:,:), invemass_ieig(:,:),invpolmass(:,:)

!************************************************************************

!Define Luttinger and Phonon integration parameters
!Based solely solely one value - ntheta
 ntheta   = dtset%efmas_ntheta
 nphi     = 2*ntheta
 nqdir    = nphi*ntheta

!Define constants
! 3x3 Unitary matrix
 unitary_33  = 0.0_dp
 do i=1,3
  unitary_33(i,i) = 1.0_dp
 enddo

!Define unit_kdir for Luttinger parameters
 nkdir=3
 ABI_MALLOC(lutt_unit_kdir,(3,nkdir))
 lutt_unit_kdir(:,1) = (/1,0,0/)
 lutt_unit_kdir(:,2) = 1/sqrt(2.0)*(/1,1,0/)
 lutt_unit_kdir(:,3) = 1/sqrt(3.0)*(/1,1,1/)
!These are for testing purpose only
! lutt_unit_kdir(:,4) = (/0,1,0/)
! lutt_unit_kdir(:,5) = (/0,0,1/)

!Compute effective masses, and integrate the Frohlich model
 do ikpt=1,dtset%nkpt

   kpt(:)=dtset%kptns(:,ikpt)
   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)

     deg_dim    = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1

     ABI_MALLOC(eig2_diag_cart,(3,3,deg_dim,deg_dim))

     !Convert eig2_diag to cartesian coordinates
     do iband=1,deg_dim
       do jband=1,deg_dim
         eig2_diag_cart(:,:,iband,jband)=efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)
         eig2_diag_cart(:,:,iband,jband)=&
&           matmul(matmul(cryst%rprimd,eig2_diag_cart(:,:,iband,jband)),transpose(cryst%rprimd))/two_pi**2
       enddo ! jband
    enddo ! iband

    ABI_MALLOC(leigenval,(deg_dim))

    !Initializations for the diagonalization routine
    if(deg_dim>1)then

      ABI_MALLOC(leigenvec,(deg_dim,deg_dim))
      lwork=-1
      ABI_MALLOC(rwork,(3*deg_dim-2))
      ABI_MALLOC(work,(1))
      call zheev('V','U',deg_dim,leigenvec,deg_dim,leigenval,work,lwork,rwork,info)
      lwork=int(work(1))
      ABI_FREE(work)
      ABI_MALLOC(work,(lwork))

    endif

    !Compute the Luttinger parameters for the cubic case (deg_dim=3)
     if(deg_dim==3) then

       ABI_MALLOC(lutt_eigenval, (3,deg_dim))
       ABI_MALLOC(lutt_dij, (deg_dim, deg_dim))

       !Degeneracy problems warning
       lutt_warn=(/.false.,.false.,.false./)

       !Inverse effective mass tensor eigenvalues in lutt_unit_kdir directions
       do idir=1,3
         do iband=1,deg_dim
           do jband=1,deg_dim
             lutt_dij(iband,jband)=&
&             DOT_PRODUCT(lutt_unit_kdir(:,idir),MATMUL(eig2_diag_cart(:,:,iband,jband),lutt_unit_kdir(:,idir)))
           enddo
         enddo

         leigenvec=lutt_dij ; lutt_eigenval(idir,:)=zero
         work=zero     ; rwork=zero
         call zheev('V','U',deg_dim,leigenvec,deg_dim,lutt_eigenval(idir,:),work,lwork,rwork,info)
         ABI_CHECK(info == 0, sjoin("zheev returned info:", itoa(info)))
       enddo

       !Check degeneracies in (100) direction, and evaluate A and B.
       !Eigenvalues are 2*A (d=1), 2*B (d=2)
       if(abs(lutt_eigenval(1,2)-lutt_eigenval(1,3))<tol5) then
         lutt_params(2)=0.5*((lutt_eigenval(1,2)+lutt_eigenval(1,3))/2)
         lutt_params(1)=0.5*lutt_eigenval(1,1)
       else if(abs(lutt_eigenval(1,2)-lutt_eigenval(1,1))<tol5) then
         lutt_params(2)=0.5*((lutt_eigenval(1,2)+lutt_eigenval(1,1))/2)
         lutt_params(1)=0.5*lutt_eigenval(1,3)
       else
         lutt_warn(1)=.true.
       endif

       !Check degeneracies in (111) direction and evaluate C
       !Eigenvalues are 2/3*(A+2B-C) (d=2), 2/3*(A+2B+2C) (d=1)
       if(abs(lutt_eigenval(3,2)-lutt_eigenval(3,3))<tol5) then
         lutt_params(3)=lutt_params(1)+2*lutt_params(2)-1.5*(0.5*(lutt_eigenval(3,2)+lutt_eigenval(3,3)))
       else if(abs(lutt_eigenval(3,2)-lutt_eigenval(3,1))<tol5) then
         lutt_params(3)=lutt_params(1)+2*lutt_params(2)-1.5*(0.5*(lutt_eigenval(3,2)+lutt_eigenval(3,1)))
       else
         lutt_warn(2)=.true.
       endif

       !Verify that the (110) direction eigenvalues are coherent with Luttinger parameters
       !Eigenvalues are 2B, A+B-C, A+B+C
       lutt_found=(/.false.,.false.,.false./)
       do ipar=1,deg_dim
         if(abs(lutt_eigenval(2,ipar)-2*lutt_params(2))<tol4) then
           lutt_found(1)=.true.
         else if(abs(lutt_eigenval(2,ipar)-(lutt_params(1)+lutt_params(2)-lutt_params(3)))<tol4) then
           lutt_found(2)=.true.
         else if(abs(lutt_eigenval(2,ipar)-(lutt_params(1)+lutt_params(2)+lutt_params(3)))<tol4) then
           lutt_found(3)=.true.
         endif
       enddo

       if(.not. (all(lutt_found))) then
         lutt_warn(3)=.true.
       endif

       ABI_FREE(lutt_eigenval)
       ABI_FREE(lutt_dij)

     endif !Luttinger parameters

     if(deg_dim>1)then
       ABI_FREE(leigenvec)
       ABI_FREE(rwork)
       ABI_FREE(work)
     endif

     !Print the Luttinger for the cubic case (deg_dim=3)
     if(deg_dim==3) then
         if(any(lutt_warn)) then
           ! Warn for degeneracy breaking in inverse effective mass tensor eigenvalues
           write(ab_out, '(2a)') ch10, ' Luttinger parameters could not be determined:'
           if (lutt_warn(1)) then
             write(ab_out, '(a)') '     Predicted degeneracies for deg_dim = 3 are not met for (100) direction.'
           endif
           if (lutt_warn(2)) then
             write(ab_out, '(a)') '     Predicted degeneracies for deg_dim = 3 are not met for (111) direction.'
           endif
           if (lutt_warn(3)) then
             write(ab_out, '(a)') '     Predicted inverse effective mass tensor eigenvalues for direction (110) are not met.'
           endif
           write(ab_out, '(a)') ch10
         else
           write(ab_out, '(a,3f14.6)') '   Luttinger parameters (A, B, C) (a.u.): ',lutt_params(:)
         endif
     endif

     ABI_FREE(eig2_diag_cart)
     ABI_FREE(leigenval)

   enddo ! ideg
 enddo ! ikpt

!Build inverse electronic effective mass for different directions from Luttinger params
 deg_dim = 3
 ABI_MALLOC(invemass,(deg_dim,nkdir))
 ABI_MALLOC(invemass_ieig,(deg_dim,nkdir))

! 100 -direction
 invemass(1,1) = two*lutt_params(1) ! 2A
 invemass(2,1) = two*lutt_params(2) ! 2B
 invemass(3,1) = two*lutt_params(2) ! 2B
! 110 -direction
 invemass(1,2) = lutt_params(1) + lutt_params(2) + lutt_params(3) ! A + B + C
 invemass(2,2) = MIN(two*lutt_params(2), lutt_params(1) + lutt_params(2) - lutt_params(3) ) ! 2B
 invemass(3,2) = MAX(two*lutt_params(2), lutt_params(1) + lutt_params(2) - lutt_params(3) ) ! A + B - C
! 111 -direction
 invemass(1,3) = two*( lutt_params(1) + two*lutt_params(2) + two*lutt_params(3) ) / three ! 2(A + 2B + 2C)/3
 invemass(2,3) = two*( lutt_params(1) + two*lutt_params(2) -     lutt_params(3) ) / three ! 2(A + 2B - 2C)/3
 invemass(3,3) = two*( lutt_params(1) + two*lutt_params(2) -     lutt_params(3) ) / three ! 2(A + 2B - 2C)/3

!END Diagonalize 3x3 Luttinger-Kohn Hamiltonian 

 ABI_MALLOC(unit_qdir,(3,nqdir))
 ABI_MALLOC(polarity_qdir,(3,3*cryst%natom,nqdir))
 ABI_MALLOC(phfrq_qdir,(3*cryst%natom,nqdir))
 ABI_MALLOC(dielt_qdir,(nqdir))

 ABI_MALLOC(gq_points_th,(ntheta))
 ABI_MALLOC(gq_weights_th,(ntheta))
 ABI_MALLOC(gq_points_cosph,(nphi))
 ABI_MALLOC(gq_points_sinph,(nphi))
 ABI_MALLOC(weight_qdir,(nqdir))

 call cgqf(ntheta,1,zero,zero,-one,one,gq_points_th,gq_weights_th)
 weight_phi=two*pi/real(nphi,dp)
 do iphi=1,nphi
   angle_phi=weight_phi*(iphi-1)
   gq_points_cosph(iphi)=cos(angle_phi)
   gq_points_sinph(iphi)=sin(angle_phi)
 enddo
 nqdir=0
 do itheta=1,ntheta
   costh=gq_points_th(itheta)
   sinth=sqrt(one-costh**2)
   weight=gq_weights_th(itheta)*weight_phi
   do iphi=1,nphi
     cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)
     nqdir=nqdir+1

     unit_qdir(1,nqdir)=sinth*cosph
     unit_qdir(2,nqdir)=sinth*sinph
     unit_qdir(3,nqdir)=costh
     weight_qdir(nqdir)=weight

   enddo
 enddo

 ABI_FREE(gq_points_th)
 ABI_FREE(gq_weights_th)
 ABI_FREE(gq_points_cosph)
 ABI_FREE(gq_points_sinph)

!In the following, compute invepsilonstar(imode) benefiting from the cubic symmetry.
!Retrieve IR active phonon frequencies
!Compute phonon frequencies and mode-polarity for each qdir (actually, only the first would suffice)
 call ifc%calcnwrite_nana_terms(cryst, nqdir, unit_qdir, phfrq2l=phfrq_qdir, polarity2l=polarity_qdir)

!Calculate inverse epsilon* for each optical phonon mode
! (epsilon*)**(-1) = 4Pi/Omega_0 (p_j0/(dielt_qdir*omeja_j0))**2
 ABI_MALLOC(invepsilonstar,(3*cryst%natom))
 ABI_MALLOC(omega_zero,(3*cryst%natom))
 invepsilonstar = zero

 do imode = 1,3*cryst%natom
   !For ease of treatment loop over all phonon branches but...
   !Avoid the acoustic branches
   if(imode > 3) then
     invepsilonstar(imode) = &
       four*pi/cryst%ucvol*(dot_product(unit_qdir(:,1),polarity_qdir(:,imode,1)) &
       /( ifc%dielt(1,1)*phfrq_qdir(imode,1)) )**two
   endif
 enddo ! imode

!Set nkgrid and krange
!For the time being we need 3 points for the finite difference
 nkgrid = 3
!Set material dependent length scale for the finite difference
!Lowest optical phonon frequency to be used
 minelecmass=1.0_dp/maxval(abs(invemass))
 krange = sqrt(two*minelecmass*phfrq_qdir(4,1)/1000.0)

!Diagonalize 3x3 Luttinger-Kohn Hamiltonian 

!Initializations for the diagonalization routine
 ABI_MALLOC(eigenval,(deg_dim,nkgrid,nkdir))
 ABI_MALLOC(eigenvec,(deg_dim,deg_dim,nkgrid,nkdir))
 lwork=-1
 ABI_MALLOC(rwork,(3*deg_dim-2))
 ABI_MALLOC(work,(1))

!Initialize eigenval
 eigenval = zero
 do idir = 1,nkdir
   do ik = 1, nkgrid
     kpoint(:) = (ik -1.0_dp)*krange*lutt_unit_kdir(:,idir)
     eham = hamiltonian(lutt_params, kpoint)
     call dsyev('V','U',3,eham,3,eigenval(1:3,ik,idir),work,lwork,info)
     eigenvec(:,:,ik,idir) = eham
     lwork=int(work(1))
     ABI_FREE(work)
     ABI_MALLOC(work,(lwork))
   enddo
 enddo

 ABI_FREE(rwork)
 ABI_FREE(work)

 ABI_MALLOC(intsum,(deg_dim,deg_dim,nkgrid,nkdir))
 ABI_MALLOC(sigma,(nkgrid,deg_dim,nkdir))

!main loop
!Dummy value for omega_zero
 omega_zero  = phfrq_qdir(:,1)
!Establish VB or CB sign 
 signpm = 1
 if(lutt_params(1) .LT. zero) then
   signpm = -1
 endif

!Define Luttinger and Phonon integration parameters
!One value input - effmass_ntheta input variable
 nxi      = ntheta

!Normalization factor for the integral
 nq_factor = nxi*two/pi*(four*pi)

!Initialize self-energy
 sigma = zero

!Summation over IR active phonon modes
 do iphon = 1, 3*cryst%natom
   if( invepsilonstar(iphon) > 1E-10 ) then
     !Summation over relevant directions (100,110,111)
     do idir=1,nkdir
       !Summation over electronic eigenvalues 
       do ieig=1,3
         !Summation over k-range 
         !Here consider a small BZ region around Gamma
         !A convergene study might be performed around this value
         do ik = 1,nkgrid
           !Define k-point vector around which one integrates
           !Needs at least 2 considering TRS and finite difference employed later
           !A three point formula is implemented later.
           kpoint(:) = ( ik - 1.0_dp) * krange*lutt_unit_kdir(:,idir)
           !Perform hyperbolic tangent integration for the semi-infinite q domain
           !Use a mapping to a tangent function - faster convergence wrt qpt sampling
           do ixi = 0,nxi
             xi = ixi*pi/( two * nxi)
             if ( ixi .EQ. nxi ) xi = xi - tol8
             !Wave-vector length
             qpt = ( omega_zero(iphon)/abs(lutt_params(1)) )**half*tan(xi)
!            XG 20211106 : the best integration scheme for a finite interval is based on Gauss-Legendre approach,
!            which was set up previously with nqdir point. So replaced the itheta and iphi loop by the loop over nqdir
             do iqdir=1,nqdir
               k_plus_q = kpoint + qpt*unit_qdir(:,iqdir)
               intsum(:,:,ik, idir) = &
                      abs(eigenval(ieig,ik,idir))*unitary_33 - &
                      ( signpm*hamiltonian(lutt_params, k_plus_q ) + omega_zero(iphon)*unitary_33 )
               temporary1 = invmat3( intsum(:,:,ik, idir) )
               !Use the eigenvector from the second k point along the direction.
               temporary2 = matmul( temporary1(1:3,1:3), eigenvec(1:3,ieig,2,idir) )
               temporary3 = dot_product( eigenvec(1:3,ieig,2,idir), temporary2(1:3) )
               sigma(ik,ieig,idir) = sigma(ik,ieig,idir) &
                   + signpm*piinv*invepsilonstar(iphon)*omega_zero(iphon) &
                   *temporary3*(omega_zero(iphon)/abs(lutt_params(1)))**half/(cos(xi))**two * weight_qdir(iqdir)
             enddo  ! iqdir
           enddo ! ixi
         enddo ! ik  
       enddo ! ieig 
     enddo ! idir
   endif
 enddo ! iphon

!Normalize self-energy integral
 sigma = sigma/nq_factor

 ABI_MALLOC(d2sigmadk2,(deg_dim,nkdir))

 d2sigmadk2 = zero

 do idir = 1,nkdir
   do ieig=1,3
     !Compute effective mass for the correct eigenvector
     invemass_ieig(ieig,idir)=2.0_dp*eigenval(ieig,2,idir)/krange**2.0_dp
     !Summation over k-range 
     !Here consider a small BZ region around Gamma
     !A convergene study might be performed around this value
     !Finite difference for the 2nd derivative of the self energy         
     d2sigmadk2(ieig,idir) = 2.0_dp/krange**2.0_dp* (four/three)* &
&      ( sigma(2,ieig,idir) - sigma(1,ieig,idir) - (sigma(3,ieig,idir) - sigma(1,ieig,idir))/16.0_dp )
   enddo
 enddo

 ABI_MALLOC(invpolmass,(deg_dim,nkdir))

 do idir = 1,nkdir
   do ieig=1,3
    !Calculate the inverse polaron mass
    invpolmass(ieig,idir) =  invemass_ieig(ieig,idir) + d2sigmadk2( ieig, idir)
!DEBUG
!     write(std_out,*)' idir, ieig, invemass(ieig,idir),invemass_ieig(ieig,idir), d2sigmadk2( ieig, idir), invpolmass(ieig,idir)=',&
!&     idir, ieig, invemass(ieig,idir), invemass_ieig(ieig,idir), d2sigmadk2( ieig, idir), invpolmass(ieig,idir)
!ENDDEBUG
   enddo
 enddo

!Print inverse electronic effective masses in the output
 write(ab_out,'(a)')'--------------------------------------------------------------------------------'
 write(ab_out,'(a)')'   Polaron properties from the generalized Froehlich model'
 write(ab_out,'(a)')'--------------------------------------------------------------------------------'
 write(ab_out,'(a)')'   Polar modes'
 write(ab_out,'(a)')'   ##      Frequency(meV)            Epsilon*'
 do imode = 1,3*cryst%natom
   if(invepsilonstar(imode) > tol10) then
    write(ab_out,'(2x,i3,5x,f15.6,5x,f15.6)')imode,omega_zero(imode)*Ha_eV*1000.0_dp,1.0_dp/invepsilonstar(imode)
   endif
 enddo
 write(ab_out,'(a)')' '
 write(ab_out,'(a,f10.2)')'   ZPR (meV): ',sigma(1,1,1)*Ha_eV*1000.0_dp
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')'   Electronic effective mass (a.u.) along 3 directions'
 write(ab_out,'(a, 3f15.6)')'    Direction 100:         ',one/invemass_ieig(:,1)
 write(ab_out,'(a, 3f15.6)')'    Direction 110:         ',one/invemass_ieig(:,2)
 write(ab_out,'(a, 3f15.6)')'    Direction 111:         ',one/invemass_ieig(:,3)

!Print inverse polaron effective masses in the output
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')'   Polaron effective mass (a.u.) along 3 directions'
 write(ab_out,'(a, 3f15.6)')'    Direction 100:         ',one/invpolmass(:,1)
 write(ab_out,'(a, 3f15.6)')'    Direction 110:         ',one/invpolmass(:,2)
 write(ab_out,'(a, 3f15.6)')'    Direction 111:         ',one/invpolmass(:,3)
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')'   Sum rule of inverse polaron masses check-up (for convergence purposes):'
 write(ab_out,'(a, 3f15.6)')'    Direction 100:         ',SUM(invpolmass(:,1))
 write(ab_out,'(a, 3f15.6)')'    Direction 110:         ',SUM(invpolmass(:,2))
 write(ab_out,'(a, 3f15.6)')'    Direction 111:         ',SUM(invpolmass(:,3))
 

 ABI_FREE(lutt_unit_kdir)
 ABI_FREE(eigenvec)
 ABI_FREE(eigenval)

 ABI_FREE(weight_qdir)
 ABI_FREE(unit_qdir)
 ABI_FREE(polarity_qdir)
 ABI_FREE(phfrq_qdir)
 ABI_FREE(dielt_qdir)
 ABI_FREE(invepsilonstar)
 ABI_FREE(omega_zero)

 ABI_FREE(intsum)
 ABI_FREE(sigma)
 ABI_FREE(d2sigmadk2)
 ABI_FREE(invpolmass)
 ABI_FREE(invemass)
 ABI_FREE(invemass_ieig)

 end subroutine polaronmass
!!***

function invmat3(A) result(B)
    !Invert a 3x3 matrix  directly - Faster than calling LAPACK 
    real(dp), intent(in) :: A(3,3)   !! Matrix
    real(dp)             :: B(3,3)   !! Inverse matrix
    real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function invmat3

function hamiltonian(luttin, V) result(H)
    !Build the Hamiltonian that enters the self-energy integral
    real(dp), intent(in) :: luttin(3)   !! Luttinger parameters
    real(dp), intent(in) :: V(3)   !! wave-vector
    real(dp)             :: H(3,3)

    H(1,1) = luttin(1)*V(1)**2.0_dp + luttin(2)*( V(2)**2.0_dp + V(3)**2.0_dp )
    H(1,2) = luttin(3)*V(1)*V(2)
    H(1,3) = luttin(3)*V(1)*V(3)
    H(2,1) = H(1,2)
    H(2,2) = luttin(1)*V(2)**2.0_dp + luttin(2)*( V(1)**2.0_dp + V(3)**2.0_dp ) 
    H(2,3) = luttin(3)*V(2)*V(3)
    H(3,1) = H(1,3) 
    H(3,2) = H(2,3) 
    H(3,3) = luttin(1)*V(3)**2.0_dp + luttin(2)*( V(1)**2.0_dp + V(2)**2.0_dp ) 
    
end function hamiltonian

end module m_frohlichmodel
!!***
