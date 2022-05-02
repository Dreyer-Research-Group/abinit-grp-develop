!!****m* ABINIT/opernlb_ylm_met
!! NAME
!! opernlb_ylm_met
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   from projected scalars to reciprocal space.
!! * Operate with the non-local projectors and the overlap matrix,
!!   from projected scalars to reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module m_opernlb_ylm_met

 use defs_basis
 use m_profiling_abi
 use m_errors
#if defined HAVE_OPENMP
 use OMP_LIB
#endif

 implicit none

 private
!!***

 public :: opernlb_ylm_met
 !integer,public,save :: opernla_counter = -1

contains
!!***

!!****f* ABINIT/opernlb_ylm_met
!! NAME
!! opernlb_ylm_met
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   from projected scalars to reciprocal space.
!! * Operate with the non-local projectors and the overlap matrix,
!!   from projected scalars to reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=chooses possible output (see below)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_dgxdt(ndgxdt_fac) = used only when cplex = 1
!!    cplex_dgxdt(i)=1 if dgxdt(1,i,:,:) is real, 2 if it is pure imaginary
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfacrelated to Sij (overlap)
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5, 51, 52 and signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  kpg(npw,nkpg)=(k+G) components (if nkpg=3)
!!  matblk=dimension of the array ph3d
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0 or 3)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_in>
!!      if (choice=2) <G|dV_nonlocal/d(atm. pos)|vect_in>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_in>
!!      if (choice=5) <G|dV_nonlocal/d(k)|vect_in>
!!      if (choice=51) <G|d(right)V_nonlocal/d(k)|vect_in>
!!      if (choice=52) <G|d(left)V_nonlocal/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)V_nonlocal/d(k)|vect_in>
!!  if (paw_opt=2)
!!    vectout(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_in> (note: not including <G|I|c>)
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm. pos)|vect_in>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=5) <G|d[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1) <G|I+S|vect_in> (note: not including <G|I|c>)
!!      if (choice=2) <G|dS/d(atm. pos)|vect_in>
!!      if (choice=3) <G|dS/d(strain)|vect_in>
!!      if (choice=5) <G|dS/d(k)|vect_in>
!!      if (choice=51) <G|d(right)S/d(k)|vect_in>
!!      if (choice=52) <G|d(left)S/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)S/d(k)|vect_in>
!!      if (choice=7) <G|sum_i[p_i><p_i]|vect_in>
!!
!! NOTES
!! 1-The openMP version is different from the standard version:
!!   the standard version is more effifient on one CPU core.
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!! 
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

subroutine opernlb_ylm_met(choice,cplex,cplex_dgxdt,cplex_fac,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
&                      ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
&                      nspinor,paw_opt,ph3d,svect,ucvol,vect,&
&                      dgxdtqfac,gxqfac,qpt,ffnlout,dimffnlout,npwout,kpgout,nkpgout,ph3dout,vectk)

 implicit none

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nincat
 integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
!AMSrev
 integer,intent(in) :: dimffnlout,npwout,nkpgout 
!arrays
 integer,intent(in) ::  cplex_dgxdt(ndgxdtfac),indlmn(6,nlmn),nloalg(3)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(2,npw*nspinor*(paw_opt/3)),vect(2,npw*nspinor)
!AMSrev
 real(dp),intent(in) :: dgxdtqfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxqfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: qpt(3),ffnlout(npwout,dimffnlout,nlmn)
 real(dp),intent(in) :: kpgout(npwout,nkpgout),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectk(2,npwout*nspinor)
!Local variables-------------------------------
!Arrays
!scalars
 integer :: fdb,fdf,ia,iaph3d,ic,ii,il,ilmn,ipw,ipwshft,ispinor,jc,nthreads
 real(dp) :: scale,wt
 logical :: parity
!arrays
 integer,parameter :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 real(dp),allocatable :: dgxdtfac_(:,:,:),dgxdtfacs_(:,:,:),gxfac_(:,:),gxfacs_(:,:)
 complex(dpc),allocatable :: ztab(:)
!AMSrev
 real(dp),allocatable :: dgxdtqfac_(:,:,:),gxqfac_(:,:)
 complex(dpc),allocatable :: ztab1(:)
 integer :: ipwshftk

! *************************************************************************

 DBG_ENTER("COLL")

!Nothing to do when choice=4, 6 or 23
 if (choice==4.or.choice==6.or.choice==23) return

!DDK not compatible with istwkf > 1 
 if(cplex==1.and.any(cplex_dgxdt(:)==2))then 
   ABI_BUG("opernlb_ylm+ddk not compatible with istwfk>1")
 end if

!Inits
 wt=four_pi/sqrt(ucvol)
 nthreads=1
#if defined HAVE_OPENMP
 nthreads=OMP_GET_NUM_THREADS()
#endif

 if (paw_opt/=3) then
   ABI_MALLOC(gxfac_,(2,nlmn))
   gxfac_(:,:)=zero
   if (choice>1) then
     ABI_MALLOC(dgxdtfac_,(2,ndgxdtfac,nlmn))
     if(ndgxdtfac>0) dgxdtfac_(:,:,:)=zero
   end if
 end if
!AMSrev[
  if (paw_opt/=3) then
   ABI_MALLOC(gxqfac_,(2,nlmn))
   gxqfac_(:,:)=zero
   if (choice>1) then
     ABI_MALLOC(dgxdtqfac_,(2,ndgxdtfac,nlmn))
     if(ndgxdtfac>0) dgxdtqfac_(:,:,:)=zero
   end if
 end if
!AMSrev]
 if (paw_opt>=3) then
   ABI_MALLOC(gxfacs_,(2,nlmn))
   gxfacs_(:,:)=zero
   if (choice>1) then
     ABI_MALLOC(dgxdtfacs_,(2,ndgxdtfac,nlmn))
     if (ndgxdtfac>0) dgxdtfacs_(:,:,:)=zero
   end if
 end if

!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 if (nthreads==1) then

!  Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw
! AMSrev
     ipwshftk=(ispinor-1)*npwout
   

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!      Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic =  cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do                 
               end if
             end if 
           end do
         end if
       end if

!AMSrev [
!      Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxqfac_(1:cplex_fac,ilmn)=scale*gxqfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxqfac_(2,ilmn)=zero
           else
             gxqfac_(2,ilmn)=-scale*gxqfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxqfac_(1,ilmn)=scale*gxqfac(2,ilmn,ia,ispinor)
             else
               gxqfac_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtqfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtqfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtqfac_(ic,ii,ilmn)=scale*dgxdtqfac(1,ii,ilmn,ia,ispinor)
                   dgxdtqfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtqfac_(2,ii,ilmn)=-scale*dgxdtqfac(1,ii,ilmn,ia,ispinor)
                   dgxdtqfac_(1,ii,ilmn)= scale*dgxdtqfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic =  cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtqfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtqfac_(jc,ii,ilmn)=-scale*dgxdtqfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtqfac_(jc,ii,ilmn)= scale*dgxdtqfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

!AMSrev]

!      Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

!AMSrev
       ABI_MALLOC(ztab,(npw))
       ABI_MALLOC(ztab1,(npwout))

!      Compute <g|Vnl|c> (or derivatives) for each plane wave:

       if (paw_opt/=3) then

!AMSrev
         ztab(:)=czero
         ztab1(:)=czero

!        ------
         if (choice==1) then ! <g|Vnl|c>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!AMSrev[
! In eq. 47 (Max notes):
!- the complex cognugate of \zeta functions does not matter because these functions are real(real spherical harmonics)
!- the factor 2*pi is already present only in dgxdtqfac_ and dgxdtfac_ 
!- a factor '+i' is already present in dgxdtqfac_ and dgxdtfac_ (as for the phonons, eq. 55 Gonze) 
!        ------
         if (choice==22) then ! derivative w.r.t. metric wave
           !write(*,*) 'ams: ... : npw, idir', npw, npwout, idir



           do ilmn=1,nlmn
             ztab(:)=ztab(:)+two_pi*qpt(idir)/two*ffnl(:,1,ilmn)*cmplx(-gxfac_(2,ilmn),gxfac_(1,ilmn),kind=dp)
           end do

           do ilmn=1,nlmn
             ztab1(:)=ztab1(:)+two_pi*qpt(idir)/two*ffnlout(:,1,ilmn)*cmplx(-gxqfac_(2,ilmn),gxqfac_(1,ilmn),kind=dp)

!             ztab(:)=ztab(:)+two_pi*qpt(idir)/two*ffnl(:,1,ilmn)*cmplx(-gxqfac_(2,ilmn),gxqfac_(1,ilmn),kind=dp)
           end do

           do ilmn=1,nlmn
             do ipw=1,npw
               ztab(ipw)=ztab(ipw)+two_pi*kpg(ipw,idir)*ffnl(ipw,1,ilmn)*cmplx(-gxfac_(2,ilmn),gxfac_(1,ilmn),kind=dp)
             end do
           end do

           do ilmn=1,nlmn
             do ipw=1,npwout
               ztab1(ipw)=ztab1(ipw)+two_pi*kpgout(ipw,idir)*ffnlout(ipw,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
             end do
           end do

           do ilmn=1,nlmn
             ztab1(:)=ztab1(:)-ffnlout(:,1,ilmn)*cmplx(dgxdtqfac_(1,1,ilmn),dgxdtqfac_(2,1,ilmn),kind=dp)

!             ztab(:)=ztab(:)-ffnl(:,1,ilmn)*cmplx(dgxdtqfac_(1,1,ilmn),dgxdtqfac_(2,1,ilmn),kind=dp)
           end do

           do ilmn=1,nlmn
             ztab1(:)=ztab1(:)+ffnlout(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do


         end if
!AMSrev]

!        ------
         if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&               *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp)&
&               -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&               -ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
         end if

!        ------
         if (choice==5) then ! full derivative w.r.t. k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&             +ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
!                                                -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
&             ffnl(:,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
&             ffnl(:,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if


!        ------
!AMSrev
!         write(211,*) 'aa', iaph3d
!         do ipw=1,npw
!           write(211,*) ph3d(:,ipw,iaph3d)
!         end do
!         write(210,*) 'aa', iaph3d
!         do ipw=1,npwout
!           write(210,*) ph3dout(:,ipw,iaph3d)
!         end do
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         vect(1,1+ipwshft:npw+ipwshft)=vect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         vect(2,1+ipwshft:npw+ipwshft)=vect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
!AMSrev
         ztab1(:)=ztab1(:)*cmplx(ph3dout(1,:,iaph3d),-ph3dout(2,:,iaph3d),kind=dp) !!the phase must be changed!!
         vectk(1,1+ipwshftk:npwout+ipwshftk)=vectk(1,1+ipwshftk:npwout+ipwshftk)+real(ztab1(:))
         vectk(2,1+ipwshftk:npwout+ipwshftk)=vectk(2,1+ipwshftk:npwout+ipwshftk)+aimag(ztab1(:))
       end if

!      Compute <g|S|c> (or derivatives) for each plane wave:

       if (paw_opt>=3) then

         ztab(:)=czero

!        ------
         if (choice==1) then ! <g|S|c>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&               *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
&               -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&               -ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end if
         end if

!        ------
         if (choice==5) then ! full derivative w.r.t. k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&             +ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
!                                                -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
&             ffnl(:,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
&             ffnl(:,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if


!        ------
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
       end if

       ABI_FREE(ztab)
       ABI_FREE(ztab1)

!      End loop on atoms
     end do
   end do !  End loop on spinors


!  ==========================================================================
!  ========== OPENMP VERSION ================================================
!  ==========================================================================
 else

! AMSrev
   write(*,*) 'ams: opermlb_ylm_met: ATT!! entering in OPENMP loop, which was not modified, insted of standard loop'

!  Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!      Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic =  cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do                 
               end if
             end if 
           end do
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

!      Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do        
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

       ABI_MALLOC(ztab,(npw))

!      Compute <g|Vnl|c> (or derivatives) for each plane wave:
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ipw,ilmn,fdf,fdb)

!        ------
         if (choice==1) then ! <g|Vnl|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
&                 *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp) &
&                 -ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
&                 -ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if


!        ------
         else if (choice==5) then ! full derivative w.r.t. k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
&               +ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,2,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
!                                                     -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) &
&               -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE 
         end if

!        ------
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw+ipwshft)=vect(1,ipw+ipwshft)+real(ztab(ipw))
           vect(2,ipw+ipwshft)=vect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO

!$OMP END PARALLEL
       end if

!      Compute <g|S|c> (or derivatives) for each plane wave:
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,ipw,fdf,fdb)

!        ------
         if (choice==1) then ! <g|S|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==3) then ! derivative w.r.t. strain
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
&                 *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
&                 -ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
&                 -ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if

!        ------
         else if (choice==5) then ! full derivative w.r.t. k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
&               +ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,2,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi> -
!          <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) &
&               -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE 
         end if


!        ------
!        The OMP WORKSHARE directive doesn't have a good performance with Intel Compiler
!        !$OMP WORKSHARE
!        ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
!        !$OMP END WORKSHARE
!        !$OMP WORKSHARE
!        svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
!        svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
!        !$OMP END WORKSHARE
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw+ipwshft)=svect(1,ipw+ipwshft)+real(ztab(ipw))
           svect(2,ipw+ipwshft)=svect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO
!$OMP END PARALLEL
       end if

       ABI_FREE(ztab)

!      End loop on atoms
     end do
!    End loop on spinors
   end do

!  ==========================================================================
 end if

 if (paw_opt/=3) then
   ABI_FREE(gxfac_)
   if (choice>1) then
     ABI_FREE(dgxdtfac_)
   end if
 end if
!AMSrev[
 if (paw_opt/=3) then
   ABI_FREE(gxqfac_)
   if (choice>1) then
     ABI_FREE(dgxdtqfac_)
   end if
 end if
!AMSrev]
 if (paw_opt>=3) then
   ABI_FREE(gxfacs_)
   if (choice>1) then
     ABI_FREE(dgxdtfacs_)
   end if
 end if

 DBG_EXIT("COLL")

#if !defined HAVE_OPENMP
!Fake use of unused variable
 if (.false.) write(std_out,*) ipw
#endif
 
end subroutine opernlb_ylm_met
!!***

end module m_opernlb_ylm_met
!!***
