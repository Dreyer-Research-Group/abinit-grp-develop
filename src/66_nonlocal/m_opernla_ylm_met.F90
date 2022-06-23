!!****m* ABINIT/m_opernla_ylm_met
!! NAME
!! m_opernla_ylm
!!
!! FUNCTION
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
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

module m_opernla_ylm_met

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
#if defined HAVE_OPENMP
 use OMP_LIB
#endif
 use m_time,        only : timab


 implicit none

 private
!!***

 public :: opernla_ylm_met
 !integer,public,save :: opernla_counter = -1

contains
!!***

!!****f* ABINIT/opernla_ylm_met
!! NAME
!! opernla_ylm_met
!!
!! FUNCTION
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=chooses possible output:
!!         if choice>=0: compute projected scalars
!!         if choice<0: same as choice>0 but use already computed projected scalars
!!         if ABS(choice)>1, then compute additional quantities:
!!           2: compute projected scalars and derivatives wrt atm pos.
!!           3: compute projected scalars and derivatives wrt strains
!!           23: compute projected scalars, derivatives wrt atm pos. and derivatives wrt strains
!!           4, 24: compute projected scalars, derivatives wrt atm pos.
!!                  and 2nd derivatives wrt atm pos.
!!           5,51,52: compute projected scalars and derivatives wrt wave vector k
!!           53: compute projected scalars and derivatives wrt wave vector k in direction idir+1 and idir-1
!!           54: compute projected scalars, deriv. wrt atm pos., deriv. wrt wave vector k
!!               and 2nd derivatives wrt left wave vector k and atm pos.
!!           6: compute projected scalars, derivatives wrt atm pos., derivatives wrt strains,
!!              2nd derivatives wrt 2 strains and derivatives wrt strain and atm pos.
!!           7: not available
!!           8: compute projected scalars, derivatives wrt wave vector k
!!              and 2nd derivatives wrt 2 wave vectors k
!!  cplex=1 if <p_lmn|c> scalars are real or pure imaginary (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpg(npw,nkpg)=(k+G) components          for ikpg=1...3   (if nkpg=3 or 9)
!!       [(k+G)_a].[(k+G)_b] quantities for ikpg=4...9   (if nkpg=9)
!!  matblk=dimension of the array ph3d
!!  mpi_enreg=informations about MPI parallelization
!!  ndgxdt=second dimension of dgxdt
!!  nd2gxdt=second dimension of d2gxdt
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0, 3 or 9)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  signs=chooses possible output:
!!   signs=1: compute derivatives in all directions
!!   signs=2: compute derivative in direction IDIR only
!!            compatible only with 1st-order derivatives and "single" derivatives
!!  ucvol=unit cell volume (bohr^3)
!!  vect(2,npw*my_nspinor)=starting vector in reciprocal space
!!  AMSrev added:
!!  kpgout(npwout,nkpgout)=k+q+G
!!  ffnlout(npwout,dimffnlout,nlmn)=ffnl(k+q+G)
!!  dimffnlout
!!  npwout
!!  nkpgout
!!
!! OUTPUT
!!  if (choice>1) dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=
!!     gradients of projected scalars wrt coords  (choice=2, 23, 4, 54 or 6)
!!                                    wrt strains (choice=3 or 23)
!!                                    wrt k wave vect. (choice=5, 51, 52, 53, 54, 8)
!!                AMSrev -->          wrt metric perturbation (choice=99):
!!                                        In this case I have also as output
!!                                        gxq(cplex,nlmn,nincat,nspinor)=\sum_G f(k+q+G) c(G) exp(..)
!!                                        dgxdtq(cplex,nlmn,nincat,nspinor)=\sum_G f(k+q+G) (k+G+q) c(G) exp(..)
!!  if (choice=4, 24, 54, 6 or 8) d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=
!!     2nd grads of projected scalars wrt 2 coords (choice=4 or 24)
!!                                    wrt coords & k wave vect. (choice=54)
!!                                    wrt coords & strains (choice=6)
!!                                    wrt 2 k wave vect. (choice=8)
!!     only compatible with signs=1
!!  cplex_dgxdt(ndgxdt) = used only when cplex = 1
!!             cplex_dgxdt(i) = 1 if dgxdt(1,i,:,:)   is real, 2 if it is pure imaginary
!!  cplex_d2gxdt(nd2gxdt) = used only when cplex = 1
!!             cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
!!
!! SIDE EFFECTS
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars - input if choice<0, output if choice>=0
!!
!! NOTES
!! 1-The openMP version is different from the standard version:
!!   the standard version is more effifient on one CPU core.
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine opernla_ylm_met(choice,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnl,d2gxdt,dgxdt,ffnl,gx,&
&       ia3,idir,indlmn,istwf_k,kpg,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg,nlmn,&
&       nloalg,npw,nspinor,ph3d,signs,ucvol,vect,&
&       gxq,dgxdtq,kpgout,npwout,nkpgout,ffnlout,dimffnlout,vectk,ph3dout)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,dimffnl,ia3,idir,istwf_k,matblk,nd2gxdt
 integer,intent(in) :: ndgxdt,nincat,nkpg,nlmn,npw,nspinor,signs
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!! AMSrev
 integer,intent(in) :: npwout,dimffnlout,nkpgout
!arrays
 integer,intent(in) :: indlmn(6,nlmn),nloalg(3) !AMSrev: nloalg had extent 5 in ab7...
 integer,intent(out) :: cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(in) :: vect(2,npw*nspinor)
 real(dp),intent(out) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor),gx(cplex,nlmn,nincat,nspinor)
!!AMSrev
 real(dp),intent(out) :: gxq(cplex,nlmn,nincat,nspinor),dgxdtq(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: kpgout(npwout,nkpgout),ffnlout(npwout,dimffnlout,nlmn)
 real(dp),intent(in) :: vectk(2,npwout*nspinor)
 real(dp),intent(in) :: ph3dout(2,npwout,matblk)
!Local variables-------------------------------
!no_abirules
 integer :: choice_,ffnl_dir1,ffnl_dir2,gama,gamb,gamc,gamd,i1,i2,i3,i4,i5,i6,ia,iaph3d
 integer :: ierr,il,ilmn,ipw,ipw0,ipwshft,ishift,ispinor,jpw,mu
 integer :: mua,mub,nthreads,nua1,nua2,nub1,nub2
 real(dp), parameter :: two_pi2=two_pi*two_pi
 real(dp) :: aux_i,aux_i2,aux_i3,aux_i4
 real(dp) :: aux_r,aux_r2,aux_r3,aux_r4
 real(dp) :: buffer_i,buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6
 real(dp) :: buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if
 real(dp) :: buffer_r,buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6
 real(dp) :: buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf
 real(dp) :: kpga,kpgb,kpgc,kpgd,scale,scale2,wt
 logical :: check,parity
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 real(dp) :: tsec(2)
 real(dp),allocatable :: scali(:),scalr(:),scalari(:,:),scalarr(:,:)
!AMSrev 
 integer :: ipwshftout
 real(dp),allocatable :: scalki(:),scalkr(:)

! *************************************************************************

 if (choice==-1) return

!Useful variables
 choice_=abs(choice)
 wt=four_pi/sqrt(ucvol);if (cplex==1) wt=2.d0*wt
 ipw0=1;if (istwf_k==2.and.mpi_enreg%me_g0==1) ipw0=2
 cplex_dgxdt(:)  = 0 ; if (cplex == 1) cplex_dgxdt(:)  = 1
 cplex_d2gxdt(:) = 0 ; if (cplex == 1) cplex_d2gxdt(:) = 1
 nthreads=1

!Check compatibility of options
!signs=1, almost all choices
 if (signs==1) then
   check=(choice_==0.or.choice_==1.or.choice_==2.or.choice_==3.or.choice_==4.or.&
&   choice_==23.or.choice_==24.or.choice_==5.or.choice_==51.or.choice_==52.or.&
&   choice_==53.or.choice_==54.or.choice_==6.or.choice_==8)
   ABI_CHECK(check,'BUG: choice not compatible (for signs=1)')
 end if
!signs=2, only "single derivative" choices
!AMSrev
 if (signs==2) then
   check=(choice_==99.or.choice_==0.or.choice_==1.or.choice_==2.or.choice_==3.or.&
&   choice_==5.or.choice_==51.or.choice_==52.or.choice_==53)
   ABI_CHECK(check,'BUG: signs=2 not compatible with this choice')
 end if
!1<=idir<=6 is  required when choice=3 and signs=2
 if (choice_==3.and.signs==2) then
   check=(idir>=1.and.idir<=6)
   ABI_CHECK(check,'BUG: choice=3 and signs=2 requires 1<=idir<=6')
!  idir=0 is required when choice=3 and signs=2
 else if (choice_==6.and.signs==1) then
   check=(idir==0)
   ABI_CHECK(check,'BUG: choice=6 and signs=1 requires idir = 0')
 else
!  signs=2 requires 1<=idir<=3 when choice>1
   check=(signs/=2.or.choice<=1.or.(idir>=1.and.idir<=3))
   ABI_CHECK(check,'BUG: signs=2 requires 1<=idir<=3')
 end if

#if defined HAVE_OPENMP
 nthreads=OMP_GET_NUM_THREADS()
#endif

!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
!Parallelization OpenMP on nlmn loops
 if (nthreads==1) then

!$OMP PARALLEL PRIVATE(il,ilmn,ipw,jpw,parity,buffer_r,buffer_i, &
!$OMP ffnl_dir1,ffnl_dir2,scale,scale2, &
!$OMP buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6, &
!$OMP buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6, &
!$OMP aux_i,aux_i2,aux_i3,aux_i4,aux_r,aux_r2,aux_r3,aux_r4, &
!$OMP buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if, &
!$OMP buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf), &
!$OMP PRIVATE(ispinor,ipwshft,ia,iaph3d,i1,i2,i3)

!  Loop on spinorial components
   do ispinor =1,nspinor
     ipwshft=(ispinor-1)*npw
!AMSrev
     ipwshftout=(ispinor-1)*npwout

!    Allocate work space
!$OMP SECTIONS
!$OMP SECTION
     ABI_MALLOC(scali,(npw))
!AMSrev
     ABI_MALLOC(scalki,(npwout))
!$OMP SECTION
     ABI_MALLOC(scalr,(npw))
!AMSrev
     ABI_MALLOC(scalkr,(npwout))
!$OMP END SECTIONS

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!      Compute Sum_g[c(g).exp(2pi.i.g.R)]
!$OMP DO
       do ipw=ipw0,npw
         jpw=ipw+ipwshft
         scalr(ipw)=(vect(1,jpw)*ph3d(1,ipw,iaph3d)-vect(2,jpw)*ph3d(2,ipw,iaph3d))
         scali(ipw)=(vect(2,jpw)*ph3d(1,ipw,iaph3d)+vect(1,jpw)*ph3d(2,ipw,iaph3d))
       end do
!$OMP END DO
!$OMP SINGLE
       if (ipw0==2) then
         scalr(1)=half*vect(1,1+ipwshft)*ph3d(1,1,iaph3d)
         scali(1)=half*vect(1,1+ipwshft)*ph3d(2,1,iaph3d)
       end if
!$OMP END SINGLE
!AMSrev[
       do ipw=ipw0,npwout
         jpw=ipw+ipwshftout
         scalkr(ipw)=(vectk(1,jpw)*ph3dout(1,ipw,iaph3d)-vectk(2,jpw)*ph3dout(2,ipw,iaph3d))
         scalki(ipw)=(vectk(2,jpw)*ph3dout(1,ipw,iaph3d)+vectk(1,jpw)*ph3dout(2,ipw,iaph3d))
       end do
!AMSrev]
!AMSrev[
       if (ipw0==2) then
         scalkr(1)=half*vectk(1,1+ipwshftout)*ph3dout(1,1,iaph3d)
         scalki(1)=half*vectk(1,1+ipwshftout)*ph3dout(2,1,iaph3d)
       end if
!AMSrev]

!      --------------------------------------------------------------------
!      ALL CHOICES:
!      Accumulate Gx
!      --------------------------------------------------------------------
       if (choice>=0) then
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r = zero ; buffer_i = zero
             do ipw=1,npw
               buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
               buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
             end do
             if (parity) then
               gx(1,ilmn,ia,ispinor) = scale*buffer_r ; gx(2,ilmn,ia,ispinor) = scale*buffer_i
             else
               gx(1,ilmn,ia,ispinor) =-scale*buffer_i ; gx(2,ilmn,ia,ispinor) = scale*buffer_r
             end if
           else
             if (parity) then
               buffer_r =  zero
               do ipw=1,npw
                 buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
               end do
               gx(1,ilmn,ia,ispinor) = scale*buffer_r
             else
               buffer_i = zero
               do ipw=1,npw
                 buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
               end do
               gx(1,ilmn,ia,ispinor) =-scale*buffer_i
             end if
           end if
         end do
!$OMP END DO
       end if

! AMSrev[
!      --------------------------------------------------------------------
!      ALL CHOICES:
!      Accumulate Gxp
!      --------------------------------------------------------------------

       if (choice>=0) then
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r = zero ; buffer_i = zero
             do ipw=1,npwout
! AMSrev
               buffer_r = buffer_r + scalkr(ipw)*ffnlout(ipw,1,ilmn)
               buffer_i = buffer_i + scalki(ipw)*ffnlout(ipw,1,ilmn)
             end do
             if (parity) then
               gxq(1,ilmn,ia,ispinor) = scale*buffer_r ; gxq(2,ilmn,ia,ispinor) = scale*buffer_i
             else
               gxq(1,ilmn,ia,ispinor) =-scale*buffer_i ; gxq(2,ilmn,ia,ispinor) = scale*buffer_r
             end if
           else
             if (parity) then
               buffer_r =  zero
               do ipw=1,npwout
                 buffer_r = buffer_r + scalkr(ipw)*ffnlout(ipw,1,ilmn)
               end do
               gxq(1,ilmn,ia,ispinor) = scale*buffer_r
             else
               buffer_i = zero
               do ipw=1,npwout
                 buffer_i = buffer_i + scalki(ipw)*ffnlout(ipw,1,ilmn)
               end do
               gxq(1,ilmn,ia,ispinor) =-scale*buffer_i
             end if
           end if
         end do
!$OMP END DO
       end if
! AMSrev ]

!      --------------------------------------------------------------------
!      CHOICE= 2, 4, 6, 23, 24, 54  --  SIGNS= 1
!      Accumulate dGxdt --- derivative wrt atm pos. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and. &
&       (choice_==2.or.choice_==23.or.choice_==24.or.choice_==4.or.choice_==54.or.choice_==6)) then
         i1=1;i2=2;i3=3
         if (choice_==23.or.choice_==6) then
           i1=7;i2=8;i3=9
         end if
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_i1 = zero
             buffer_i2 = zero ; buffer_i3 = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
               buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
               buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_i3
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_r3
             end if
           else
             if (parity) then
               buffer_r1 = zero ; buffer_r2 = zero ; buffer_r3 = zero
               do ipw=1,npw
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
                 buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
                 buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_r3
             else
               buffer_i1 = zero ; buffer_i2 = zero ; buffer_i3 = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_i3
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 2  --  SIGNS= 2
!      Accumulate dGxdt --- derivative wrt atm pos. --- for direction IDIR
!      --------------------------------------------------------------------
! AMSrev
       if ((signs==2).and.((choice_==2).or.(choice_==99))) then
         i1=1
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,idir)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,idir)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scali(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scalr(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
             end if
           end if
         end do
!$OMP END DO
       end if

!AMSrev [ here there is factor i
!      --------------------------------------------------------------------
!      CHOICE= 99  --  SIGNS= 2
!      Accumulate dGxdtp --- derivative wrt atm pos. --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==99)) then
         i1=1
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npwout
               aux_r = scalkr(ipw)*ffnlout(ipw,1,ilmn)
               aux_i = scalki(ipw)*ffnlout(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpgout(ipw,idir)
               buffer_i1 = buffer_i1 + aux_r*kpgout(ipw,idir)
             end do
             if (parity) then
               dgxdtq(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdtq(2,i1,ilmn,ia,ispinor) = scale2*buffer_i1
             else
               dgxdtq(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdtq(2,i1,ilmn,ia,ispinor) = scale2*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npwout
                 buffer_r1 = buffer_r1 - scalki(ipw)*ffnlout(ipw,1,ilmn)*kpgout(ipw,idir)
               end do
               dgxdtq(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npwout
                 buffer_i1 = buffer_i1 + scalkr(ipw)*ffnlout(ipw,1,ilmn)*kpgout(ipw,idir)
               end do
               dgxdtq(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
             end if
           end if
         end do
!$OMP END DO
       end if

! AMSrev ]

!      --------------------------------------------------------------------
!      CHOICE= 3, 23 or 6  -- SIGNS= 1
!      Accumulate dGxdt --- derivative wrt strain --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==3.or.choice_==23.or.choice_==6)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=half*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_r4 = zero
             buffer_r5 = zero ; buffer_r6 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero ; buffer_i4 = zero
             buffer_i5 = zero ; buffer_i6 = zero
             do ipw=1,npw
               aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
               aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
               aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
               aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
               aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
               aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
               buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
               buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
               buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
               buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
               buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
               buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
               buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
               buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
               buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(2,i4,ilmn,ia,ispinor) =-scale2*buffer_i4
               dgxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(2,i5,ilmn,ia,ispinor) =-scale2*buffer_i5
               dgxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
               dgxdt(2,i6,ilmn,ia,ispinor) =-scale2*buffer_i6
             else
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(2,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(2,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_i6
               dgxdt(2,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
             end if
           else
             if (parity) then
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero ; buffer_r4 = zero
               buffer_r5 = zero ; buffer_r6 = zero
               do ipw=1,npw
                 aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
                 aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
                 aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
                 buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
                 buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
                 buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
                 buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
                 buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
                 buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
             else
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero ; buffer_i4 = zero
               buffer_i5 = zero ; buffer_i6 = zero
               do ipw=1,npw
                 aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
                 aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
                 aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
                 buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
                 buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
                 buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
                 buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_i6
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 3  --  SIGNS= 2
!      Accumulate dGxdt --- derivative wrt strain --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==3)) then
         i1=1
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 4, 24  -- SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 atm pos. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==4.or.choice_==24)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi2*scale
           if (cplex==2) then
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
               buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
               buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
               buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
               buffer_re = buffer_re - aux_r*kpg(ipw,8)
               buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
               buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
               buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
               buffer_id = buffer_id - aux_i*kpg(ipw,7)
               buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
               buffer_if = buffer_if - aux_i*kpg(ipw,9)
             end do
             if (parity) then
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale2*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale2*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_rf
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale2*buffer_if
             else
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_ia
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_ib
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_ic
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_id
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_ie
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_if
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale2*buffer_rf
             end if
           else
             if (parity) then
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
                 buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
                 buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
                 buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
                 buffer_re = buffer_re - aux_r*kpg(ipw,8)
                 buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               end do
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_rf
             else
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
                 buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
                 buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
                 buffer_id = buffer_id - aux_i*kpg(ipw,7)
                 buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
                 buffer_if = buffer_if - aux_i*kpg(ipw,9)
               end do
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_if
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 5, 53, 54, 8  --  SIGNS= 1
!      Accumulate dGxdt --- derivative wrt k --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==5.or.choice_==53.or.choice_==54.or.choice_==8)) then
         i1=1;i2=2;i3=3
         if (choice_==54) then
           i1=4;i2=5;i3=6
         end if
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
               buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
               buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
               buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_i3
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_r3
             end if
           else
             cplex_dgxdt(i1:i3) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
                 buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
                 buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
             else
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
                 buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
                 buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_r3
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 5, 51, 52 -- SIGNS= 2
!      Accumulate dGxdt --- derivative wrt k --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==5.or.choice_==51.or.choice_==52)) then
         i1=1
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
           else
             cplex_dgxdt(i1) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 53 -- SIGNS= 2
!      Accumulate dGxdt --- derivative wrt k --- for directions IDIR-1 & IDIR+1
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==53)) then
         i1=1;i2=2
!        Case 53 though need multiple directions and here ffnl will contain the
!        derivative information in locations 2,3, and 4 corresponding to idir = 1, 2,
!        and 3. Moreover, choice 53 needs the derivatives in direction idir+1 and idir-1.
!        The parameter vector ffnl_dir_dat contains the necessary translations in
!        locations 1,2 for idir=1; 3,4 for idir=2; and 5,6 for idir=3.
         ffnl_dir1 = ffnl_dir_dat(2*idir-1) ! idir+1 derivative
         ffnl_dir2 = ffnl_dir_dat(2*idir)   ! idir-1 derivative
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             buffer_r2 = zero ; buffer_i2 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
               buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
             end do
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_i2
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_r2
             end if
           else
             cplex_dgxdt(i1:i2) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
               buffer_i1 = zero ; buffer_i2 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
                 buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
             else
               buffer_r1 = zero ; buffer_r2 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
                 buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               end do
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE= 54,  --  SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt k and atm. pos --- in all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==54)) then
         ishift=0
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 mu=ishift+mub+3*(mua-1)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                 end if
               end do
             end do
           else
             cplex_d2gxdt(ishift+1:ishift+9) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
               do mua=1,3 ! atm. pos
                 do mub=1,3 ! k
                   mu=ishift+mub+3*(mua-1)
                   buffer_r = zero
                   do ipw=1,npw
                     buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             else
               do mua=1,3 ! atm. pos
                 do mub=1,3 ! k
                   mu=ishift+mub+3*(mua-1)
                   buffer_i = zero
                   do ipw=1,npw
                     buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = -scale2*buffer_i
                 end do
               end do
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE 6:
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 strains --- for all directions
!      Accumulate d2Gxdt --- 2nd derivative wrt strain and atm. pos --- for all directions
!        --------------------------------------------------------------------
       if ((signs==1).and.(choice_==6)) then
!$OMP SECTIONS
!$OMP SECTION
         ABI_MALLOC(scalarr,(npw,10))
!$OMP SECTION
         ABI_MALLOC(scalari,(npw,10))
!$OMP END SECTIONS
!$OMP DO PRIVATE(mu,mua,mub,ishift,nua1,nua2,nub1,nub2), &
!$OMP PRIVATE(il,scale,scale2,parity), &
!$OMP PRIVATE(gama,gamb,gamc,gamd), &
!$OMP PRIVATE(kpga,kpgb,kpgc,kpgd)
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           do mu=1,10
             do ipw=1,npw
               scalarr(ipw,mu)=scalr(ipw)*ffnl(ipw,mu,ilmn)
               scalari(ipw,mu)=scali(ipw)*ffnl(ipw,mu,ilmn)
             end do
           end do

!          ===== Accumulate 2nd derivative of Gx wrt two strains =====
           ishift=0;scale2=quarter*scale
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,6
                 mu=ishift+mua+6*(mub-1)
                 nua1=alpha(mua);nua2=beta(mua)
                 gama=gamma(nub1,nua1)
                 gamb=gamma(nub2,nua1)
                 gamc=gamma(nub1,nua2)
                 gamd=gamma(nub2,nua2)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   kpga=kpg(ipw,nub1)
                   kpgb=kpg(ipw,nub2)
                   kpgc=kpg(ipw,nua1)
                   kpgd=kpg(ipw,nua2)
                   buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                   + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                   + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_r = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                     + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_i = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                     + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                 end do
               end do
             end if
           end if

!            ===== Accumulate 2nd derivative of Gx wrt strain and atm pos. =====
           ishift=36;scale2=pi*scale
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,3
                 mu=ishift+mua+3*(mub-1)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   buffer_r = buffer_r + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                   buffer_i = buffer_i + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalari(ipw,1+nub2)*kpg(ipw,nub1))
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) =-scale2*buffer_r
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_i = zero
                   do ipw=1,npw
                     buffer_i = buffer_i + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalari(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_r = zero
                   do ipw=1,npw
                     buffer_r = buffer_r + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             end if
           end if
         end do
!$OMP END DO
!$OMP SECTIONS
!$OMP SECTION
         ABI_FREE(scalarr)
!$OMP SECTION
         ABI_FREE(scalari)
!$OMP END SECTIONS
       end if

!      --------------------------------------------------------------------
!      CHOICE= 8  --  SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 k wave vect. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==8)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
             do ipw=1,npw
               buffer_ra = buffer_ra + scalr(ipw)*ffnl(ipw,5,ilmn)
               buffer_rb = buffer_rb + scalr(ipw)*ffnl(ipw,6,ilmn)
               buffer_rc = buffer_rc + scalr(ipw)*ffnl(ipw,7,ilmn)
               buffer_rd = buffer_rd + scalr(ipw)*ffnl(ipw,8,ilmn)
               buffer_re = buffer_re + scalr(ipw)*ffnl(ipw,9,ilmn)
               buffer_rf = buffer_rf + scalr(ipw)*ffnl(ipw,10,ilmn)
               buffer_ia = buffer_ia + scali(ipw)*ffnl(ipw,5,ilmn)
               buffer_ib = buffer_ib + scali(ipw)*ffnl(ipw,6,ilmn)
               buffer_ic = buffer_ic + scali(ipw)*ffnl(ipw,7,ilmn)
               buffer_id = buffer_id + scali(ipw)*ffnl(ipw,8,ilmn)
               buffer_ie = buffer_ie + scali(ipw)*ffnl(ipw,9,ilmn)
               buffer_if = buffer_if + scali(ipw)*ffnl(ipw,10,ilmn)
             end do
             if (parity) then
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale*buffer_rf
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale*buffer_if
             else
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_ia
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_ib
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_ic
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale*buffer_id
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale*buffer_ie
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale*buffer_if
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale*buffer_rf
             end if
           else
             if (parity) then
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
               do ipw=1,npw
                 buffer_ra = buffer_ra + scalr(ipw)*ffnl(ipw,5,ilmn)
                 buffer_rb = buffer_rb + scalr(ipw)*ffnl(ipw,6,ilmn)
                 buffer_rc = buffer_rc + scalr(ipw)*ffnl(ipw,7,ilmn)
                 buffer_rd = buffer_rd + scalr(ipw)*ffnl(ipw,8,ilmn)
                 buffer_re = buffer_re + scalr(ipw)*ffnl(ipw,9,ilmn)
                 buffer_rf = buffer_rf + scalr(ipw)*ffnl(ipw,10,ilmn)
               end do
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale*buffer_rf
             else
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
               do ipw=1,npw
                 buffer_ia = buffer_ia + scali(ipw)*ffnl(ipw,5,ilmn)
                 buffer_ib = buffer_ib + scali(ipw)*ffnl(ipw,6,ilmn)
                 buffer_ic = buffer_ic + scali(ipw)*ffnl(ipw,7,ilmn)
                 buffer_id = buffer_id + scali(ipw)*ffnl(ipw,8,ilmn)
                 buffer_ie = buffer_ie + scali(ipw)*ffnl(ipw,9,ilmn)
                 buffer_if = buffer_if + scali(ipw)*ffnl(ipw,10,ilmn)
               end do
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale*buffer_if
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      END CHOICES
!      --------------------------------------------------------------------

     end do ! End loop on atoms

!    Deallocate temporary space
!$OMP SECTIONS
!$OMP SECTION
     ABI_FREE(scali)
!$OMP SECTION
     ABI_FREE(scalr)
!$OMP END SECTIONS

   end do !  End loop on spinorial components
!$OMP END PARALLEL


!  ==========================================================================
!  ========== OPENMP VERSION ================================================
!  ==========================================================================
!  Parallelization OpenMP on npw loops
 else
!AMSrev
   write(*,*) 'ams: opermla_ylm_met: ATT!! entering in OPENMP loop, which was not modified, insted of standard loop'

!$OMP PARALLEL PRIVATE(il,ilmn,ipw,jpw,parity,ffnl_dir1,ffnl_dir2), &
!$OMP PRIVATE(scale,scale2), &
!$OMP PRIVATE(ispinor,ipwshft,ia,iaph3d)

!  Loop on spinorial components
   do ispinor =1,nspinor
     ipwshft=(ispinor-1)*npw

!    Allocate work space
!$OMP SECTIONS
!$OMP SECTION
     ABI_MALLOC(scali,(npw))
!$OMP SECTION
     ABI_MALLOC(scalr,(npw))
!$OMP END SECTIONS

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!      Compute Sum_g[c(g).exp(2pi.i.g.R)]
!$OMP DO
       do ipw=ipw0,npw
         jpw=ipw+ipwshft
         scalr(ipw)=(vect(1,jpw)*ph3d(1,ipw,iaph3d)-vect(2,jpw)*ph3d(2,ipw,iaph3d))
         scali(ipw)=(vect(2,jpw)*ph3d(1,ipw,iaph3d)+vect(1,jpw)*ph3d(2,ipw,iaph3d))
       end do
!$OMP END DO
!$OMP SINGLE
       if (ipw0==2) then
         scalr(1)=half*vect(1,1+ipwshft)*ph3d(1,1,iaph3d)
         scali(1)=half*vect(1,1+ipwshft)*ph3d(2,1,iaph3d)
       end if
!$OMP END SINGLE

!      --------------------------------------------------------------------
!      ALL CHOICES
!      Accumulate Gx
!      --------------------------------------------------------------------
       if (choice>=0) then ! JWZ to check: I dont think 53 needs this
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r = zero ; buffer_i = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r,buffer_i) !!, &
!            !$OMP ORDERED
             do ipw=1,npw
               buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
               buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               gx(1,ilmn,ia,ispinor) = scale*buffer_r ; gx(2,ilmn,ia,ispinor) = scale*buffer_i
             else
               gx(1,ilmn,ia,ispinor) =-scale*buffer_i ; gx(2,ilmn,ia,ispinor) = scale*buffer_r
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_r =  zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r)
               do ipw=1,npw
                 buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               gx(1,ilmn,ia,ispinor) = scale*buffer_r
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_i = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i)
               do ipw=1,npw
                 buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               gx(1,ilmn,ia,ispinor) =-scale*buffer_i
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 2, 4, 6, 23, 24, 54  --  SIGNS= 1
!      Accumulate dGxdt --- derivative wrt atm pos. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and. &
&       (choice_==2.or.choice_==23.or.choice_==24.or.choice_==4.or.choice_==54.or.choice_==6)) then
         i1=1;i2=2;i3=3
         if (choice_==23.or.choice_==6) then
           i1=7;i2=8;i3=9
         end if
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_i1 = zero
             buffer_i2 = zero ; buffer_i3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3,buffer_i1,buffer_i2,buffer_i3),&
!$OMP PRIVATE(aux_r,aux_i)
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
               buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
               buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_i3
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_r3
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_r1 = zero ; buffer_r2 = zero ; buffer_r3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3), &
!$OMP PRIVATE(aux_r,aux_i)
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
                 buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
                 buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_r3
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_i1 = zero ; buffer_i2 = zero ; buffer_i3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1,buffer_i2,buffer_i3), &
!$OMP PRIVATE(aux_r,aux_i)
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
               end do
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_i3
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 2  --  SIGNS= 2
!      Accumulate dGxdt --- derivative wrt atm pos. --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==2)) then
         i1=1
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_i1), &
!$OMP PRIVATE(aux_r,aux_i)
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,idir)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,idir)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_r1
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_r1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1)
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scali(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_r1
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1)
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scalr(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_i1
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 3, 23 or 6  -- SIGNS= 1
!      Accumulate dGxdt --- derivative wrt strain --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==3.or.choice_==23.or.choice_==6)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=half*scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_r4 = zero
             buffer_r5 = zero ; buffer_r6 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero ; buffer_i4 = zero
             buffer_i5 = zero ; buffer_i6 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6), &
!$OMP REDUCTION(+:buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6), &
!$OMP PRIVATE(aux_r2,aux_r3,aux_r4,aux_i2,aux_i3,aux_i4)
             do ipw=1,npw
               aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
               aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
               aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
               aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
               aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
               aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
               buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
               buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
               buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
               buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
               buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
               buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
               buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
               buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
               buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(2,i4,ilmn,ia,ispinor) =-scale2*buffer_i4
               dgxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(2,i5,ilmn,ia,ispinor) =-scale2*buffer_i5
               dgxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
               dgxdt(2,i6,ilmn,ia,ispinor) =-scale2*buffer_i6
             else
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(2,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(2,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_i6
               dgxdt(2,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero ; buffer_r4 = zero
               buffer_r5 = zero ; buffer_r6 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6), &
!$OMP PRIVATE(aux_r2,aux_r3,aux_r4)
               do ipw=1,npw
                 aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
                 aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
                 aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
                 buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
                 buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
                 buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
                 buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
                 buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
                 buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               end do
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_r6
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero ; buffer_i4 = zero
               buffer_i5 = zero ; buffer_i6 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6), &
!$OMP PRIVATE(aux_i2,aux_i3,aux_i4)
               do ipw=1,npw
                 aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
                 aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
                 aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
                 buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
                 buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
                 buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
                 buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_i6
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 3  --  SIGNS= 2
!      Accumulate dGxdt --- derivative wrt strain --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==3)) then
         i1=1
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_i1)
             do ipw=1,npw
               buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_r1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1)
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1)
               do ipw=1,npw
                 buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 4, 24  -- SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 atm pos. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==4.or.choice_==24)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi2*scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP PRIVATE(aux_r,aux_i), &
!$OMP REDUCTION(+:buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf),&
!$OMP REDUCTION(+:buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if)
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
               buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
               buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
               buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
               buffer_re = buffer_re - aux_r*kpg(ipw,8)
               buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
               buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
               buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
               buffer_id = buffer_id - aux_i*kpg(ipw,7)
               buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
               buffer_if = buffer_if - aux_i*kpg(ipw,9)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale2*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale2*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_rf
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale2*buffer_if
             else
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_ia
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_ib
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_ic
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_id
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_ie
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_if
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale2*buffer_rf
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP PRIVATE(aux_r,aux_i), &
!$OMP REDUCTION(+:buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf)
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
                 buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
                 buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
                 buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
                 buffer_re = buffer_re - aux_r*kpg(ipw,8)
                 buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               end do
!$OMP END DO
!$OMP SINGLE
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale2*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale2*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale2*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale2*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale2*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale2*buffer_rf
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP PRIVATE(aux_i,aux_r), &
!$OMP REDUCTION(+:buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if)
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
                 buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
                 buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
                 buffer_id = buffer_id - aux_i*kpg(ipw,7)
                 buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
                 buffer_if = buffer_if - aux_i*kpg(ipw,9)
               end do
!$OMP END DO
!$OMP SINGLE
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale2*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale2*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale2*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale2*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale2*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale2*buffer_if
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 5, 51, 52, 53, 54, 8  --  SIGNS= 1
!      Accumulate dGxdt --- derivative wrt k --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.&
&       (choice_==5.or.choice_==51.or.choice_==52.or.choice_==53.or.choice_==54.or.choice_==8)) then
         i1=1;i2=2;i3=3
         if (choice_==54) then
           i1=4;i2=5;i3=6
         end if
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3), &
!$OMP REDUCTION(+:buffer_i1,buffer_i2,buffer_i3)
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
               buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
               buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
               buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_r3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_i3
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_r3
             end if
!$OMP END SINGLE
           else
             cplex_dgxdt(i1:i3) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
!$OMP SINGLE
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1,buffer_i2,buffer_i3)
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
                 buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
                 buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_i3
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_r3)
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
                 buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
                 buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_r3
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 5, 51, 52 -- SIGNS= 2
!      Accumulate dGxdt --- derivative wrt k --- for direction IDIR
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==5.or.choice_==51.or.choice_==52)) then
         i1=1
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_i1)
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
!$OMP END SINGLE
           else
             cplex_dgxdt(i1) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
!$OMP SINGLE
               buffer_i1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1)
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_r1 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1)
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE= 53 -- SIGNS= 2
!      Accumulate dGxdt --- derivative wrt k --- for directions IDIR-1 & IDIR+1
!      --------------------------------------------------------------------
       if ((signs==2).and.(choice_==53)) then
         i1=1;i2=2
!        Case 53 though need multiple directions and here ffnl will contain the
!        derivative information in locations 2,3, and 4 corresponding to idir = 1, 2,
!        and 3. Moreover, choice 53 needs the derivatives in direction idir+1 and idir-1.
!        The parameter vector ffnl_dir_dat contains the necessary translations in
!        locations 1,2 for idir=1; 3,4 for idir=2; and 5,6 for idir=3.
         ffnl_dir1 = ffnl_dir_dat(2*idir-1) ! idir+1 derivative
         ffnl_dir2 = ffnl_dir_dat(2*idir)   ! idir-1 derivative
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_r1 = zero ; buffer_i1 = zero
             buffer_r2 = zero ; buffer_i2 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2,buffer_i1,buffer_i2)
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
               buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
             end do ! loop over npw
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_i2
             else
               dgxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_r2
             end if ! end if on parity
!$OMP END SINGLE
           else
             cplex_dgxdt(i1:i2) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
!$OMP SINGLE
               buffer_i1 = zero ; buffer_i2 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_i1,buffer_i2)
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
                 buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_i2
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_r1 = zero ; buffer_r2 = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_r1,buffer_r2)
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir1,ilmn)
                 buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,ffnl_dir2,ilmn)
               end do ! loop over npw
!$OMP END DO
!$OMP SINGLE
               dgxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_r2
!$OMP END SINGLE
             end if !end if on parity
           end if ! end if on cplx
         end do ! end loop over nlmn
       end if

!      --------------------------------------------------------------------
!      CHOICE= 54,  --  SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt k and atm. pos --- in all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==54)) then
         ishift=0
!$OMP DO &
!$OMP PRIVATE(mu,mua,mub,il,scale,scale2,parity), &
!$OMP PRIVATE(buffer_r,buffer_i)
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 mu=ishift+mub+3*(mua-1)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                 end if
               end do
             end do
           else
             cplex_d2gxdt(ishift+1:ishift+9) = 2 ! Warning dgxdt is here pure imaginary
             if (parity) then
               do mua=1,3 ! atm. pos
                 do mub=1,3 ! k
                   mu=ishift+mub+3*(mua-1)
                   buffer_r = zero
                   do ipw=1,npw
                     buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             else
               do mua=1,3 ! atm. pos
                 do mub=1,3 ! k
                   mu=ishift+mub+3*(mua-1)
                   buffer_i = zero
                   do ipw=1,npw
                     buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1+mub,ilmn)*kpg(ipw,mua)
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = -scale2*buffer_i
                 end do
               end do
             end if
           end if
         end do
!$OMP END DO
       end if

!      --------------------------------------------------------------------
!      CHOICE 6:
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 strains --- for all directions
!      Accumulate d2Gxdt --- 2nd derivative wrt strain and atm. pos --- for all directions
!        --------------------------------------------------------------------
       if ((signs==1).and.(choice_==6)) then
!$OMP SECTIONS
!$OMP SECTION
         ABI_MALLOC(scalarr,(npw,10))
!$OMP SECTION
         ABI_MALLOC(scalari,(npw,10))
!$OMP END SECTIONS
!$OMP DO &
!$OMP PRIVATE(mu,mua,mub,ishift,nua1,nua2,nub1,nub2), &
!$OMP PRIVATE(il,scale,scale2,parity), &
!$OMP PRIVATE(gama,gamb,gamc,gamd), &
!$OMP PRIVATE(kpga,kpgb,kpgc,kpgd), &
!$OMP PRIVATE(buffer_r,buffer_i)
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           do mu=1,10
             do ipw=1,npw
               scalarr(ipw,mu)=scalr(ipw)*ffnl(ipw,mu,ilmn)
               scalari(ipw,mu)=scali(ipw)*ffnl(ipw,mu,ilmn)
             end do
           end do

!          ===== Accumulate 2nd derivative of Gx wrt two strains =====
           ishift=0;scale2=quarter*scale
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,6
                 mu=mua+6*(mub-1)
                 nua1=alpha(mua);nua2=beta(mua)
                 gama=gamma(nub1,nua1)
                 gamb=gamma(nub2,nua1)
                 gamc=gamma(nub1,nua2)
                 gamd=gamma(nub2,nua2)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   kpga=kpg(ipw,nub1)
                   kpgb=kpg(ipw,nub2)
                   kpgc=kpg(ipw,nua1)
                   kpgd=kpg(ipw,nua2)
                   buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                   + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                   + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_r = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                     + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_i = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                     + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale2*buffer_i
                 end do
               end do
             end if
           end if

!          ===== Accumulate 2nd derivative of Gx wrt strain and atm pos. =====
           ishift=36;scale2=pi*scale
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,3
                 mu=ishift+mua+3*(mub-1)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   buffer_i = buffer_i + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalari(ipw,1+nub2)*kpg(ipw,nub1))
                   buffer_r = buffer_r + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) =-scale2*buffer_r
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_i = zero
                   do ipw=1,npw
                     buffer_i = buffer_i + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalari(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_i
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_r = zero
                   do ipw=1,npw
                     buffer_r = buffer_r + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale2*buffer_r
                 end do
               end do
             end if
           end if
         end do
!$OMP END DO
!$OMP SECTIONS
!$OMP SECTION
         ABI_FREE(scalarr)
!$OMP SECTION
         ABI_FREE(scalari)
!$OMP END SECTIONS
       end if

!      --------------------------------------------------------------------
!      CHOICE= 8  --  SIGNS= 1
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 k wave vect. --- for all directions
!      --------------------------------------------------------------------
       if ((signs==1).and.(choice_==8)) then
         i1=1;i2=2;i3=3;i4=4;i5=5;i6=6
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
!$OMP SINGLE
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_ra,buffer_rb,buffer_rc), &
!$OMP REDUCTION(+:buffer_rd,buffer_re,buffer_rf), &
!$OMP REDUCTION(+:buffer_ia,buffer_ib,buffer_ic), &
!$OMP REDUCTION(+:buffer_id,buffer_ie,buffer_if)
             do ipw=1,npw
               buffer_ra = buffer_ra + scalr(ipw)*ffnl(ipw,5,ilmn)
               buffer_rb = buffer_rb + scalr(ipw)*ffnl(ipw,6,ilmn)
               buffer_rc = buffer_rc + scalr(ipw)*ffnl(ipw,7,ilmn)
               buffer_rd = buffer_rd + scalr(ipw)*ffnl(ipw,8,ilmn)
               buffer_re = buffer_re + scalr(ipw)*ffnl(ipw,9,ilmn)
               buffer_rf = buffer_rf + scalr(ipw)*ffnl(ipw,10,ilmn)
               buffer_ia = buffer_ia + scali(ipw)*ffnl(ipw,5,ilmn)
               buffer_ib = buffer_ib + scali(ipw)*ffnl(ipw,6,ilmn)
               buffer_ic = buffer_ic + scali(ipw)*ffnl(ipw,7,ilmn)
               buffer_id = buffer_id + scali(ipw)*ffnl(ipw,8,ilmn)
               buffer_ie = buffer_ie + scali(ipw)*ffnl(ipw,9,ilmn)
               buffer_if = buffer_if + scali(ipw)*ffnl(ipw,10,ilmn)
             end do
!$OMP END DO
!$OMP SINGLE
             if (parity) then
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale*buffer_rf
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale*buffer_if
             else
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_ia
               d2gxdt(2,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_ib
               d2gxdt(2,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_ic
               d2gxdt(2,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale*buffer_id
               d2gxdt(2,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale*buffer_ie
               d2gxdt(2,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale*buffer_if
               d2gxdt(2,i6,ilmn,ia,ispinor) = scale*buffer_rf
             end if
!$OMP END SINGLE
           else
             if (parity) then
!$OMP SINGLE
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_ra,buffer_rb,buffer_rc), &
!$OMP REDUCTION(+:buffer_rd,buffer_re,buffer_rf)
               do ipw=1,npw
                 buffer_ra = buffer_ra + scalr(ipw)*ffnl(ipw,5,ilmn)
                 buffer_rb = buffer_rb + scalr(ipw)*ffnl(ipw,6,ilmn)
                 buffer_rc = buffer_rc + scalr(ipw)*ffnl(ipw,7,ilmn)
                 buffer_rd = buffer_rd + scalr(ipw)*ffnl(ipw,8,ilmn)
                 buffer_re = buffer_re + scalr(ipw)*ffnl(ipw,9,ilmn)
                 buffer_rf = buffer_rf + scalr(ipw)*ffnl(ipw,10,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               d2gxdt(1,i1,ilmn,ia,ispinor) = scale*buffer_ra
               d2gxdt(1,i2,ilmn,ia,ispinor) = scale*buffer_rb
               d2gxdt(1,i3,ilmn,ia,ispinor) = scale*buffer_rc
               d2gxdt(1,i4,ilmn,ia,ispinor) = scale*buffer_rd
               d2gxdt(1,i5,ilmn,ia,ispinor) = scale*buffer_re
               d2gxdt(1,i6,ilmn,ia,ispinor) = scale*buffer_rf
!$OMP END SINGLE
             else
!$OMP SINGLE
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
!$OMP END SINGLE
!$OMP DO &
!$OMP REDUCTION(+:buffer_ia,buffer_ib,buffer_ic), &
!$OMP REDUCTION(+:buffer_id,buffer_ie,buffer_if)
               do ipw=1,npw
                 buffer_ia = buffer_ia + scali(ipw)*ffnl(ipw,5,ilmn)
                 buffer_ib = buffer_ib + scali(ipw)*ffnl(ipw,6,ilmn)
                 buffer_ic = buffer_ic + scali(ipw)*ffnl(ipw,7,ilmn)
                 buffer_id = buffer_id + scali(ipw)*ffnl(ipw,8,ilmn)
                 buffer_ie = buffer_ie + scali(ipw)*ffnl(ipw,9,ilmn)
                 buffer_if = buffer_if + scali(ipw)*ffnl(ipw,10,ilmn)
               end do
!$OMP END DO
!$OMP SINGLE
               d2gxdt(1,i1,ilmn,ia,ispinor) =-scale*buffer_ia
               d2gxdt(1,i2,ilmn,ia,ispinor) =-scale*buffer_ib
               d2gxdt(1,i3,ilmn,ia,ispinor) =-scale*buffer_ic
               d2gxdt(1,i4,ilmn,ia,ispinor) =-scale*buffer_id
               d2gxdt(1,i5,ilmn,ia,ispinor) =-scale*buffer_ie
               d2gxdt(1,i6,ilmn,ia,ispinor) =-scale*buffer_if
!$OMP END SINGLE
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      END CHOICES
!      --------------------------------------------------------------------

     end do ! End loop on atoms

!    Deallocate temporary space
!$OMP SECTIONS
!$OMP SECTION
     ABI_FREE(scali)
!$OMP SECTION
     ABI_FREE(scalr)
!$OMP END SECTIONS

   end do !  End loop on spinorial components
!$OMP END PARALLEL

!  ==========================================================================
 end if

!Has to reduce arrays in case of FFT parallelization
 if (mpi_enreg%nproc_fft>1) then
   call timab(48,1,tsec)
   if (choice>=0) then
     call xmpi_sum(gx,mpi_enreg%comm_fft,ierr)
   end if
   if (choice_>1) then
     call xmpi_sum(dgxdt,mpi_enreg%comm_fft,ierr)
   end if
   if (choice_==4.or.choice_==24.or.choice_==54.or.choice_==6.or.choice==8) then
     call xmpi_sum(d2gxdt,mpi_enreg%comm_fft,ierr)
   end if
   call timab(48,2,tsec)
 end if

end subroutine opernla_ylm_met
!!***

end module m_opernla_ylm_met
!!***
