!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_dfpt_indpol
!! NAME
!! m_dfpt_indpol
!!
!! FUNCTION
!! Tools for computing the induced polarization, e.g., for Flexo
!!
!! Written by Cyrus Dreyer, Rutgers 2016
!!
!!
!! INPUTS
!!  calcnl = 1 to calc NL part, 0 to just to local
!!  cg = GS wavefunction pw coefs
!!  cg1 = AD wavefunction pw coefs
!!  cryst <type(crystal)>=Unit cell and symmetries
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  kg = GS wavefunction G vectors
!!  kg = AD wavefunction G vectors 
!!  psps = Psuedopotential info  
!! 
!!
!!
!!
!!
!! OUTPUT
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

module m_dfpt_indpol

 use m_cgtools
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_dtset
 use m_dtfil
! use m_crystal,    only : crystal_init, crystal_free, crystal_t,isalchemical
 use m_pawcprj, only : pawcprj_type
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_paw_sphharm!, only : ass_leg_pol, plm_dtheta, plm_dphi, plm_coeff
 use m_time,       only : timab
 use m_fftcore,    only : sphereboundary
 use m_fft
 use m_splines
 use m_mpinfo
 use m_geometry

implicit none

private

public :: indpol ! Main subroutine for calculating induced polarization
public :: joper ! Apply current operator up to second order in q
public :: cmkffnl ! Fully cartesian and similified version of mkffnl, TODO: use original instead 
public :: cinitylmgi ! Fully cartsian and simplified version of initylmgi, TODO: use original instead

contains

!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_indpol
!! NAME
!! dfpt_indpol
!!
!! FUNCTION
!! Compute the induced polarization 
!!
!! Written by Cyrus Dreyer, Rutgers 2016
!!
!!
!! INPUTS
!!  calcnl = 1 to calc NL part, 0 to just to local
!!  cg = GS wavefunction pw coefs
!!  cg1 = AD wavefunction pw coefs
!!  cryst <type(crystal)>=Unit cell and symmetries
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  kg = GS wavefunction G vectors
!!  kg = AD wavefunction G vectors 
!!  psps = Psuedopotential info  
!! 
!!
!!
!!
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine indpol(calcden,calcnl,cg,cg1,cgp,cgq,dtfil,dtset,gmet,gprimd,idir,ipert,istwfk,kg,kg1,kpt, &
     & mcg,mcgq,mcg1,mk1mem_rbz,mkmem_rbz,mpi_enreg,mpw,mpw1,nband,nfftf,npwarr,npwar1,& 
     & nkpt,occ,psps,rhorout,rhonlrout,rmet,rprimd,ucvol,wtk)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indpol'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in), target :: dtset
 type(pseudopotential_type), intent(inout) :: psps
 type(MPI_type), intent(in) :: mpi_enreg
 integer, intent(in) :: mcg,mcg1,mcgq,mkmem_rbz,mk1mem_rbz,mpw,mpw1
 integer, intent(in) :: calcnl,calcden
 integer, intent(in) :: nkpt,nfftf
 integer, intent(in) :: istwfk(nkpt),nband(nkpt),ipert,idir
 integer, intent(in) :: kg(3,mpw*mkmem_rbz),kg1(3,mpw1*mkmem_rbz)
 integer, intent(in) :: npwarr(nkpt),npwar1(nkpt)
 real(dp), intent(in) :: cg(2,mcg), cg1(2,mcg1),cgp(2,mcgq),cgq(2,mcgq)
 real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp), intent(in) :: kpt(3,nkpt),wtk(nkpt)
 real(dp), intent(inout) :: rhonlrout(2*nfftf)
 real(dp), intent(inout) :: rhorout(2*nfftf)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),ucvol

!Local variables--------------------------------
 integer, parameter :: tim_fourwf0=0 ! timing for FFT
 integer :: mm,kk,jj,m1,k1,j1,iwrite,calcpm ! various loops
 integer :: ikpt,ibnd,ikg,ix,ir,iat,ilm,ix1,ix2,id,id1,il,ityp,fid,fid1,itypat ! various loops
 integer :: placek, placek1,placebnd,placebnd1,placept,placeout! To keep the place in input files
 integer :: nbnd,npw,npw1,mpw2,nln ! numbers of things
 integer :: ngfft(18),n4,n5,n6 ! dimensions for fft
 integer :: mgfft
 integer :: cplex
 real(dp) :: ecut_eff
 real(dp) :: weight ! scrap for fft
! arrays
 integer, allocatable :: gbd(:,:),gbd1(:,:) ! for FFT, largest g
! integer, allocatable :: npwarr(:),npwar1(:)
 integer, allocatable :: kg_k(:,:),kg1_k(:,:)
 integer, allocatable :: nlmn(:),iln(:,:)
! real(dp) :: rprimd(3,3) ! lattice vectors
! real(dp) :: kpt(3)
 real(dp) :: xcart(3,dtset%natom),qpc(3) !,gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: kplusg(3) ! G vectors, and shifted g vectors points
 real(dp) :: cpuin,wallin,cpuout,wallout
 real(dp), allocatable :: gsg(:,:),adg(:,:),gsgd(:,:,:),adgd(:,:,:), adgs(:,:)!gsgdq(:,:,:) ! reciprocal wfks and derivatives for given band 
 real(dp), allocatable :: jgsg(:,:,:),jgsg1(:,:),jgsg1nl(:,:),jgsgnl(:,:,:),jgsg1_tmp(:,:),psitil(:,:)
 real(dp), allocatable :: denpot(:,:,:) ! scap for fft    
 real(dp), allocatable :: pindr(:,:,:,:,:),pindg(:,:,:),vnlr(:,:,:,:,:),vnlg(:,:,:) ! output induced p
 real(dp), allocatable :: rhor(:,:,:,:),rhonlr(:,:,:,:) ! output density
 real(dp), allocatable :: fog(:,:)
 real(dp), allocatable :: dumg(:,:),dumg1(:,:),dumr(:,:,:,:,:)
 real(dp), allocatable :: kpq(:,:) ! for finite q
 real(dp), allocatable :: scprod(:,:)
 complex(dp) :: jpindg(3),psitil0(3)
 !complex(dp), allocatable :: ffnl(:,:,:,:),ffnl1(:,:,:,:)
 complex(dp) :: vnlg0(3),pindg0(3),pos,jpindg0(3),jvnlg0(3)
 complex(dp) :: fcta1,fcta2(3),fcta3(3),fcta4,fcta5(3,3),fcta6(3,3),fcta7(3,3,3),fcta8(3,3,3)
 complex(dp) :: rhoqnl,rhoqnl3(3)
 integer,allocatable :: sign_dyad(:,:)

 ! TEST 08/06/19
 real(dp) :: kcart(3),kmod


!************************************************
! Some inital setup
!************************************************

!VER9rev: timein to timeab
!timing
!call timeab(cpuin,wallin)

! TEST
 write(*,*) "In indpol"

! Some inital variables we will use a lot
nbnd=nband(1) ! no metals
mgfft=dtset%mgfft
ngfft=dtset%ngfft

ecut_eff=dtset%ecut*(dtset%dilatmx)**2

! calcpm: 0 is ICL path, 1 is PM
calcpm=0
if (dtset%useric>0) calcpm=dtset%useric
if (calcpm>0) write(*,*) 'WARNING: PM PATH IS ON!!!!'

! VER9rev: use dtset for ngfft
! Initializations for FFTs
!call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfft(2),ngfft(3),'all')
n4=dtset%ngfft(1)!ngfft(4)
n5=dtset%ngfft(2)!ngfft(5)
n6=dtset%ngfft(3)!ngfft(6)
cplex=0

!TEST
!write(*,*) "proc",mpi_enreg%me_kpt, "nkpt",nkpt,"mkmem_rbz",mkmem_rbz,"mk1mem_rbz",mk1mem_rbz
!write(*,*) "mpi_enreg%distribfft indpol 151",mpi_enreg%distribfft


! VER9rev: Now passed so don't need this
!Obtain dimensional translations in reciprocal space gprimd,
!metrics and unit cell volume, from rprimd. Also output rprimd, gprimd and ucvol
!call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
!call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!TEST:
!write(*,*) "GPRIMD"
!write(*,*) gprimd(1,:)
!write(*,*) gprimd(2,:)
!write(*,*) gprimd(3,:)
!stop

!write(*,*) "GPRIMD IS TRANSPOSED!!!!!!"

!Cartesian q: 
qpc(:)=two_pi*matmul(gprimd(:,:),dtset%qptn(:))

!Get catesian positions of atoms for PM nonlocal part
do iat=1,dtset%natom
   xcart(:,iat)=matmul(rprimd(:,:),dtset%xred_orig(:,iat,1))
end do

!Phase shift for atoms not at the origin
!TEST:
pos=cmplx(cos(dot_product(xcart(:,ipert),qpc(:))), &
     &   sin(dot_product(xcart(:,ipert),qpc(:)))) 

!pos=cmplx(cos(dot_product(xcart(:,dtset%rfatpol(1)),qpc(:))), &
!     &   sin(dot_product(xcart(:,dtset%rfatpol(1)),qpc(:))))
!if (dtset%rfatpol(1) /= dtset%rfatpol(2)) then
!   write(*,*) 'Can only handle one atom motion per dataset!!!'
!   stop
!end if

!************************************************
! Setup the reciprocal space grids
!************************************************
!ABI_MALLOC(npwarr,(nkpt))
!ABI_MALLOC(npwar1,(nkpt))
ABI_MALLOC(gbd,(2*mgfft+8,2))
ABI_MALLOC(gbd1,(2*mgfft+8,2))

! AB9rev: passed mpw
!  Compute maximum number of planewaves at k
!call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk,kpt,mpi_enreg,mpw,nkpt)


! AB9rev: passed kg
!ABI_MALLOC(kg,(3,mpw*mkmem_rbz))
!ABI_MALLOC(kg_k,(3,mpw))

! get G vectors for gs
!call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk,kg,dtfil%fnametmp_kg,&
!     &   kpt,mkmem_rbz,nband,nkpt,'PERS',mpi_enreg,mpw,npwarr,npwarr,dtset%nsppol,dtfil%unkg)

! k+q mesh
ABI_MALLOC(kpq,(3,nkpt))
do ikpt=1,nkpt
   kpq(:,ikpt)=dtset%qptn(:)+kpt(:,ikpt)
end do

! AB9rev: Now all passed
!  Compute maximum number of planewaves at k+q
!call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk,kpq,mpi_enreg,mpw1,nkpt)
!ABI_MALLOC(kg1,(3,mpw1*mk1mem_rbz))
! get G vectors for ad wfk at k+q                       
!call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk,kg1,dtfil%fnametmp_kg1,&
!     &   kpq,mkmem_rbz,nband,nkpt,'PERS',mpi_enreg,mpw1,npwar1,npwar1,dtset%nsppol,dtfil%unkg1)

!*****************************************************
! allocate arrays
!*****************************************************
! real space, allocate here since they don't change with kpt
if (calcden==1) then
   ABI_MALLOC(pindr,(3,2,n4,n5,n6))
   ABI_MALLOC(rhor,(2,n4,n5,n6))
   pindr=0;rhor=0;
end if

! dummy
ABI_MALLOC(dumr,(3,2,n4,n5,n6))
!TEST: for now allocate w/ cplex=1
!ABI_MALLOC(denpot,(cplex*n4,n5,n6))
ABI_MALLOC(denpot,(n4,n5,n6)) 
denpot=one
! nonlocal 
if (calcnl == 1) then
   ABI_MALLOC(vnlr,(3,2,n4,n5,n6))
   ABI_MALLOC(rhonlr,(2,n4,n5,n6))
   vnlr=0;rhonlr=0
   vnlg0=zero;jpindg=zero;pindg0=zero;rhoqnl=zero
   jpindg0=zero;jvnlg0=zero;psitil0=zero
   ABI_MALLOC(sign_dyad,(psps%lnmax,psps%ntypat))
end if


! TEST:
!do ikpt=1,nkpt
!   write(*,*) ikpt,mpi_enreg%me_kpt,proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nbnd,1,mpi_enreg%me_kpt)
!end do
!stop

!******************************************************
! Big k point loop
!******************************************************

placek=0
placek1=0
placebnd=0
placebnd1=0
placept=0
do ikpt=1,nkpt

   write(*,*) "KPT",ikpt,"/",nkpt

   ! Test if this kpoint should be treated by this process
   ! no spin, no metal
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nbnd,1,mpi_enreg%me_kpt)) cycle

   !TEST
   !write(*,*) 'ikpt',ikpt 

   npw=npwarr(ikpt)
   npw1=npwar1(ikpt)

   ! Allocate arrays
   ABI_MALLOC(gsg,(2,npw))
   ABI_MALLOC(adg,(2,npw1))
   ABI_MALLOC(fog,(2,npw1))
   ABI_MALLOC(gsgd,(3,2,npw))
   ABI_MALLOC(adgd,(3,2,npw1))
   ABI_MALLOC(dumg,(2,npw))
   ABI_MALLOC(kg_k,(3,npw))
   ABI_MALLOC(kg1_k,(3,npw1))
   If (calcnl==1) then
      ABI_MALLOC(adgs,(2,npw))
      !ABI_MALLOC(ffnl,(npw,20,Psps%lnmax*Psps%lnmax,dtset%natom))
      !ABI_MALLOC(ffnl1,(npw1,20,Psps%lnmax*Psps%lnmax,dtset%natom))
      !Joper implemantation
      ABI_MALLOC(jgsg,(3,2,npw))
      ABI_MALLOC(jgsg1,(2,npw1))
      ABI_MALLOC(jgsg1nl,(2,npw1))
      ABI_MALLOC(jgsgnl,(3,2,npw))
      jgsg=zero;jgsg1=zero;jgsg1nl=zero;jgsgnl=zero
      !Psitilde
      ABI_MALLOC(psitil,(2,npw1))
   end if

   ! TEST
   !write(*,'(a10,2i10,a20,2i10,a10,4i10)') "npw,npw1",npw,npw1,"placek,placek1",placek,placek1,"max,max1",mpw,mpw1,mkmem_rbz,mk1mem_rbz

   do ikg=1,npw
      kg_k(:,ikg)=kg(:,ikg+placek)
   end do

   do ikg=1,npw1
      kg1_k(:,ikg)=kg1(:,ikg+placek1)
   end do

   placek=placek+npw
   placek1=placek1+npw1

   !TEST:
   if (npw==0 .or. npw1==0) then
      write(*,*) 'npw=0, something is wrong...'
      stop
   end if
      
   ! gbounds for all wfks
   call sphereboundary(gbd,istwfk(ikpt),kg_k,mgfft,npw)
   call sphereboundary(gbd1,istwfk(ikpt),kg1_k,mgfft,npw1)


   ! Initialize KB potential. Have this here since npw may be different for different kpoints
   !if (calcnl==1) then
   !   ffnl=zero;ffnl1=zero;
   !   call cmkffnlindpol(0,dtfil,dtset,ffnl,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,sign_dyad)
   !   call cmkffnlindpol(0,dtfil,dtset,ffnl1,kg1_k,kpt(:,ikpt),mpi_enreg,npw1,psps,sign_dyad)
   !end if

   !loop over bands
   do ibnd=1,nbnd

      !TEST
      !write(*,*) 'ibnd',ibnd

      !make sure that this band should be treated by this processor 
      if (mpi_enreg%proc_distrb(ikpt,ibnd,1) /= mpi_enreg%me_kpt) then
         placebnd=placebnd+npw
         placebnd1=placebnd1+npw1
         cycle
      end if
      
      ! Only occupied bands
      if (occ(nbnd*(ikpt-1)+ibnd)<1d-2) then
         write(*,*) 'done with occupied bands'
         exit
      end if

      do ikg=1,npw
         gsg(:,ikg)=cg(:,ikg+placebnd)
         ! k+G in cartesian coordinates. 
         kplusg(:) = two_pi*matmul(gprimd(:,:),kpt(:,ikpt)+real(kg_k(:,ikg)))

         gsgd(1,1,ikg)=kplusg(1)*gsg(1,ikg)
         gsgd(1,2,ikg)=kplusg(1)*gsg(2,ikg)
         gsgd(2,1,ikg)=kplusg(2)*gsg(1,ikg)
         gsgd(2,2,ikg)=kplusg(2)*gsg(2,ikg)
         gsgd(3,1,ikg)=kplusg(3)*gsg(1,ikg)
         gsgd(3,2,ikg)=kplusg(3)*gsg(2,ikg)
      end do
      placebnd=placebnd+npw

      do ikg=1,npw1
         adg(:,ikg)=cg1(:,ikg+placebnd1)         
         kplusg(:) = two_pi*matmul(gprimd(:,:),kpt(:,ikpt)+dtset%qptn(:)+real(kg1_k(:,ikg)))

         adgd(1,1,ikg)=kplusg(1)*adg(1,ikg)
         adgd(1,2,ikg)=kplusg(1)*adg(2,ikg)
         adgd(2,1,ikg)=kplusg(2)*adg(1,ikg)
         adgd(2,2,ikg)=kplusg(2)*adg(2,ikg)
         adgd(3,1,ikg)=kplusg(3)*adg(1,ikg)
         adgd(3,2,ikg)=kplusg(3)*adg(2,ikg)
         
         ! For FO chg density
         fog(:,ikg)=cgp(:,ikg+placebnd1)
      end do
      placebnd1=placebnd1+npw1
      
      !********************************************************
      ! Compute the current: joper implementation             *        
      !********************************************************
      if (dtset%userib==1) then
         !do ix=1,3
                     
            !Convert mesh
            jgsg1(:,:)=adg(:,:)
            denpot=one
            call fourwf(1,denpot,jgsg1,jgsg(1,:,:),dumr(1,:,:,:,:),gbd1,gbd,&
                 &     istwfk(ikpt),kg1_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw1,&
                 &     npw,n4,n5,n6,2,0,one,one)
            jgsgnl(:,:,:)=jgsg(:,:,:)
            !Apply current operator
            !Just local part
            !call joper_sign(0,jgsg,dtfil,dtset,ffnl,gprimd,ix,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc,rprimd,sign_dyad)
            call joper(0,jgsg,dtfil,dtset,gprimd,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc)
            !Just NL part
            !if (calcnl==1) call joper_sign(2,jgsgnl,dtfil,dtset,ffnl,gprimd,ix,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc,rprimd,sign_dyad)
            if (calcnl==1) call joper(2,jgsgnl,dtfil,dtset,gprimd,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc)
            do ikg=1,npw
               jpindg0(:)=jpindg0(:)+wtk(ikpt)*conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*cmplx(jgsg(:,1,ikg),jgsg(:,2,ikg))
               jvnlg0(:)=jvnlg0(:)+wtk(ikpt)*conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*cmplx(jgsgnl(:,1,ikg),jgsgnl(:,2,ikg))
            end do
            !  end if
         !end do
         !TEST
         write(*,'(a5,2i4,12e12.4e2)') 'IND:',ibnd,ikpt,jpindg0,jvnlg0
      end if
      !********************************************************
      ! Diamagnetic Susceptability                            *
      !********************************************************
      if (dtset%userib==2) then

         !do ix=1,3
            jgsg1(:,:)=adg(:,:)

            ! Need to convert adg to kg_k mesh
            denpot=one
            call fourwf(1,denpot,jgsg1,jgsg(1,:,:),dumr(1,:,:,:,:),gbd1,gbd,&
                 &     istwfk(ikpt),kg1_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw1,&
                 &     npw,n4,n5,n6,2,0,one,one)
            !Apply current operator 
            !call joper_sign(1,jgsg,dtfil,dtset,ffnl,gprimd,ix,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc,rprimd,sign_dyad)
            call joper(1,jgsg,dtfil,dtset,gprimd,kg_k,kpt(:,ikpt),mpi_enreg,npw,psps,qpc)
            do ikg=1,npw
               !Factor of 1/4 due to the fact that there is an extra factor of 4 in J from solving Sternheimer Equation
               jpindg(:)=jpindg(:)+0.25*wtk(ikpt)*conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*cmplx(jgsg(:,1,ikg),jgsg(:,2,ikg))
               !jpindg(ix)=jpindg(ix)+4*wtk(ikpt)*conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*cmplx(jgsg(1,ikg),jgsg(2,ikg))
            end do
         !end do !ix
      end if !userib==2
      

      !*************************************************************
      ! Compute the current densities, and G=0 using old method    *
      !*************************************************************

      ! AB9rev: Not sure whether to keep this. I think it is useful for q=0 local polarization and testing, TODO
      !if (calcden==1 .and. dtset%userib == 1) then
      !   call polden(adg,adgd,calcnl,dtfil,dtset,gbd,gbd1,gsgd,gsg,istwfk(ikpt),kg_k,kg1_k,kpt(:,ikpt), &
      !        & mpi_enreg,ngfft,nkpt,npw,npw1,pindg0,pindr,psps,qpc,rhoqnl,ucvol,vnlg0,vnlr,rhonlr,wtk(ikpt)) 
      !end if

   end do !band

   !Deallocate arrays
   ABI_FREE(gsg)
   ABI_FREE(adg)
   ABI_FREE(fog)
   ABI_FREE(gsgd)
   ABI_FREE(adgd)
   ABI_FREE(dumg)
   ABI_FREE(kg_k)
   ABI_FREE(kg1_k)
   If (calcnl==1) then
      ABI_FREE(adgs)
      !ABI_FREE(ffnl)
      !ABI_FREE(ffnl1)
      ABI_FREE(jgsg)
      ABI_FREE(jgsg1)
      ABI_FREE(jgsg1nl)
      ABI_FREE(jgsgnl)
      ABI_FREE(psitil)
   end if

   !Need place counter for PsiTilde that is start for a given kpoint
   placept=placept+npw1*nbnd

end do !kpt
 
!***************************************************************
! Write out quantities
!***************************************************************
jvnlg0(:)=matmul(jvnlg0(:),gprimd(:,:))
jpindg0(:)=matmul(jpindg0(:),gprimd(:,:))


!Phase shift for atoms not at the origin (phonons only)
if (dtset%userib==1.and.dtset%userid==0) then
   jpindg0(:)=pos*jpindg0(:) !Joper imp
   jvnlg0(:)=pos*jvnlg0(:)
   if (mpi_enreg%me_kpt==0) write(*,*) "e^iq\tau phase added"
end if

! TEST
write(*,*) "Before summation"

!MPI case:
!#ifdef HAVE_MPI
call xmpi_barrier(mpi_enreg%comm_kpt)

! TEST
write(*,*) "Before sum_master"

call xmpi_sum_master(jpindg0,0,mpi_enreg%comm_kpt,jj)
if (jj==1) write(*,*) 'ERROR: Summation of jpindg0'
call xmpi_sum_master(jvnlg0,0,mpi_enreg%comm_kpt,jj)
if (jj==1) write(*,*) 'ERROR: Summation of jvnlg0'

if (dtset%userib==2) then
   call xmpi_sum_master(jpindg,0,mpi_enreg%comm_kpt,jj)
   if (jj==1) write(*,*) 'ERROR: Summation of jpindg'
end if

! TEST
write(*,*) "Before writeout"


if (mpi_enreg%me_kpt==0) then
   
   !AB9rev: For now, lets output to stout
   !open (unit=21, file='pindg.dat', position = 'append',Status='old')
   !open (unit=22, file='vnlg.dat', position = 'append',Status='old')
   !open (unit=23, file='totg0.dat', position = 'append',Status='old')

!   write(21,'(8e20.10e2)') jpindg0(:)
   write(*,'(a10,2i5,3e16.6e2,8e20.10e2)') "PINDG",ipert,idir,qpc(:), jpindg0(:)

   if (dtset%userib==2) then
      write(*,'(a10,2i5,3e16.6e2,6e20.10e2)') "DIAMAG",ipert,idir,qpc(:),jpindg(:) !'Diamag Susc,is :',jpindg(:)
   else
      !write(23,'(8e20.10e2)') jpindg0(:)+jvnlg0(:)
      write(*,'(a10,2i5,3e16.6e2,8e20.10e2)') "TOTG",ipert,idir,qpc(:),jpindg0(:)+jvnlg0(:)
   end if
   !write(22,'(a30,6e20.10e2)') 'nonlocal:',jvnlg0(:)
   write(*,'(a10,2i5,3e16.6e2,8e20.10e2)') "VNLG",ipert,idir,qpc(:), jvnlg0(:)
end if
   !close (unit=21)
   !close (unit=22)
   !close (unit=23)
!TEST
!write(*,*) 'after close'

!iwrite=1
!#endif

! Sequential case
!!$if (iwrite==0) then
!!$   if (dtset%userib==3) then
!!$      write(*,'(a15,6e16.5e2)') 'PsiTilde',psitil0(:)
!!$   end if
!!$
!!$   open (unit=21, file='pindg.dat', position = 'append',Status='old')
!!$   open (unit=22, file='vnlg.dat', position = 'append',Status='old')
!!$   open (unit=23, file='totg0.dat', position = 'append',Status='old')
!!$   write(21,'(8e20.10e2)') jpindg0(:)
!!$
!!$   if (dtset%userib==2) then
!!$      write(23,'(a30,6e20.10e2)') 'Diamag Susc is :',jpindg(:)
!!$   else
!!$      write(23,'(8e20.10e2)') jpindg0(:)+jvnlg0(:)
!!$   end if
!!$   write(22,'(a30,6e20.10e2)') 'nonlocal',jvnlg0(:)
!!$end if
!!$close (unit=21)
!!$close (unit=22)
!!$close (unit=23)



! Write densities and G=0 using polden method

! Skip this for now
!write(*,*) "mpi_enreg%distribfft indpol 606",mpi_enreg%distribfft

!if (calcden==1 .and. dtset%userib == 1) then
    
   
   !TEST
   !rhoqnl3(:)= zero; rhoqnl3(1)=rhoqnl
   !call xmpi_sum_master(rhoqnl3,0,mpi_enreg%comm_kpt,jj)
   !if (jj==1) write(*,*) 'ERROR: Summation of rhoqnl'
   !! For now, just one processor
   !if (mpi_enreg%me_kpt==0) write(*,'(a30,2e20.10e2)') 'rhoqnl:',rhoqnl3(1)

   
   !call poldenout(dtfil,dtset,ecut_eff,gmet,gprimd,istwfk,kpt,mkmem_rbz,mpi_enreg,nband, &
   !     & nfftf,ngfft,nkpt,pindg0,pindr,pos,qpc,rhorout,rprimd,ucvol,vnlg0,vnlr,rhonlr,rhonlrout)
!end if



!timing
!call timein(cpuout,wallout)
!write(*,*) 'wall time in indpol:',wallout-wallin


end subroutine indpol
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/joper
!! NAME
!! joper 
!!
!! FUNCTION
!! Applies the local and nonlocal expansion up to 2nd order  of the 
!! current density operator to a wavefunction
!!
!! Written by Cyrus Dreyer, Rutgers 2017
!!
!!
!! INPUTS
!!  calcnl = 1 to calc NL part, 0 to just to local
!!  cg = GS wavefunction pw coefs
!!  cg1 = AD wavefunction pw coefs
!!  cryst <type(crystal)>=Unit cell and symmetries
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  kg = GS wavefunction G vectors
!!  kg = AD wavefunction G vectors 
!!  psps = Psuedopotential info  
!! 
!!
!!
!!
!!
!! OUTPUT
!!
!! PARENTS
!! 
!!
!!
!! CHILDREN
!!
!! SOURCE

subroutine joper(calcnl,cwave0_npw1,dtfil,dtset,gprimd,kg1_k,kpt,mpi_enreg,npw1,psps,qpc)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'joper'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in), target :: dtset
 type(pseudopotential_type), intent(in) :: psps
 type(MPI_type), intent(in) :: mpi_enreg
! type(crystal_t), intent(in) :: crystal
 integer,intent(in) :: npw1,calcnl
 real(dp), intent(inout) :: cwave0_npw1(3,2,npw1)
 real(dp), intent(in) :: kpt(3),qpc(3),gprimd(3,3)
 integer, intent(in) :: kg1_k(3,npw1)


!Local variables-------------------------------
 real(dp) :: gsg(2,npw1),gsgd(3,2,npw1)
 real(dp) :: kplusg(3),qpc_(3)
 real(dp) :: cpuin,wallin,cpuout,wallout
 complex(dp) :: vnlg0(3,npw1), vnlg1(3,npw1),vnlg2(3,npw1)
 complex(dp) :: fcta1,fcta2(3),fcta3(3),fcta4,fcta5(3,3),fcta6(3,3),fcta7(3,3,3),fcta8(3,3,3)
 complex(dp) :: ffnl1(3),ffnl2(3,3),ffnl3(3,3,3)
 complex(dp) :: ffnl(npw1,20,psps%lnmax*psps%lnmax)
 integer :: sign_dyad(psps%lnmax,psps%ntypat)
 integer, allocatable :: nlmn(:),iln(:,:)
 integer :: ikg,ix,ix1,ix2,iat,ilm,il,ityp,nln,fid,tr!,irfdir
 integer, parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer, parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer, parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer, parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
!------------------------------------------------

!************************************************
! Setup
!************************************************
 
! if (calcnl==0) then
!    write(*,*)'LOCAL joper will be applied in direction',irfdir
! else if (calcnl==2) then
!    write(*,*)'NONLOCAL joper will be applied in direction',irfdir
! else if (calcnl==1) then
!    write(*,*)'DIAMAG joper will be applied in direction',irfdir
! end if

!timing
!call timein(cpuin,wallin)

!TESTs
tr=0 !Only TR systyems right now
!sign_dyadj(:,:)=sign_dyad(:,:) !Test removing dyadic sign
!write(*,*) 'joper WITHOUT SD, really'

!************************************************
! Local part
!************************************************
  
 ! Apply local part of current density operator (momentum) to wavefunction
 do ikg=1,npw1
    gsg(:,ikg)=cwave0_npw1(1,:,ikg)
    
    if (calcnl<2) then
       kplusg(:)=two_pi*matmul(gprimd(:,:),kpt(:)+real(kg1_k(:,ikg)))+0.5*qpc(:)

       gsgd(1,1,ikg)=kplusg(1)*gsg(1,ikg)
       gsgd(1,2,ikg)=kplusg(1)*gsg(2,ikg)
       gsgd(2,1,ikg)=kplusg(2)*gsg(1,ikg)
       gsgd(2,2,ikg)=kplusg(2)*gsg(2,ikg)
       gsgd(3,1,ikg)=kplusg(3)*gsg(1,ikg)
       gsgd(3,2,ikg)=kplusg(3)*gsg(2,ikg)

       !    write(*,'(3i5,7e16.5e2)') kg1_k(:,ikg),gsg(:,ikg),kplusg(:),gsgd(irfdir,:,ikg)
    end if
 end do
 ! stop
 !************************************************
 ! Nonlocal potential
 !************************************************

 if (calcnl>0) then
    
    vnlg0=zero;vnlg1=zero;vnlg2=zero;

    ! determine number of nonlocal projectors to loop over  
    ABI_MALLOC(nlmn,(dtset%natom))
    ABI_MALLOC(iln,(dtset%natom,Psps%lnmax*Psps%lnmax))
    nlmn(:)=0
    do iat=1,dtset%natom
       ityp=dtset%typat(iat)
       nln=count(Psps%indlmn(3,:,ityp)>0)
       do ilm=1,nln
          il=Psps%indlmn(1,ilm,ityp)
          !Need this for sign_dyad
          iln(iat,nlmn(iat)+1:nlmn(iat)+2*il+1)=Psps%indlmn(5,ilm,ityp)
          nlmn(iat)= nlmn(iat)+2*il+1
       end do
    end do

    !TEST
    !write(*,*) 'joper nlmn',nlmn(1)
    !stop

    !loop over atoms and ang momentum channels
    do iat=1,dtset%natom
       ityp=dtset%typat(iat)

       ! Calculate ffnl here. This takes more time, but avoinds storing a 
       ! very large array which is prohibative for large cells
       ! TODO: ALSO CALC FFNL ON THE FLY FOR iln!!!!
       ffnl=zero
       call cmkffnl(0,dtfil,dtset,ffnl,iat,kg1_k,kpt,mpi_enreg,npw1,psps,sign_dyad)

       ! Loop over angular momentum channels
       do ilm=1,nlmn(iat)
 
          ! =====Calculate overlaps=====
          ! For now have them all. Since we have TR, we only need half.
          fcta1=0;fcta2=0;fcta3=0;fcta4=0;fcta5=0;fcta6=0;fcta7=0;fcta8=0

          if (tr==1) then
             do ikg=1,npw1

                ! For order in q >=0
                fcta1=fcta1+conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*ffnl(ikg,1,ilm)
                fcta3(:)=fcta3(:)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                     & *ffnl(ikg,2:4,ilm)

                ! For order in q >=1
                fid=4
                do ix=1,3 ! alpha
                   do ix1=ix,3 ! beta
                      fid=fid+1

                      fcta5(ix,ix1)=fcta5(ix,ix1)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                           & *ffnl(ikg,fid,ilm)
                      fcta5(ix1,ix)= fcta5(ix,ix1)

                   end do
                end do
                ! For order in q >=2
                fid=10
                do ix=1,3
                   do ix1=ix,3
                      do ix2=ix1,3
                         fid=fid+1

                         fcta7(ix,ix1,ix2)=fcta7(ix,ix1,ix2)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                              & *ffnl(ikg,fid,ilm)
                         fcta7(ix,ix2,ix1)=fcta7(ix,ix1,ix2)
                         fcta7(ix1,ix,ix2)=fcta7(ix,ix1,ix2)
                         fcta7(ix1,ix2,ix)=fcta7(ix,ix1,ix2)
                         fcta7(ix2,ix,ix1)=fcta7(ix,ix1,ix2)
                         fcta7(ix2,ix1,ix)=fcta7(ix,ix1,ix2)

                      end do
                   end do
                end do
             end do !ikg
          end if !tr

         do ikg=1,npw1

            ! For order in q >=0
            fcta2(:)=fcta2(:)+conjg(ffnl(ikg,2:4,ilm)) &
                 & *cmplx(gsg(1,ikg),gsg(2,ikg))
            fcta4=fcta4+conjg(ffnl(ikg,1,ilm))*cmplx(gsg(1,ikg),gsg(2,ikg))
            
            !TEST
            !write(*,'(a10,4e20.10)') "FTA4",conjg(ffnl(ikg,1,ilm)),gsg(1,ikg),gsg(2,ikg)

            ! For order in q >=1
            fid=4
            do ix=1,3 ! alpha
               do ix1=ix,3 ! beta
                  fid=fid+1

                  fcta6(ix,ix1)=fcta6(ix,ix1)+conjg(ffnl(ikg,fid,ilm)) &
                       & *cmplx(gsg(1,ikg),gsg(2,ikg))
                  fcta6(ix1,ix)=fcta6(ix,ix1)

               end do
            end do
            ! For order in q >=2
            fid=10
            do ix=1,3
               do ix1=ix,3
                  do ix2=ix1,3
                     fid=fid+1

                     fcta8(ix,ix1,ix2)=fcta8(ix,ix1,ix2)+conjg(ffnl(ikg,fid,ilm)) &
                          & *cmplx(gsg(1,ikg),gsg(2,ikg))
                     fcta8(ix,ix2,ix1)=fcta8(ix,ix1,ix2)
                     fcta8(ix1,ix,ix2)=fcta8(ix,ix1,ix2)
                     fcta8(ix1,ix2,ix)=fcta8(ix,ix1,ix2)
                     fcta8(ix2,ix,ix1)=fcta8(ix,ix1,ix2)
                     fcta8(ix2,ix1,ix)=fcta8(ix,ix1,ix2)

                  end do
               end do
            end do
         end do !ikg

          
          !TEST: Check for NaN
          if (fcta4 /= fcta4) then
             write(*,*) 'ilm',ilm,'iat',iat,'fta1',fcta1
             write(*,*) 'ilm',ilm,'iat',iat,'fta4',fcta4
             stop
          end if
          do ix=1,3
             if (fcta2(ix) /= fcta2(ix)) then
                write(*,*) 'ilm',ilm,'iat',iat,'fta3',fcta3
                stop
             end if
             do ix1=1,3
                if (fcta6(ix,ix1) /= fcta6(ix,ix1)) then
                   write(*,*) 'ilm',ilm,'iat',iat,'fta5',fcta5
                   stop
                end if
                do ix2=1,3
                   if (fcta8(ix,ix1,ix2) /= fcta8(ix,ix1,ix2)) then
                      write(*,*) 'ilm',ilm,'iat',iat,'fcta7',fcta7
                      stop
                   end if
                end do
             end do
          end do
          
          !write(*,*) 'No NaN in ftas'
          !write (*,'(a15,6e12.2e2)') 'joper fcta',fcta5(3,:)
          !write (*,'(a10,6e16.4e2)') 'joper fcta3',fcta3
          !write (*,'(a10,2e16.4e2)') 'fcta5',fcta5
          !write (*,'(a10,2e16.4e2)') 'fcta7',fcta7
          
          
          ! Correction term expansion
          do ikg=1,npw1

             ! Convert ffnl to explict x,y,z indicies
             do ix=1,3
                ffnl1(ix)=ffnl(ikg,ix+1,ilm)
              end do
             do ix=1,6
                ffnl2(alpha(ix),beta(ix))=ffnl(ikg,ix+4,ilm)
                ffnl2(beta(ix),alpha(ix))=ffnl(ikg,ix+4,ilm)
             end do
             do ix=1,10
                ffnl3(alpha3(ix),beta3(ix),gamma3(ix))=ffnl(ikg,ix+10,ilm)
                ffnl3(alpha3(ix),gamma3(ix),beta3(ix))=ffnl(ikg,ix+10,ilm)
                ffnl3(beta3(ix),alpha3(ix),gamma3(ix))=ffnl(ikg,ix+10,ilm)
                ffnl3(beta3(ix),gamma3(ix),alpha3(ix))=ffnl(ikg,ix+10,ilm)
                ffnl3(gamma3(ix),alpha3(ix),beta3(ix))=ffnl(ikg,ix+10,ilm)
                ffnl3(gamma3(ix),beta3(ix),alpha3(ix))=ffnl(ikg,ix+10,ilm)
             end do

             ! zeroth order nl induced polarization
             vnlg0(:,ikg)=vnlg0(:,ikg)-4*( &
                  & ffnl(ikg,1,ilm)*fcta2(:)+ffnl(ikg,2:4,ilm)*fcta4)

             do ix=1,3 ! alpha
                do ix1=1,3 !beta
                   
                   !ICL path
                   if (dtset%useric==0) then
                      vnlg1(ix,ikg)=vnlg1(ix,ikg)-2.*qpc(ix1)*( &
                           &  ffnl1(ix)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl1(ix1)*fcta2(ix)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl(ikg,1,ilm)*fcta6(ix,ix1) &
                           & +ffnl2(ix,ix1)*fcta4)
                     
                      !PM path
                   else if (dtset%useric==1) then
                      vnlg1(ix,ikg)=vnlg1(ix,ikg)-2.*qpc(ix1)*( &
                           &  2.*ffnl1(ix)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl(ikg,1,ilm)*fcta6(ix,ix1) &
                           & +ffnl2(ix,ix1)*fcta4)

                   end if



                   do ix2=1,3 ! gamma
                      if (dtset%useric==0) then
                         ! ICL path
                         vnlg2(ix,ikg)=vnlg2(ix,ikg)-(2./3.)*qpc(ix1)*qpc(ix2)*( &
                              &  ffnl(ikg,1,ilm)*fcta8(ix,ix1,ix2) &
                              & +ffnl1(ix)*fcta6(ix1,ix2)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl1(ix1)*fcta6(ix2,ix)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl1(ix2)*fcta6(ix,ix1)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl2(ix,ix1)*fcta2(ix2)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl2(ix,ix2)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl2(ix1,ix2)*fcta2(ix)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl3(ix,ix1,ix2)*fcta4)            
                         ! PM path
                      else if (dtset%useric==1) then
                         vnlg2(ix,ikg)=vnlg2(ix,ikg)-(2./3.)*qpc(ix1)*qpc(ix2)*( &
                              &  ffnl(ikg,1,ilm)*fcta8(ix,ix1,ix2) &
                              & +3.*ffnl1(ix)*fcta6(ix1,ix2)*sign_dyad(iln(iat,ilm),ityp) &
                              & +3.*ffnl2(ix,ix2)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                              & +ffnl3(ix,ix1,ix2)*fcta4)
                      end if
                   end do !ix2
                end do !ix1
             end do ! ix

          end do !ikg
       end do !ilm
    end do !iat

 end if !calcnl
 
!TEST: JUST CHECK VNLG0
!write(*,*) 'JUST vnlg0 and vnlg1!!!!!'
!vnlg1=zero
!vnlg2=zero


 !Sum contributions. For now just in rfdir direction
 cwave0_npw1=zero
 do ikg=1,npw1
    !TEST:
    !if (calcnl<2) cwave0_npw1(:,ikg)=-gsgd(irfdir,:,ikg) 
    if (calcnl<2) then
       cwave0_npw1(:,1,ikg)=-4*gsgd(:,1,ikg)
       cwave0_npw1(:,2,ikg)=-4*gsgd(:,2,ikg)
    end if

    if (calcnl>0) then
       cwave0_npw1(:,1,ikg)= cwave0_npw1(:,1,ikg) &
            & +real(real(vnlg0(:,ikg)+vnlg1(:,ikg)+vnlg2(:,ikg)))
       cwave0_npw1(:,2,ikg)= cwave0_npw1(:,2,ikg) &
            & +real(aimag(vnlg0(:,ikg)+vnlg1(:,ikg)+vnlg2(:,ikg)))
    end if
 end do


!timing

!call timab(101,4,tottim)
!call timein(cpuout,wallout)
!write(*,*) 'time in joper:',wallout-wallin
!if (mpi_enreg%me_kpt==0) write(*,*) 'time in joper:',wallout-wallin

end subroutine joper
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkffnl
!! NAME
!! mkffnl
!!
!! FUNCTION
!! Make FFNL, nonlocal form factors, for each type of atom up to ntypat
!! and for each angular momentum.
!! When Legendre polynomials are used in the application of the
!!   nonlocal operator, FFNLs depend on (l,n) components; in this
!!   case, form factors are real and divided by |k+G|^l;
!! When spherical harmonics are used, FFNLs depend on (l,m,n)
!!   components; in this case, form factors are multiplied by Ylm(k+G).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, MT, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=described below
!!
!!
!! NOTES
!!
!! CEDrev: FULL CARTESIAN VERSION
!!
!! TODO
!!  Some parts can be rewritten with BLAS1 calls.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cmkffnl(calcpm,dtfil,dtset,ffnl,iat,kg,kpt,mpi_enreg,npw,Psps,sign_dyad)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cmkffnl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,calcpm,iat
! type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(in) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
!arrays
 integer,intent(in) :: kg(3,npw)
 integer,intent(out) :: sign_dyad(Psps%lnmax,Psps%ntypat)
 real(dp),intent(in) :: kpt(3) 
 complex(dp),intent(out) :: ffnl(npw,20,Psps%lnmax**2)
!type(kb_potential),intent(out) :: KBgrad_k

!Local variables-------------------------------
!scalars
 integer :: iffnl,ig,ig0,il,ilm,ilmn,iln,iln0,im,itypat,mu,mua,mub,muc,nlmn,imm
 real(dp),parameter :: renorm_factor=0.5d0/pi**2,tol_norm=tol10
 real(dp) :: ecut,ecutsm,effmass,fact,kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3,rmetab,yp1,factor
 logical :: testnl=.false.
 character(len=500) :: message
!arrays

!CEDrev: FOR TESTING, NEED TO CHANGE BACK WHEN USING INDPOL
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer,parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer,parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer,parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
! integer,allocatable :: sign_dyad(:,:)
 real(dp) :: gprimd(3,3),rprimd(3,3),ucvol


 real(dp) :: gmet(3,3),rmet(3,3),gam(3,3,npw),rred(3,dtset%natom),rcart(3)
 real(dp) :: tsec(2)
 real(dp) :: xdotg(npw),gcart(3,npw)
 real(dp),allocatable :: kpgc(:,:),kpgn(:,:),kpgnorm(:),kpgnorm_inv(:),wk_ffnl1(:)
 real(dp),allocatable :: wk_ffnl2(:),wk_ffnl3(:),wk_ffspl(:,:),wk_ffnl11(:),wk_ffnl12(:),wk_ffnl13(:),wk_ffspl1(:,:)
 complex(dp) :: sfac(npw),ffnl2(3,3),ylm2(3,3)
 complex(dp) :: ffnlpm(npw,20,Psps%lnmax**2,dtset%natom)
 complex(dp),allocatable :: ylm(:,:),ylm_gr(:,:,:)

!CEDrev:
logical :: filexist

! *************************************************************************

! DBG_ENTER("COLL")

!Keep track of time spent in mkffnl
 call timab(16,1,tsec)

!AB9rev: read this in

! Get metric, etc.
call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Get complex spherical harmonics:
ABI_MALLOC(ylm,(npw,Psps%mpsang*Psps%mpsang))
ABI_MALLOC(ylm_gr,(npw,27,Psps%mpsang*Psps%mpsang))
call cinitylmgi(gprimd,kg,kpt,1,mpi_enreg,Psps%mpsang,npw,&
&  dtset%nband(1),1,(/npw/),dtset%nsppol,8,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylm_gr)


! Get dyadic sign
do itypat=1,Psps%ntypat
   nlmn=count(Psps%indlmn(3,:,itypat)>0)
   iln0=0
   do ilmn=1,nlmn
      iln=Psps%indlmn(5,ilmn,itypat)
      if (iln>iln0) then
         iln0=iln
         sign_dyad(iln,itypat)=nint(DSIGN(one,Psps%ekb(ilmn,itypat)))

         ! TEST:
         !write(*,*) "SIGN",iln,itypat,Psps%ekb(ilmn,itypat),sign_dyad(iln,itypat)

      end if
   end do
end do

!Get (k+G) and |k+G|:
 ABI_MALLOC(kpgnorm,(npw))
 ABI_MALLOC(kpgnorm_inv,(npw))
 ig0=-1 ! index of |k+g|=0 vector
 ABI_MALLOC(kpgc,(npw,3))
 ABI_MALLOC(kpgn,(npw,3))
 do ig=1,npw
    kpg1=kpt(1)+dble(kg(1,ig))
    kpg2=kpt(2)+dble(kg(2,ig))
    kpg3=kpt(3)+dble(kg(3,ig))

    ! leave out 2pi here for splfit, then add in to gam and kpgnorm_inv in the second order term
    kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
    kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
    kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
    kpgc(ig,1)=kpgc1
    kpgc(ig,2)=kpgc2
    kpgc(ig,3)=kpgc3
    kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
    if (kpgnorm(ig)<=tol_norm) ig0=ig
    kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)

    !All cartesian, all the time
    kpgn(ig,1)=kpgc1*kpgnorm_inv(ig)
    kpgn(ig,2)=kpgc2*kpgnorm_inv(ig)
    kpgn(ig,3)=kpgc3*kpgnorm_inv(ig)

    ! For structure factor
    gcart(:,ig)=real(kg(1,ig))*gprimd(:,1)+real(kg(2,ig))*gprimd(:,2)+real(kg(3,ig))*gprimd(:,3)

    ! For second derivative:
    do mua=1,3
       do mub=1,3
          !If EVERYTHING is cartesian, rmet is delta function
          if (mua==mub) then
             gam(mua,mub,ig)=(one-kpgn(ig,mua)*kpgn(ig,mub))*kpgnorm_inv(ig)/two_pi  
          else             
             gam(mua,mub,ig)=(-kpgn(ig,mua)*kpgn(ig,mub))*kpgnorm_inv(ig)/two_pi
          end if
       end do
    end do
  
 end do

! Arrays for radial part of ff and RADIAL derivatives 
 ABI_MALLOC(wk_ffnl1,(npw))
 ABI_MALLOC(wk_ffnl2,(npw))
 ABI_MALLOC(wk_ffnl3,(npw))
 ABI_MALLOC(wk_ffspl,(Psps%mqgrid_ff,2))
 ABI_MALLOC(wk_ffnl11,(npw))
 ABI_MALLOC(wk_ffnl12,(npw))
 ABI_MALLOC(wk_ffnl13,(npw))
 ABI_MALLOC(wk_ffspl1,(Psps%mqgrid_ff,2))

! NO LOOP OVER ATOMS!
! do ia=1,dtset%natom
    itypat=dtset%typat(iat)

    ! Structure factor
    rcart(:)=matmul(rprimd(:,:),dtset%xred_orig(:,iat,1))
    xdotg(:)=two_pi*(gcart(1,:)*rcart(1)+gcart(2,:)*rcart(2)+gcart(3,:)*rcart(3))
    sfac(:)=cmplx(cos(-1.*xdotg(:)),sin(-1.*xdotg(:)))

    !  Loop over (l,m,n) values
    iln0=0;nlmn=count(Psps%indlmn(3,:,itypat)>0)
    ! zero index for ffnl
    iffnl=0
    do ilmn=1,nlmn
       il=Psps%indlmn(1,ilmn,itypat)+1
       !     im=Psps%indlmn(2,ilmn,itypat)
       ilm =Psps%indlmn(4,ilmn,itypat)
       iln =Psps%indlmn(5,ilmn,itypat)

       factor=four_pi*sqrt(abs(Psps%ekb(iln,itypat))/ucvol)

       ! Store form factors (from ffspl)
       do ig=1,Psps%mqgrid_ff
          wk_ffspl(ig,1)=PsPs%ffspl(ig,1,iln,itypat)
          wk_ffspl(ig,2)=Psps%ffspl(ig,2,iln,itypat)
          ! For third derivative, 1st and 3rd of ff
          wk_ffspl1(ig,1)=PsPs%ffspl1(ig,1,iln,itypat)
          wk_ffspl1(ig,2)=Psps%ffspl1(ig,2,iln,itypat)
       end do
       call splfit(Psps%qgrid_ff,wk_ffnl2,wk_ffspl,1,kpgnorm,wk_ffnl1,Psps%mqgrid_ff,npw)
       call splfit(Psps%qgrid_ff,wk_ffnl3,wk_ffspl,2,kpgnorm,wk_ffnl1,Psps%mqgrid_ff,npw)

       call splfit(Psps%qgrid_ff,wk_ffnl12,wk_ffspl1,1,kpgnorm,wk_ffnl11,Psps%mqgrid_ff,npw)
       call splfit(Psps%qgrid_ff,wk_ffnl13,wk_ffspl1,2,kpgnorm,wk_ffnl11,Psps%mqgrid_ff,npw)
      

       ! To account for the 1/2pi in d/dK
       wk_ffnl2=wk_ffnl2/two_pi
       wk_ffnl3=wk_ffnl3/two_pi**2
       wk_ffnl13=wk_ffnl13/two_pi**3

       ! imm is index for ylm
       imm=(il-1)*(il-1)
       ! Loop over m
       do im=1,2*(il-1)+1

          ! Indicies
          iffnl=iffnl+1
          imm=imm+1
          
          !TEST: Check the indicies
          !write(*,*) 'MKFFNL il',il,'ilm', ilm,'iln',iln,'iffnl',iffnl, 'im',im, 'imm',imm

          do ig=1,npw
             ! Direction of derivatives:
             do mu=1,10
                
                ! Potential (only do once):
                if (mu==1) ffnl(ig,1,iffnl)=factor*sfac(ig)*ylm(ig,imm)*wk_ffnl1(ig)*real(sign_dyad(iln,itypat))
                
                ! TEST:
                !write(*,'(a5,3i5)') "MU1",iln,itypat,sign_dyad(iln,itypat)  
                !factor,sfac(ig),ylm(ig,imm),wk_ffnl1(ig),sign_dyad(iln,itypat)

                ! First derivative (three components):
                if (mu <= 3) ffnl(ig,mu+1,iffnl)=factor*sfac(ig)*(ylm(ig,imm)*wk_ffnl2(ig)*kpgn(ig,mu)+ylm_gr(ig,mu,imm)*wk_ffnl1(ig))
                
                ! Second derivative (6 components)
                if (mu <= 6) then
                   mua=alpha(mu); mub=beta(mu)
                   
                   ffnl(ig,mu+4,iffnl)=factor*sfac(ig)*(&
                        &   ylm_gr(ig,3+mu,imm)*wk_ffnl1(ig) &
                        & + gam(mua,mub,ig)*ylm(ig,imm)*wk_ffnl2(ig) &
                        & + ylm(ig,imm)*kpgn(ig,mua)*kpgn(ig,mub)*wk_ffnl3(ig) &
                        & + ylm_gr(ig,mua,imm)*kpgn(ig,mub)*wk_ffnl2(ig) &
                        & + ylm_gr(ig,mub,imm)*kpgn(ig,mua)*wk_ffnl2(ig))

                end if
                ylm2(1,1)=ylm_gr(ig,4,imm)
                ylm2(1,2)=ylm_gr(ig,5,imm)
                ylm2(2,1)=ylm_gr(ig,5,imm)
                ylm2(1,3)=ylm_gr(ig,6,imm)
                ylm2(3,1)=ylm_gr(ig,6,imm)
                ylm2(2,2)=ylm_gr(ig,7,imm)
                ylm2(2,3)=ylm_gr(ig,8,imm)
                ylm2(3,2)=ylm_gr(ig,8,imm)
                ylm2(3,3)=ylm_gr(ig,9,imm)


                mua=alpha3(mu); mub=beta3(mu); muc=gamma3(mu)

                ffnl(ig,mu+10,iffnl) = factor*sfac(ig)*( &
                     &  kpgn(ig,mua)*kpgn(ig,mub)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl13(ig) & !1
                     & +ylm_gr(ig,mu+9,imm)*wk_ffnl1(ig) & !4
                     !
                     & +kpgn(ig,mua)*kpgn(ig,mub)*ylm_gr(ig,muc,imm)*wk_ffnl3(ig) & !1
                     & +kpgn(ig,mua)*kpgn(ig,muc)*ylm_gr(ig,mub,imm)*wk_ffnl3(ig) & !2
                     & +kpgn(ig,mub)*kpgn(ig,muc)*ylm_gr(ig,mua,imm)*wk_ffnl3(ig) & !3
                     !
                     & +kpgn(ig,mua)*ylm2(mub,muc)*wk_ffnl2(ig) & !2
                     & +kpgn(ig,muc)*ylm2(mua,mub)*wk_ffnl2(ig) & !4
                     & +kpgn(ig,mub)*ylm2(mua,muc)*wk_ffnl2(ig) & !3
                     ! Remember, this is all in cartesian coordinates
                     & +gam(mua,muc,ig)*kpgn(ig,mub)*ylm(ig,imm)*wk_ffnl3(ig) & !1
                     & +gam(mub,muc,ig)*kpgn(ig,mua)*ylm(ig,imm)*wk_ffnl3(ig) & !1
                     & +gam(mua,mub,ig)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl3(ig) & !5
                     !
                     & +gam(mua,muc,ig)*ylm_gr(ig,mub,imm)*wk_ffnl2(ig) & !2
                     & +gam(mub,muc,ig)*ylm_gr(ig,mua,imm)*wk_ffnl2(ig) & !3
                     & +gam(mua,mub,ig)*ylm_gr(ig,muc,imm)*wk_ffnl2(ig) & !5
                     ! WHY 1/2pi????
                     & -gam(mua,muc,ig)*kpgn(ig,mub)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi & !5.2
                     & -gam(mub,muc,ig)*kpgn(ig,mua)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi & !5.3
                     & -gam(mua,mub,ig)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi) !5


             end do !mu
             
          end do !ig          
       end do !im
    end do !ilmn
! end do !itypat

 call timab(16,2,tsec)

! DBG_EXIT("COLL")

end subroutine cmkffnl
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/initylmg
!! NAME
!! initylmg
!!
!! FUNCTION
!! Calculate the real spherical harmonics Ylm (and gradients)
!! over a set of (reciprocal space) (k+G) vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive
!!              translations (b^-1)
!!  kg(3,mpw)=integer coordinates of G vectors in basis sphere
!!  kptns(3,nkpt)=k points in terms of reciprocal translations
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere
!!         (large number)
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt  =number of k points
!!  npwarr(nkpt)=array holding npw for each k point
!!  nsppol=1 for unpolarized, 2 for polarized
!!  optder= 0=compute Ylm(K)
!!          1=compute Ylm(K) and dYlm/dKi
!!          2=compute Ylm(K), dYlm/dKi and d2Ylm/dKidKj
!!         -1=compute only dYlm/dKi
!!  rprimd(3,3)=dimensional primitive translations in real space
!!              (bohr)
!!  unkg=unit number for (k+G) sphere data
!!  unylm=unit number for storage of Ylm on disk
!!
!! OUTPUT
!!  if (optder>=0)
!!    ylm(mpw*mkmem,mpsang*mpsang) = real spherical harmonics
!!    for each G and k point
!!  if (optder>=1 or optder==-1)
!!    ylm_gr(mpw*mkmem,1:3,mpsang*mpsang)= gradients of real
!!    spherical harmonics wrt (G+k) in reduced coordinates
!!  if (optder>=2)
!!    ylm_gr(mpw*mkmem,4:9,mpsang*mpsang)= second gradients of
!!    real spherical harmonics wrt (G+k) in reduced coordinates
!!
!! NOTES
!! Remember the expression of complex spherical harmonics:
!! $Y_{lm}(%theta ,%phi)=sqrt{{(2l+1) over (4 %pi)}
!! {fact(l-m) over fact(l+m)} } P_l^m(cos(%theta))
!! func e^{i m %phi}$
!! Remember the expression of real spherical harmonics as
!!   linear combination of complex spherical harmonics:
!! $Yr_{lm}(%theta ,%phi)=(Re{Y_{l-m}}+(-1)^m Re{Y_{lm}})/sqrt{2}
!! $Yr_{l-m}(%theta ,%phi)=(Im{Y_{l-m}}-(-1)^m Im{Y_{lm}})/sqrt{2}
!!
!! CEDrev: This is altered version of initylmg to give COMPLEX Ylm 
!!        and CARTESIAN derivative up to 3rd
!!
!! PARENTS
!!      debug_tools,gstate,ks_ddiago,loop3dte,loper3,m_paw_pwij,m_shirley,m_wfs
!!      mover,nstpaw3,partial_dos_fractions,respfn,scfcv,wffile
!!
!! CHILDREN
!!      plm_coeff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cinitylmgi(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,&
&  nband,nkpt,npwarr,nsppol,optder,rprimd,unkg,unylm,ylm,ylm_gr)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cinitylmgi3'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mkmem,mpsang,mpw,nkpt,nsppol,optder
 integer,intent(in) :: unkg,unylm
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(in) :: gprimd(3,3),kptns(3,nkpt),rprimd(3,3)
 complex(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
 complex(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)
!Local variables ------------------------------
!scalars
 integer :: dimgr,ia,ib,ii,ikg,ikpt,ilang,ipw
 integer :: jj,kk,l0,ll
 integer :: me_distrb,mm,npw_k
 real(dp),parameter :: tol=1.d-10
 real(dp) :: cphi,ctheta,fact,onem,rr,sphi,stheta,work1,work2
 real(dp) :: xx,ylmcst,ylmcst2
 real(dp) :: yy,zz
 !character(len=500) :: message
!arrays

!CEDrev: 
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/)
! integer,parameter :: beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: dphi(3),dtheta(3),iphase(mpsang-1),kpg(3)
 real(dp) :: rphase(mpsang-1)
 real(dp) :: blm(5,mpsang*mpsang)
 real(dp) :: ylmgr2_cart(3,3,2),ylmgr2_tmp(3,3)
 real(dp) :: ylmgr_cart(3,2)
 real(dp) :: ylmgr_red(10,2)

 !CEDrev: for thrid derivative
 integer,parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer,parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer,parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
 integer :: pp,qq,ss,aa,igam
 real(dp) :: plm_d3t
 real(dp) :: plm_d2t(mpsang*mpsang)
 real(dp) :: dydth(2),d2ydth2(2),d3ydth3(2),dydph(2),d2ydph2(2)
 real(dp) :: d3ydph3(2),d2ydthdph(2),d3ydth2dph(2),d3ydthdph2(2)
 real(dp) :: ylmgr3_cart(3,3,3,2),ylmgr3_tmp(3,3,3)

 !TEST
 real(dp) :: ylmtstconst,ylmgr2_cart_tst(3,3,2)

!*****************************************************************

!Begin executable
 me_distrb=mpi_enreg%me_kpt
!Initialisation of spherical harmonics (and gradients)
 if (optder>=0) ylm(:,:)  =cmplx(zero,zero)
 if (optder/=0) ylm_gr(:,:,:)=cmplx(zero,zero)

!CEDrev: make room for third derivatives:
if (optder<3) then
   dimgr=3+6*(optder/2)
else if (optder==8) then
   dimgr=27
end if
! dimgr=3+6*(optder/2)

!Allocate some memory
! if (optder/=0) then
!   ABI_MALLOC(ylmgr_cart,(3,2))
! end if
! if (optder/=0.and.optder/=2) then
!   ABI_MALLOC(ylmgr_red,(3,2))
! end if
! if (optder==2) then
!   ABI_MALLOC(ylmgr2_cart,(3,3,2))
!   ABI_MALLOC(ylmgr2_tmp,(3,3))
!   ABI_MALLOC(ylmgr_red,(6,2))
!   ABI_MALLOC(blm,(5,mpsang*mpsang))
! end if

!CEDrev: For third derivative
! if (optder==8) then
!   ABI_MALLOC(ylmgr2_cart,(3,3,2))
!   ABI_MALLOC(ylmgr2_tmp,(3,3))
!   ABI_MALLOC(ylmgr_red,(10,2))
!   ABI_MALLOC(blm,(5,mpsang*mpsang))
! end if

!Loop over k-points:
 ikg=0
 do ikpt=1,nkpt

   ! Don't need this since its only for one kpoint
   !if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband(ikpt),-1,me_distrb)) cycle 

!  Get k+G-vectors, for this k-point:
   npw_k=npwarr(ikpt)
   ABI_MALLOC(kg_k,(3,npw_k))
   kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

!  Special case for l=0
   if (optder>=0) ylm(1+ikg:npw_k+ikg,1)=cmplx(1._dp/sqrt(four_pi),0)
   if (optder/=0) ylm_gr(1+ikg:npw_k+ikg,1:dimgr,1)=zero

   if (mpsang>1) then
!    Loop over all k+G
      do ipw=1,npw_k
        
!      Load k+G
       kpg(1)=kptns(1,ikpt)+real(kg_k(1,ipw),dp)
       kpg(2)=kptns(2,ikpt)+real(kg_k(2,ipw),dp)
       kpg(3)=kptns(3,ikpt)+real(kg_k(3,ipw),dp)

!      Calculate mod of k+G
       xx=two_pi*(gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3))
       yy=two_pi*(gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3))
       zz=two_pi*(gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3))
       rr=sqrt(xx**2+yy**2+zz**2)
       
       !TEST
       !if (ipw==57)then
       !   write(*,*) xx,yy,zz
       !end if

!      Continue only for k+G<>0
       if (rr>tol) then

!        Determine theta and phi
         cphi=one
         sphi=zero
         ctheta=zz/rr
         stheta=sqrt(abs((one-ctheta)*(one+ctheta)))
         if (stheta>tol) then
           cphi=xx/(rr*stheta)
           sphi=yy/(rr*stheta)
         end if
         do mm=1,mpsang-1
           rphase(mm)=dreal(dcmplx(cphi,sphi)**mm)
           iphase(mm)=aimag(dcmplx(cphi,sphi)**mm)
         end do

!        Determine gradients of theta and phi
         if (optder/=0) then
           dtheta(1)=ctheta*cphi
           dtheta(2)=ctheta*sphi
           dtheta(3)=-stheta
           dphi(1)=-sphi
           dphi(2)=cphi
           dphi(3)=zero
         end if

!        COMPUTE Ylm(K) 
!        ============================================
         if (optder>=0) then
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)
!            Special case m=0
             ylm(ikg+ipw,l0)=cmplx(ylmcst*ass_leg_pol(ll,0,ctheta),0)
!            Compute for m>0
             onem=one
             do mm=1,ll
               onem=-onem
               work1=ylmcst*sqrt(fact)*ass_leg_pol(ll,mm,ctheta)!*sqrt(2._dp)
               ylm(ikg+ipw,l0+mm)=cmplx(work1*rphase(mm),work1*iphase(mm))
               !TEST
               ylm(ikg+ipw,l0-mm)=-conjg(ylm(ikg+ipw,l0+mm))
!               ylm(ikg+ipw,l0-mm)=onem*conjg(cmplx(work1*rphase(mm),work1*iphase(mm)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
             end do ! End loop over m
           end do  ! End loop over l
         end if

!        COMPUTE dYlm/dKi
!        ============================================
         if (optder/=0) then
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/rr
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             work1=ylmcst*plm_dtheta(ll,0,ctheta)
             ylmgr_cart(1:3,1)=work1*dtheta(1:3)
!            2-Transfer gradients into reduced coordinates
!             do ii=1,3
!               ylmgr_red(ii,1)=(rprimd(1,ii)*ylmgr_cart(1,1)+&
!&               rprimd(2,ii)*ylmgr_cart(2,1)+&
!&               rprimd(3,ii)*ylmgr_cart(3,1))
!             end do
!            3-Store gradients
             ylm_gr(ikg+ipw,1:3,l0) =cmplx(ylmgr_cart(1:3,1),0)
!             ylm_gr(ikg+ipw,1:3,l0) =cmplx(ylmgr_red(1:3,1),0)
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll
               onem=-onem
!              1-compute gradients in cartesian coordinates
               work1=ylmcst*sqrt(fact)*plm_dtheta(ll,mm,ctheta)
               work2=ylmcst*sqrt(fact)*plm_dphi  (ll,mm,ctheta)
               ylmgr_cart(1:3,1)=rphase(mm)*work1*dtheta(1:3)-iphase(mm)*work2*dphi(1:3)
               ylmgr_cart(1:3,2)=iphase(mm)*work1*dtheta(1:3)+rphase(mm)*work2*dphi(1:3)
               
               !TEST
               !if (ll==1.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=1,mm=1,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !else if (ll==2.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=1,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !else if (ll==2.and.mm==2.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=2,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !end if

!              2-Transfer gradients into reduced coordinates
!               do kk=1,2
!                 do ii=1,3
!                   ylmgr_red(ii,kk)=(rprimd(1,ii)*ylmgr_cart(1,kk)+&
!&                   rprimd(2,ii)*ylmgr_cart(2,kk)+&
!&                   rprimd(3,ii)*ylmgr_cart(3,kk))
!                 end do
!               end do
!              3-Store gradients
               ylm_gr(ikg+ipw,1:3,l0+mm) =cmplx(ylmgr_cart(1:3,1),ylmgr_cart(1:3,2))
               ylm_gr(ikg+ipw,1:3,l0-mm) =-conjg(ylm_gr(ikg+ipw,1:3,l0+mm))
!               ylm_gr(ikg+ipw,1:3,l0-mm) =onem*conjg(cmplx(ylmgr_cart(1:3,1),ylmgr_cart(1:3,2)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
            end do ! End loop over m
          end do  ! End loop over l
        end if

!        COMPUTE d2Ylm/dKidKj
!        ============================================
         ! CEDrev:
        if (optder>=2) then 
!         if (optder==2) then
           call plm_coeff(blm,mpsang,ctheta)
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/(rr**2)
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             ylmgr2_cart(1,1,1)=ylmcst*(-blm(3,l0)*sphi*sphi+blm(4,l0)*cphi*cphi)
             ylmgr2_cart(2,2,1)=ylmcst*(-blm(3,l0)*cphi*cphi+blm(4,l0)*sphi*sphi)
             ylmgr2_cart(3,3,1)=ylmcst*blm(1,l0)
             ylmgr2_cart(3,1,1)=ylmcst*blm(2,l0)*cphi
             ylmgr2_cart(3,2,1)=ylmcst*blm(2,l0)*sphi
             ylmgr2_cart(2,1,1)=ylmcst*(blm(3,l0)+blm(4,l0))*sphi*cphi
             ylmgr2_cart(1,3,1)=ylmgr2_cart(3,1,1)
             ylmgr2_cart(1,2,1)=ylmgr2_cart(2,1,1)
             ylmgr2_cart(2,3,1)=ylmgr2_cart(3,2,1)

             !TEST: Check against my derivatives in cartesian coordinates
             !if (ll==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
             !   ylmtstconst=-sqrt(3./pi)/(2*rr**5)
             !   ylmgr2_cart_tst(1,1,1)=ylmtstconst*zz*(-2.*xx**2+yy**2+zz**2)
             !   ylmgr2_cart_tst(2,2,1)=ylmtstconst*zz*(xx**2-2.*yy**2+zz**2)
             !   ylmgr2_cart_tst(3,3,1)=ylmtstconst*3.*zz*(xx**2+yy**2)
             !   ylmgr2_cart_tst(3,1,1)=ylmtstconst*xx*(xx**2+yy**2-2.*zz**2)
             !   ylmgr2_cart_tst(3,2,1)=ylmtstconst*yy*(xx**2+yy**2-2.*zz**2)
             !   ylmgr2_cart_tst(2,1,1)=-ylmtstconst*(xx*yy*zz)
             !   ylmgr2_cart_tst(1,3,1)=ylmgr2_cart_tst(3,1,1)
             !   ylmgr2_cart_tst(1,2,1)=ylmgr2_cart_tst(2,1,1)
             !   ylmgr2_cart_tst(2,3,1)= ylmgr2_cart_tst(3,2,1)
                
             !   write(*,'(3i5)') kg_k(:,ipw)
                !write(*,'(a15,2e16.4e2)') 'ylm2 11',ylmgr2_cart(1,1,1),ylmgr2_cart_tst(1,1,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 22', ylmgr2_cart(2,2,1),ylmgr2_cart_tst(2,2,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 33', ylmgr2_cart(3,3,1),ylmgr2_cart_tst(3,3,1)
             !   write(*,'(a15,2e16.4e2)') 'ylm2 31', ylmgr2_cart(3,1,1),ylmgr2_cart_tst(3,1,1)
             !   write(*,*) 'l0',l0,'cphi',cphi,'blm(2,l0)',blm(2,l0)
             !   write(*,*) 'rr',rr,'xx',xx
                !write(*,'(a15,2e16.4e2)') 'ylm2 32', ylmgr2_cart(3,2,1),ylmgr2_cart_tst(3,2,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 21', ylmgr2_cart(2,1,1),ylmgr2_cart_tst(2,1,1)
                !stop
             !end if

!            2-Transfer gradients into reduced coordinates
!             do jj=1,3
!               do ii=1,3
!                 ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,1)+&
!&                 rprimd(2,jj)*ylmgr2_cart(2,ii,1)+&
!&                 rprimd(3,jj)*ylmgr2_cart(3,ii,1))
!               end do
          !end do
          do ii=1,6
               ia=alpha(ii);ib=beta(ii)
!               ylm_gr(ikg+ipw,4:9,l0) =cmplx(ylmgr2_cart(ia,ib,1),0)
               ylm_gr(ikg+ipw,3+ii,l0) =cmplx(ylmgr2_cart(ia,ib,1),0) 

!               ylmgr_red(ii,1)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
!&               rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
!&               rprimd(3,ia)*ylmgr2_tmp(3,ib))
            end do
             !CEDrev: DO I NEED A SQRT 2 HERE????
!             ylm_gr(ikg+ipw,4:9,l0) =cmplx(ylmgr_cart(1:6,1),0)
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll                
                onem=-onem
               ylmcst2=ylmcst*sqrt(fact)!*sqrt(two)

               ylmgr2_cart(1,1,1)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*rphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(1,1,2)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*iphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(2,2,1)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*rphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(2,2,2)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*iphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(3,3,1)=ylmcst2*blm(1,l0+mm)*rphase(mm)
               ylmgr2_cart(3,3,2)=ylmcst2*blm(1,l0+mm)*iphase(mm)
               ylmgr2_cart(3,1,1)=ylmcst2*(blm(2,l0+mm)*cphi*rphase(mm)-&
&               mm*iphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,1,2)=ylmcst2*(blm(2,l0+mm)*cphi*iphase(mm)+&
&               mm*rphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,2,1)=ylmcst2*(blm(2,l0+mm)*sphi*rphase(mm)+&
&               mm*iphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
              ylmgr2_cart(3,2,2)=ylmcst2*(blm(2,l0+mm)*sphi*iphase(mm)-&
&               mm*rphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(2,1,1)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*rphase(mm)-&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*iphase(mm))
               ylmgr2_cart(2,1,2)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*iphase(mm)+&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*rphase(mm))
               ylmgr2_cart(1,3,:)=ylmgr2_cart(3,1,:)
               ylmgr2_cart(1,2,:)=ylmgr2_cart(2,1,:)
               ylmgr2_cart(2,3,:)=ylmgr2_cart(3,2,:)
!              2-Transfer gradients into reduced coordinates
!               do kk=1,2
!                 do jj=1,3
!                   do ii=1,3
!                     ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,kk)+&
!&                     rprimd(2,jj)*ylmgr2_cart(2,ii,kk)+&
!&                     rprimd(3,jj)*ylmgr2_cart(3,ii,kk))
!                   end do
               !                 end do

               !TEST
               !if (ll==1.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   ylmtstconst=sqrt(3./two_pi)/(2*rr**5)
               !   ylmgr2_cart_tst(3,1,1)=ylmtstconst*zz*(-2.*xx**2+yy**2+zz**2)
               !   ylmgr2_cart_tst(3,1,2)=ylmtstconst*zz*(-3.*xx*yy)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 31', onem*ylmgr2_cart(3,1,1),ylmgr2_cart_tst(3,1,1),onem*ylmgr2_cart(3,1,2),ylmgr2_cart_tst(3,1,2)
               !   write(*,'(a15,4e16.4e2)') 'll=1,mm=1 ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 12', ylmgr2_cart(2,1,1),ylmgr2_cart(2,1,2)
                  !write(*,*) 'l0+mm',l0+mm,'cphi',cphi,'blm(2,l0+mm)',blm(2,l0+mm)
                  !write(*,*) 'rr',rr,'xx',xx
                  !stop
               !end if
               !if (ll==2.and.mm==2.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
                  !write(*,'(a15,4e16.4e2)') 'ylm2 31', ylmgr2_cart(3,1,1),ylmgr2_cart(3,1,2)
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=2 ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 12', ylmgr2_cart(2,1,1),ylmgr2_cart(2,1,2)
                  !stop
               !end if
               !if (ll==2.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=1,ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
               !end if


               do ii=1,6
                  ia=alpha(ii);ib=beta(ii)
                  
                  !TEST
                  ylm_gr(ikg+ipw,3+ii,l0+mm) =onem*cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2))
                  ylm_gr(ikg+ipw,3+ii,l0-mm) =-conjg(ylm_gr(ikg+ipw,3+ii,l0+mm)) 

                  !ylm_gr(ikg+ipw,3+ii,l0+mm) =cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2))
                  !ylm_gr(ikg+ipw,3+ii,l0-mm) =onem*conjg(cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2)))

                  ! ylmgr_red(ii,kk)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
                  !& rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
                  !& rprimd(3,ia)*ylmgr2_tmp(3,ib))
               end do
               !              end do
               !ylm_gr(ikg+ipw,4:9,l0+mm) =cmplx(ylmgr_cart(1:6,1),ylmgr_cart(1:6,2))
               !ylm_gr(ikg+ipw,4:9,l0-mm) =onem*conjg(cmplx(ylmgr_cart(1:6,1),ylmgr_cart(1:6,2)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
            end do ! End loop over m
         end do  ! End loop over l
      end if

!CEDrev: Since I'm not converting to reduced coordinates, I can clean this up significantly
!        COMPUTE d3Ylm/dKidKjKn
!        ============================================
        if (optder==8) then
           ! Loop over angular momentum l
           do ilang=2,mpsang
              ll=ilang-1
              l0=ll**2+ll+1
              ! === Special case m=0 ===
              
              ylmgr3_cart=zero

              if (ll==1) then
                 ylmcst=sqrt(3./(pi*4.))/(rr**7)
                 !xxx
                 ylmgr3_cart(1,1,1,1)=ylmcst*(-6.*zz*xx**3+9.*xx*zz*(yy**2+zz**2))
                 !xxy
                 ylmgr3_cart(1,1,2,1)=ylmcst*(3.*yy*zz*(-4.*xx**2+yy**2+zz**2))
                 ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                 ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                 !xxz
                 ylmgr3_cart(1,1,3,1)=ylmcst*(2.*xx**4-yy**4+(yy**2)*(zz**2)+2.*zz**4+(xx**2)*(yy**2-11.*zz**2))
                 ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                 ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                 !xyy
                 ylmgr3_cart(1,2,2,1)=ylmcst*(3.*xx*zz*(xx**2-4.*yy**2+zz**2))
                 ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                 ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                 !xyz
                 ylmgr3_cart(1,2,3,1)=ylmcst*(3.*xx*yy*(xx**2+yy**2-4.*zz**2))
                 ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                 !yyy
                 ylmgr3_cart(2,2,2,1)=ylmcst*(3.*yy*zz*(3.*xx**2-2.*yy**2+3.*zz**2))
                 !yyz
                 ylmgr3_cart(2,2,3,1)=ylmcst*(-xx**4+2.*yy**4-11.*(yy**2)*(zz**2)+2.*zz**4+(xx**2)*(yy**2+zz**2))
                 ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                 ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                 !yzz
                 ylmgr3_cart(2,3,3,1)=ylmcst*(3.*zz*yy*(3.*xx**2+3.*yy**2-2.*zz**2))
                 ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                 ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                 !zzx
                 ylmgr3_cart(3,3,1,1)=ylmcst*(3.*xx*zz*(3.*xx**2+3.*yy**2-2.*zz**2))
                 ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                 ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                 !zzz
                 ylmgr3_cart(3,3,3,1)=-ylmcst*3.*(xx**2+yy**2)*(xx**2+yy**2-4.*zz**2)

              else if (ll==2) then
                 ylmcst=sqrt(5./(pi*16.))/(rr**8)
                 !xxx
                 ylmgr3_cart(1,1,1,1)=ylmcst*72.*xx*(zz**2)*(-xx**2+yy**2+zz**2)
                 !xxy
                 ylmgr3_cart(1,1,2,1)=ylmcst*24.*(-5.*xx**2+yy**2+zz**2)
                 ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                 ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                 !xxz
                 ylmgr3_cart(1,1,3,1)=ylmcst*12.*zz*(3.*xx**4-yy**4+zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                 ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                 !xyy
                 ylmgr3_cart(1,2,2,1)=ylmcst*24.*xx*(zz**2)*(xx**2-5.*yy**2+zz**2)
                 ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                 ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                 !xyz
                 ylmgr3_cart(1,2,3,1)=ylmcst*48.*xx*yy*zz*(xx**2+yy**2-2.*zz**2)
                 ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                 !yyy
                 ylmgr3_cart(2,2,2,1)=ylmcst*72.*yy*(zz**2)*(xx**2-yy**2+zz**2)
                 !yyz
                 ylmgr3_cart(2,2,3,1)=ylmcst*12.*zz*(-xx**4+2.*(xx**2)*(yy**2)+3.*(yy**4)-8.*(yy**2)*(zz**2)+zz**4)
                 ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                 ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                 !yzz
                 ylmgr3_cart(2,3,3,1)=-ylmcst*12.*yy*(xx**4+yy**4-8.*(yy**2)*(zz**2)+3.*zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                 ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                 !zzx
                 ylmgr3_cart(3,3,1,1)=-ylmcst*12.*xx*(xx**4+yy**4-8.*(yy**2)*(zz**2)+3.*zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                 ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                 !zzz
                 ylmgr3_cart(3,3,3,1)=-ylmcst*72.*zz*(xx**2+yy**2)*(xx**2+yy**2-zz**2)

              else if (ll==3) then
                 write(*,*) 'Third deriv of ylm not implemented for f electrons'
                 stop
              end if

              ! 2-Transfer gradients into reduced coordinates
              !ylmgr3_tmp=zero
              !do ii=1,3
              !   do jj=1,3
              !      do kk=1,3
              !         do pp=1,3
              !            do qq=1,3
              !               do ss=1,3
              !                  ylmgr3_tmp(ii,jj,kk)=ylmgr3_tmp(ii,jj,kk) &
              !                       & +rprimd(pp,ii)*rprimd(qq,jj)*rprimd(ss,kk)*ylmgr3_cart(pp,qq,ss,1)
              !               end do
              !            end do
              !         end do
              !      end do
              !   end do
              !end do
              do ii=1,10
                 ia=alpha3(ii);ib=beta3(ii);igam=gamma3(ii)
                 ylm_gr(ikg+ipw,9+ii,l0)=cmplx(ylmgr3_cart(ia,ib,igam,1),0)
                 !ylmgr_cart(ii,1)=ylmgr3_tmp(ia,ib,igam) 
                 !ylmgr_red(ii,1)=ylmgr3_tmp(ia,ib,igam)
              end do
              !ylm_gr(ikg+ipw,10:19,l0)=cmplx(ylmgr_red(1:10,1),0)
              ! === Compute for m>0 ===
              onem=one
              do mm=1,ll
                 onem=-onem

                 ! === elements ===
                 ylmgr3_cart=zero
                 
                 if (ll==1 .and. mm==1) then

                    ylmcst=sqrt(3./(pi*8.))/(rr**7)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*3.*(-4.*(xx**2)*(yy**2+zz**2)+(yy**2+zz**2)**2)
                    ylmgr3_cart(1,1,1,2)=ylmcst*3.*(2.*(xx**3)*yy-3.*xx*yy*(yy**2+zz**2))
                    !xxy
                    ylmgr3_cart(1,1,2,1)=-ylmcst*(-6.*yy*xx**3+9.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,1,2,2)=-ylmcst*(2.*xx**4+2.*yy**4+(yy**2)*(zz**2)-zz**4+(xx**2)*(-11.*yy**2+zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*3.*zz*(-2.*xx**3+3.*xx*(yy**2+zz**2))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*3.*zz*(-4.*yy*xx**2+yy*(yy**2+zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*(-2.*xx**4-2.*yy**4-(yy**2)*(zz**2)+zz**4+(xx**2)*(11.*yy**2-zz**2))
                    ylmgr3_cart(1,2,2,2)=ylmcst*(-9.*yy*xx**3+3.*xx*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=-ylmcst*(-12.*zz*yy*xx**2+3.*zz*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,2,3,2)=-ylmcst*3.*zz*(xx**3+xx*(-4.*yy**2+zz**2))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*3.*(-3.*yy*xx**3+xx*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=ylmcst*3.*(xx**4-4.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-4*yy**2+2.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=-ylmcst*3.*zz*(xx**3+xx*(-4.*yy**2+zz**2))
                    ylmgr3_cart(2,2,3,2)=-ylmcst*3.*zz*(3.*yy*xx**2-2.*yy**3+3.*yy*zz**2)
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*(-3.*yy*xx**3-3.*xx*(yy**3-4.*yy*zz**2))
                    ylmgr3_cart(2,3,3,2)=ylmcst*(xx**4-2.*yy**4+11.*(yy**2)*(zz**2)-2.*zz**4-(xx**2)*(yy**2+zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*(-2.*xx**4+yy**4-(yy**2)*(zz**2)-2.*zz**4-(xx**2)*(yy**2-11.*zz**2))
                    ylmgr3_cart(3,3,1,2)=ylmcst*(-3.*yy*xx**3-3.*xx*(yy**3-4.*yy*zz**2))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=-ylmcst*3.*xx*zz*(3.*xx**2+3.*yy**2-2.*zz**2)
                    ylmgr3_cart(3,3,3,2)=-ylmcst*3.*yy*zz*(3.*xx**2+3.*yy**2-2.*zz**2)

                 else if (ll==2 .and. mm==1) then

                    ylmcst=sqrt(15./(pi*8.))/(rr**8)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*6.*zz*(xx**4-6.*(xx**2)*(yy**2+zz**2)+(yy**2+zz**2)**2)
                    ylmgr3_cart(1,1,1,2)=ylmcst*6.*zz*(4.*yy*xx**3-4.*xx*yy*(yy**2+zz**2))
                    !xxy
                    ylmgr3_cart(1,1,2,1)=-ylmcst*2.*zz*(-12.*yy*xx**3+12.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,1,2,2)=-ylmcst*2.*zz*(3.*xx**4+3.*yy**4+2.*(yy**2)*(zz**2)-zz**4+2.*(xx**2)*(-9*yy**2+zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+7.*zz**2)+xx*(-3.*yy**4+6.*(yy**2)*(zz**2)+9.*zz**4))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*2.*(3.*yy*xx**4+2.*(xx**2)*(yy**3-9.*yy*zz**2)-(yy**5-2.*(yy**3)*(zz**2)-3*yy*zz**4))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*2.*zz*(-3.*xx**4-3.*yy**4-2.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(9.*yy**2-zz**2))
                    ylmgr3_cart(1,2,2,2)=ylmcst*2.*zz*(-12.*yy*xx**3+12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*2.*(-3.*yy*xx**4-2.*(xx**2)*(yy**3-9*yy*zz**2)+yy**5-2.*(yy**3)*(zz**2)-3.*yy*zz**4)
                    ylmgr3_cart(1,2,3,2)=ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+zz**2)-3*xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*6.*zz*(-4.*yy*xx**3+4.*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=ylmcst*6.*zz*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-6.*yy**2+2.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+zz**2)-3.*xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(2,2,3,2)=ylmcst*2.*(3.*yy*xx**4+2.*(xx**2)*(yy**3-3.*yy*zz**2)-(yy**5-14*(yy**3)*(zz**2)+9.*yy*zz**4))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz  
                    ylmgr3_cart(2,3,3,1)=ylmcst*2.*zz*(-12.*yy*xx**3-12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,3,3,2)=ylmcst*2.*zz*(3.*xx**4-9.*yy**4+14.*(yy**2)*(zz**2)-zz**4+(xx**2)*(-6.*yy**2+2.*zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=-ylmcst*2.*zz*(9.*xx**4-3.*yy**4-2.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(3.*yy**2-7.*zz**2))
                    ylmgr3_cart(3,3,1,2)=-ylmcst*2.*zz*(12.*yy*xx**3+12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=ylmcst*6.*xx*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(yy**2-3.*zz**2))
                    ylmgr3_cart(3,3,3,2)=ylmcst*6.*yy*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(yy**2-3.*zz**2))

                 else if (ll==2 .and. mm==2) then

                    ylmcst=sqrt(15./(pi*32.))/(rr**8)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=-ylmcst*12.*(-2.*(xx**3)*(2.*yy**2+zz**2)+2.*xx*(2.*yy**4+3.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(1,1,1,2)=-ylmcst*12.*(yy*xx**4-6.*yy*(xx**2)*(yy**2+zz**2)+yy*(yy**2+zz**2)**2)
                    !xxy
                    ylmgr3_cart(1,1,2,1)=ylmcst*4.*(-6.*yy*xx**4-2.*(yy**3)*(yy**2+zz**2)+2.*(xx**2)*(8.*yy**3+3.*yy*zz**2))
                    ylmgr3_cart(1,1,2,2)=ylmcst*4.*(xx**5-2.*(xx**3)*(7.*yy**2+zz**2)+xx*(9.*yy**4+6.*(yy**2)*(zz**2)-3.*zz**4))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*4.*zz*(3.*xx**4+3.*yy**4+4.*(yy**2)*(zz**2)+zz**4-2.*(xx**2)*(9.*yy**2+4.*zz**2))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*4.*zz*(12.*yy*xx**3-12.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*4.*(2.*xx**5+2.*(xx**3)*(-8.*yy**2+zz**2)+6.*xx*(yy**4-(yy**2)*(zz**2)))
                    ylmgr3_cart(1,2,2,2)=ylmcst*4.*(9.*yy*xx**4-2.*(xx**2)*(7.*yy**3-3.*yy*zz**2)+yy**5-2.*(yy**3)*(zz**2)-3.*yy*zz**4)
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*4.*zz*(-12.*yy*xx**3+12.*xx*yy**3)
                    ylmgr3_cart(1,2,3,2)=ylmcst*4.*zz*(3.*xx**4+3.*yy**4+2.*(yy**2)*(zz**2)-zz**4+2.*(xx**2)*(-9.*yy**2+zz**2))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=-ylmcst*12.*(-4.*yy*xx**4+2.*yy*(zz**2)*(yy**2-zz**2)+2.*(xx**2)*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=-ylmcst*12.*(xx**5+(xx**3)*(-6.*yy**2+2.*zz**2)+xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*4.*zz*(3.*xx**4+3.*yy**4-8.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-18*yy**2+4.*zz**2))
                    ylmgr3_cart(2,2,3,2)=ylmcst*4.*zz*(12*yy*xx**3-12.*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*4.*(yy*(xx**4-yy**4-2.*(xx**2)*(zz**2)+8.*(yy**2)*(zz**2)-3.*zz**4)+xx*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(2,3,3,2)=ylmcst*4.*(-xx*(xx**4-yy**4-2.*(xx**2)*(zz**2)+8.*(yy**2)*(zz**2)-3.*zz**4)+yy*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*4.*(xx*(xx**4-yy**4-8.*(xx**2)*(zz**2)+2.*(yy**2)*(zz**2)+3.*zz**4)+yy*(-2.*yy*xx**3-2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,3,1,2)=ylmcst*4.*(yy*(xx**4-yy**4-8.*(xx**2)*(zz**2)+2.*(yy**2)*(zz**2)+3.*zz**4)+xx*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=ylmcst*24.*zz*xx*(xx**2+yy**2-zz**2)
                    ylmgr3_cart(3,3,3,2)=ylmcst*24.*zz*yy*(xx**2+yy**2-zz**2)
                 else
                    write(*,*) 'Illegal mm and ll in initylmg. Remember, ll=3 not implemented yet.'
                    write(*,*) 'll',ll,'mm',mm
                    stop
                 end if

                 ! 2-Transfer gradients into reduced coordinates
                 !ylmgr3_tmp=zero
                 !do aa=1,2
                 !   do ii=1,3
                 !      do jj=1,3
                 !         do kk=1,3
                 !            do pp=1,3
                 !               do qq=1,3
                 !                  do ss=1,3
                 !                     ylmgr3_tmp(ii,jj,kk)=ylmgr3_tmp(ii,jj,kk) &
                 !                          & +rprimd(pp,ii)*rprimd(qq,jj)*rprimd(ss,kk)*ylmgr3_cart(pp,qq,ss,aa)
                 !                  end do
                 !               end do
                 !            end do
                 !         end do
                 !      end do
                 !   end do
                 do ii=1,10
                    ia=alpha3(ii);ib=beta3(ii);igam=gamma3(ii)

                    ylm_gr(ikg+ipw,9+ii,l0+mm)=cmplx(ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2))
                    ylm_gr(ikg+ipw,9+ii,l0-mm)=-conjg(ylm_gr(ikg+ipw,9+ii,l0+mm))

                    !ylm_gr(ikg+ipw,9+ii,l0+mm)=cmplx(ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2))
                    !ylm_gr(ikg+ipw,9+ii,l0-mm)=conjg(cmplx(-1*ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2)))      
                    !ylmgr_red(ii,aa)=ylmgr3_tmp(ia,ib,igam)
                 end do
                 if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
              end do !mm
                 !ylm_gr(ikg+ipw,10:19,l0+mm)=cmplx(ylmgr_red(1:10,1),ylmgr_red(1:10,2))
                 !ylm_gr(ikg+ipw,10:19,l0-mm)=onem*conjg(cmplx(ylmgr_red(1:10,1),ylmgr_red(1:10,2)))

                 !TEST: Check for NaN
!                 do ii=10,19
!                    if (real(ylm_gr(ikg+ipw,ii,l0+mm)) /= real(ylm_gr(ikg+ipw,ii,l0+mm))) then
                       !write(*,*) 'error: ylm_gr is NaN'
                       !write(*,*) 'ii',ii,'mm:',mm,'l0:',l0,'kg:',kg_k(:,ipw)
                      ! stop
!                    end if
                 !if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)              
          ! end do ! loop over m
           end do ! loop over l

        end if !optder=8

!        End condition r<>0
     end if
      
!      End loop over k+G
    end do

!    End condition l<>0
   end if

   !************************************************************* 
   !TEST
!   open (unit=15, file='ylmgr.test', status='replace')
!   write(*,*) kptns(:,1)
!   do ipw=1,npw_k

      !      Load k+G 
!      kpg(1)=kptns(1,ikpt)+real(kg_k(1,ipw),dp)
!      kpg(2)=kptns(2,ikpt)+real(kg_k(2,ipw),dp)
!      kpg(3)=kptns(3,ikpt)+real(kg_k(3,ipw),dp)


      !      Calculate mod of k+G
!      xx=two_pi*(gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3))
!      yy=two_pi*(gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3))
!      zz=two_pi*(gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3))
!      rr=sqrt(xx**2+yy**2+zz**2)

      ! TEST: Use real spherical harmonics for p
      
      !Ylm:
!      ylm(ikg+ipw,2)=cmplx(sqrt(3./four_pi)*yy/rr,zero)
!      ylm(ikg+ipw,3)=cmplx(sqrt(3./four_pi)*zz/rr,zero)
!      ylm(ikg+ipw,4)=cmplx(sqrt(3./four_pi)*xx/rr,zero)
      
      !dYlm/dKz
!      ylm_gr(ikg+ipw,3,2)=cmplx(-sqrt(3./four_pi)*yy*zz/rr**3,zero)
!      ylm_gr(ikg+ipw,3,3)=cmplx(sqrt(3./four_pi)*(xx**2+yy**2)/rr**3,zero)
!      ylm_gr(ikg+ipw,3,4)=cmplx(-sqrt(3./four_pi)*xx*zz/rr**3,zero)

      !d2Ylm/d2Kz
!      ylm_gr(ikg+ipw,9,2)=cmplx(-sqrt(3./four_pi)*yy*(xx**2+yy**2-2*zz**2)/rr**5,zero)
!      ylm_gr(ikg+ipw,9,3)=cmplx(-3.*sqrt(3./four_pi)*zz*(xx**2+yy**2)/rr**5,zero)
!      ylm_gr(ikg+ipw,9,4)=cmplx(-sqrt(3./four_pi)*xx*(xx**2+yy**2-2*zz**2)/rr**5,zero)

      !d3Ylm/d3Kz
!      ylm_gr(ikg+ipw,19,2)=cmplx(3.*sqrt(3./four_pi)*yy*zz*(3.*(xx**2+yy**2)-2*zz**2)/rr**7,zero)
!      ylm_gr(ikg+ipw,19,3)=cmplx(-3.*sqrt(3./four_pi)*(xx**2+yy**2)*(xx**2+yy**2-4.*zz**2)/rr**7,zero)
!      ylm_gr(ikg+ipw,19,4)=cmplx(3.*sqrt(3./four_pi)*xx*zz*(3.*(xx**2+yy**2)-2*zz**2)/rr**7,zero)

      


!      write(15,'(9e16.4e2)') xx,yy,zz,ylm(ipw,2),ylm(ipw,3),ylm(ipw,4)!,ylm_gr(ipw,19,2),ylm_gr(ipw,19,3),ylm_gr(ipw,19,4)


!   end do
!   close(unit=15)
!   stop
   !*************************************************************

   ABI_FREE(kg_k)

   ikg=ikg+npw_k
 end do !  End Loop over k-points

!Release the temporary memory
!Allocate some memory
!TEST
!write(*,*)'ylm(20,2):',ylm(57,4)
!write(*,*)'ylm_gr(20,1,2):',ylm_gr(57,1,4)
!write(*,*)'ylm_gr(20,4,2):',ylm_gr(57,4,4)
!write(*,*)'ylm_gr(20,10,2):',ylm_gr(57,10,4)
!stop

!write(*,*) 'here2'
! if (optder/=0) then
!   ABI_FREE(ylmgr_cart)
! end if
!TEST
!write(*,*) 'here3'
! if (optder/=0.and.optder/=2) then
!   ABI_FREE(ylmgr_red)
! end if
!TEST
!write(*,*) 'here4'
! if (optder>=2) then
!   ABI_FREE(ylmgr2_cart)
!   ABI_FREE(ylmgr2_tmp)
!   ABI_FREE(ylmgr_red)
!   ABI_FREE(blm)
! end if

end subroutine cinitylmgi
!!***





end module m_dfpt_indpol
