!!****m* ABINIT/m_eph_driver
!! NAME
!!  m_eph_driver
!!
!! FUNCTION
!!   Driver for EPH calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2022 ABINIT group (MG, MVer, GA)
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

module m_eph_driver

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_crystal
 use m_ebands
 use m_dtset
 use m_efmas_defs
 use m_dtfil
 use m_ddb
 use m_ddb_hdr
 use m_dvdb
 use m_ifc
 use m_phonons
 use m_nctk
 use m_wfk
 use m_distribfft
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : MPI_type
 use m_io_tools,        only : file_exists, open_file
 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : strcat, sjoin, ftoa, itoa
 use m_fftcore,         only : print_ngfft
 use m_frohlichmodel,   only : frohlichmodel,polaronmass
 use m_rta,             only : rta_driver, ibte_driver
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_paw_an,          only : paw_an_type, paw_an_free !, paw_an_nullify, paw_an_init,
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_symrhoij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_phgamma,         only : eph_phgamma
 use m_efmas,           only : efmasdeg_free_array, efmasval_free_array, efmas_ncread
 use m_gkk,             only : eph_gkk, ncwrite_v1qnu
 use m_phpi,            only : eph_phpi
 use m_sigmaph,         only : sigmaph
 use m_pspini,          only : pspini
 use m_ephtk,           only : ephtk_update_ebands

 implicit none

 private
!!***

 public :: eph
!!***

contains
!!***

!!****f* m_eph_driver/eph
!! NAME
!!  eph
!!
!! FUNCTION
!! Main routine to compute electron phonon coupling matrix elements and
!! calculate related properties - superconductin Tc, phonon linewidths, electronic renormalization
!! due to phonons and temperature effects...
!!
!! INPUTS
!! acell(3)=Length scales of primitive translations (bohr)
!! codvsn=Code version
!! dtfil<datafiles_type>=Variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!!   Before entering the first time in the routine, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=Dimensionless real space primitive translations.
!! xred(3,natom)=Reduced atomic coordinates.
!!
!! PARENTS
!!      m_driver
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! CHILDREN
!!      cryst%free,cryst%print,cryst_ddb%free,cwtime,cwtime_report,ddb%free
!!      ddb_from_file,ddb_hdr%free,destroy_mpi_enreg,dtset%free_nkpt_arrays
!!      dvdb%free,dvdb%interpolate_and_write,dvdb%list_perts,dvdb%open_read
!!      dvdb%print,dvdb%write_v1qavg,ebands_free,ebands_prtbltztrp,ebands_write
!!      efmas_ncread,efmasdeg_free_array,efmasval_free_array,eph_gkk
!!      eph_phgamma,eph_phpi,ephtk_update_ebands,frohlichmodel,ibte_driver
!!      ifc%free,ifc%outphbtrap,ifc%print,ifc%printbxsf,ifc_init,ifc_mkphbs
!!      init_distribfft_seq,initmpi_seq,mkphdos,ncwrite_v1qnu,pawfgr_destroy
!!      pawfgr_init,phdos%free,phdos%ncwrite,phdos%print,print_ngfft,pspini
!!      rta_driver,sigmaph,test_phrotation,wfk0_hdr%free,wfk0_hdr%vs_dtset
!!      wfq_hdr%free,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine eph(acell, codvsn, dtfil, dtset, pawang, pawrad, pawtab, psps, rprim, xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, natifc0 = 0, selectz0 = 0, nsphere0 = 0, prtsrlr0 = 0
 integer :: ii,comm,nprocs,my_rank,psp_gencond,mgfftf,nfftf
 integer :: iblock_dielt_zeff, iblock_dielt, iblock_quadrupoles, ddb_nqshift, ierr, npert_miss
 integer :: omp_ncpus, work_size, nks_per_proc, mtyp, mpert, lwsym !msize,
 integer :: iatdir, iq2dir, iq1dir, quad_unt, iatom, jj
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff
 real(dp) :: cpu,wall,gflops
 logical :: use_wfk, use_wfq, use_dvdb
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path, wfq_path, ddb_filepath, dvdb_filepath, path
 type(hdr_type) :: wfk0_hdr, wfq_hdr
 type(crystal_t) :: cryst, cryst_ddb
 type(ebands_t) :: ebands, ebands_kq
 type(ddb_type) :: ddb, ddb_lw
 type(ddb_hdr_type) :: ddb_hdr
 type(dvdb_t) :: dvdb
 type(ifc_type) :: ifc
 type(pawfgr_type) :: pawfgr
 type(mpi_type) :: mpi_enreg
 type(phonon_dos_type) :: phdos
!arrays
 integer :: ngfftc(18), ngfftf(18), count_wminmax(2)
 integer,allocatable :: dummy_atifc(:)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: wminmax(2), dielt(3,3), zeff(3,3,dtset%natom), zeff_raw(3,3,dtset%natom)
 real(dp) :: qdrp_cart(3,3,3,dtset%natom)
 real(dp),allocatable :: ddb_qshifts(:,:), kpt_efmas(:,:)
 type(efmasdeg_type),allocatable :: efmasdeg(:)
 type(efmasval_type),allocatable :: efmasval(:,:)
 !type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)

!************************************************************************

 ! This part performs the initialization of the basic objects used to perform e-ph calculations:
 !
 !     1) Crystal structure `cryst`
 !     2) Ground state band energies: `ebands`
 !     3) Interatomic force constants: `ifc`
 !     4) DVDB database with the dvscf potentials
 !     5) Pseudos and PAW basic objects.
 !
 ! Once we have these objects, we can call specialized routines for e-ph calculations.
 ! Notes:
 !
 !   * Any modification to the basic objects mentioned above should be done here (e.g. change of efermi)
 !   * This routines shall not allocate big chunks of memory. The CPU-demanding sections should be
 !     performed in the subdriver that will employ different MPI distribution schemes optimized for that particular task.

 DBG_ENTER('COLL')

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

#ifndef HAVE_MPI_IBCAST
 do ii=1,5
   ABI_WARNING("Your MPI library does not provide MPI_IBCAST. Calculations parallelized over perturbations will be slow")
 end do
#endif

 ! Initialize filenames
 wfk0_path = dtfil%fnamewffk
 wfq_path = dtfil%fnamewffq
 ddb_filepath = dtfil%filddbsin

 ! Use the ddb file as prefix if getdvdb or irddvb are not given in the input.
 dvdb_filepath = dtfil%fildvdbin
 if (dvdb_filepath == ABI_NOFILE) then
   dvdb_filepath = dtfil%filddbsin; ii = len_trim(dvdb_filepath); dvdb_filepath(ii-2:ii+1) = "DVDB"
 end if

 use_wfk = all(dtset%eph_task /= [0, 5, -5, 6, +15, -15, -16, 16])
 use_wfq = ((dtset%irdwfq /= 0 .or. dtset%getwfq /= 0 .or. dtset%getwfq_filepath /= ABI_NOFILE) .and. dtset%eph_frohlichm /= 1)

 ! If eph_task is needed and ird/get variables are not provided, assume WFQ == WFK
 if (any(dtset%eph_task == [2, -2, 3]) .and. .not. use_wfq) then
   wfq_path = wfk0_path
   use_wfq = .True.
   write(msg, "(4a)")&
     "eph_task requires WFQ but none among (irdwfq, getwfq, getwfk_filepath) is specified in input.", ch10, &
     "Will read WFQ wavefunctions from WFK file:", trim(wfk0_path)
   ABI_COMMENT(msg)
 end if

 use_dvdb = (dtset%eph_task /= 0 .and. dtset%eph_frohlichm /= 1 .and. abs(dtset%eph_task) /= 7)

 if (my_rank == master) then
   if (.not. file_exists(ddb_filepath)) ABI_ERROR(sjoin("Cannot find DDB file:", ddb_filepath))
   if (use_dvdb .and. .not. file_exists(dvdb_filepath)) ABI_ERROR(sjoin("Cannot find DVDB file:", dvdb_filepath))

   ! Accept WFK file in Fortran or netcdf format.
   if (use_wfk .and. nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find GS WFK file:", wfk0_path, ". Error:", msg))
   end if
   ! WFQ file
   if (use_wfq) then
     if (nctk_try_fort_or_ncfile(wfq_path, msg) /= 0) then
       ABI_ERROR(sjoin("Cannot find GS WFQ file:", wfq_path, ". Error:", msg))
     end if
   end if
 end if ! master

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 if (use_wfk) then
   call xmpi_bcast(wfk0_path, master, comm, ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))
 end if
 if (use_wfq) then
   call xmpi_bcast(wfq_path, master, comm, ierr)
   call wrtout(ab_out, sjoin("- Reading GS states from WFQ file:", wfq_path) )
 end if
 call wrtout(ab_out, sjoin("- Reading DDB from file:", ddb_filepath))
 if (use_dvdb) call wrtout(ab_out, sjoin("- Reading DVDB from file:", dvdb_filepath))
 if (dtset%eph_frohlichm /= 0) call wrtout(ab_out, sjoin("- Reading EFMAS information from file:", dtfil%fnameabi_efmas))
 call wrtout(ab_out, ch10//ch10)

 ! autoparal section
 ! TODO: This just to activate autoparal in AbiPy. Lot of things should be improved.
 if (dtset%max_ncpus /= 0) then
   write(ab_out,'(a)')"--- !Autoparal"
   write(ab_out,"(a)")"# Autoparal section for EPH runs"
   write(ab_out,"(a)")   "info:"
   write(ab_out,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ab_out,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ab_out,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ab_out,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ab_out,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ab_out,"(a,i0)")"    mband: ",dtset%mband
   write(ab_out,"(a,i0)")"    eph_task: ",dtset%eph_task

   work_size = dtset%nkpt * dtset%nsppol
   ! Non-scalable memory in Mb i.e. memory that is not distributed with MPI.
   nonscal_mem = zero
   max_wfsmem_mb = (two * dp * dtset%mpw * dtset%mband * dtset%nkpt * dtset%nsppol * dtset%nspinor * b2Mb) * 1.1_dp

   ! List of configurations.
   ! Assuming an OpenMP implementation with perfect speedup!
   write(ab_out,"(a)")"configurations:"

   do ii=1,dtset%max_ncpus
     nks_per_proc = work_size / ii
     nks_per_proc = nks_per_proc + mod(work_size, ii)
     eff = (one * work_size) / (ii * nks_per_proc)
     ! Add the non-scalable part and increase by 10% to account for other datastructures.
     mempercpu_mb = (max_wfsmem_mb + nonscal_mem) * 1.1_dp

     do omp_ncpus=1,1 !xomp_get_max_threads()
       write(ab_out,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
       write(ab_out,"(a,i0)")"      mpi_ncpus: ",ii
       write(ab_out,"(a,i0)")"      omp_ncpus: ",omp_ncpus
       write(ab_out,"(a,f12.9)")"      efficiency: ",eff
       write(ab_out,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
     end do
   end do
   write(ab_out,'(a)')"..."
   ABI_STOP("Stopping now!")
 end if

 call cwtime(cpu, wall, gflops, "start")

 if (use_wfk) then
   ! Construct crystal and ebands from the GS WFK file.
   ebands = wfk_read_ebands(wfk0_path, comm, out_hdr=wfk0_hdr)
   call wfk0_hdr%vs_dtset(dtset)

   cryst = wfk0_hdr%get_crystal()
   call cryst%print(header="crystal structure from WFK file")

   ! Here we change the GS bands (Fermi level, scissors operator ...)
   ! All the modifications to ebands should be done here.
   call ephtk_update_ebands(dtset, ebands, "Ground state energies")
 end if

 if (use_wfq) then
   ! Read WFQ and construct ebands on the shifted grid.
   ebands_kq = wfk_read_ebands(wfq_path, comm, out_hdr=wfq_hdr)
   ! GKA TODO: Have to construct a header with the proper set of q-shifted k-points then compare against dtset.
   !call wfq_hdr%vs_dtset(dtset)
   call wfq_hdr%free()
   call ephtk_update_ebands(dtset, ebands_kq, "Ground state energies (K+Q)")
 end if

 call cwtime_report(" eph%init", cpu, wall, gflops)

 ! =======================================
 ! Output useful info on electronic bands
 ! =======================================
 call cwtime(cpu, wall, gflops, "start")

 if (my_rank == master) then
   ! Fermi Surface
   if (dtset%prtfsurf /= 0) then
     path = strcat(dtfil%filnam_ds(4), "_BXSF")
     call wrtout(ab_out, sjoin("- Writing Fermi surface to file:", path))
     if (ebands_write_bxsf(ebands, cryst, path) /= 0) then
       msg = "Cannot produce file for Fermi surface, check log file for more info"
       ABI_WARNING(msg)
       call wrtout(ab_out, msg)
     end if
   end if

   ! Nesting factor (requires qpath)
   if (dtset%prtnest /= 0 .and. dtset%ph_nqpath > 0) then
     path = strcat(dtfil%filnam_ds(4), "_NEST")
     call wrtout(ab_out, sjoin("- Writing nesting factor to file:", path))
     if (ebands_write_nesting(ebands, cryst, path, dtset%prtnest, &
         dtset%tsmear, dtset%fermie_nest, dtset%ph_qpath(:,1:dtset%ph_nqpath), msg) /= 0) then
       ABI_WARNING(msg)
       call wrtout(ab_out,msg)
     end if
   end if

   if (use_wfk) call ebands_write(ebands, dtset%prtebands, dtfil%filnam_ds(4))
 end if

 call cwtime_report(" eph%ebands_postprocess:", cpu, wall, gflops)

 ! Read the DDB file.
 ABI_CALLOC(dummy_atifc, (dtset%natom))

 if (use_wfk) then
   call ddb_from_file(ddb, ddb_filepath, dtset%brav, dtset%natom, natifc0, dummy_atifc, ddb_hdr, cryst_ddb, comm, &
                      prtvol=dtset%prtvol)

   ! DDB cryst comes from DPPT --> no time-reversal if q /= 0
   ! Change the value so that we use the same as the GS part.
   cryst_ddb%timrev = cryst%timrev
   if (cryst%compare(cryst_ddb, header=" Comparing WFK crystal with DDB crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and DDB do not agree! Check messages above!")
   end if
   call cryst_ddb%free()
 else
   ! Get crystal from DDB.
   ! Warning: We may loose precision in rprimd and xred because DDB does not have enough significant digits.
   call ddb_from_file(ddb, ddb_filepath, dtset%brav, dtset%natom, natifc0, dummy_atifc, ddb_hdr, cryst, comm, &
                      prtvol=dtset%prtvol)
 end if

 ABI_FREE(dummy_atifc)

 ! Set the q-shift for the DDB (well we mainly use gamma-centered q-meshes)
 ddb_nqshift = 1
 ABI_CALLOC(ddb_qshifts, (3, ddb_nqshift))
 ddb_qshifts(:,1) = dtset%ddb_shiftq(:)


 ! Get Dielectric Tensor
 iblock_dielt = ddb%get_dielt(dtset%rfmeth, dielt)

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblock_dielt_zeff = ddb%get_dielt_zeff(cryst, dtset%rfmeth, dtset%chneut, selectz0, dielt, zeff, zeff_raw=zeff_raw)
 if (my_rank == master) then
   if (iblock_dielt_zeff == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_filepath))
     call wrtout(ab_out, " Values initialized with zeros.")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_filepath))
   end if
 end if

 ! Read the quadrupoles
 !iblock_quadrupoles = ddb%get_quadrupoles(1, 3, qdrp_cart)

 iblock_quadrupoles = 0
 qdrp_cart = zero

 if (my_rank == master) then
   ! MG: Temporary hack to read the quadrupole tensor from a text file
   ! Will be removed when the EPH code will be able to read Q* from the DDB.

   if (file_exists("quadrupoles_cart.out")) then
     call wrtout(std_out, " Reading quadrupoles from quadrupoles_cart.out")
     quad_unt = 71
     open(unit=quad_unt,file="quadrupoles_cart.out",action="read")
     do ii=1,2
       read(quad_unt,*) msg
       write(std_out, *)" msg: ", trim(msg)
     end do

     do ii=1,3
       do jj=1,3*3*ddb%natom
         read(quad_unt,'(4(i5,3x),2(1x,f20.10))') iq2dir,iatom,iatdir,iq1dir,qdrp_cart(iq1dir,iq2dir,iatdir,iatom)
         write(std_out, *) iq2dir,iatom,iatdir,iq1dir,qdrp_cart(iq1dir,iq2dir,iatdir,iatom)
       end do
       read(quad_unt,'(a)') msg
     end do
     close(quad_unt)
     iblock_quadrupoles = 1
   end if
 end if
 call xmpi_bcast(iblock_quadrupoles, master, comm, ierr)
 call xmpi_bcast(qdrp_cart, master, comm, ierr)

 ! Here we get the quadrupoles from the DDB file (this should become the official API).
 ! Section Copied from Anaddb.

if (iblock_quadrupoles == 0) then
 mtyp = ddb_hdr%mblktyp
 mpert = dtset%natom + MPERT_MAX
 !msize = 3*mpert*3*mpert; if (mtyp==3) msize=msize*3*mpert
 !call wrtout(std_out, sjoin(" Trying to read Q* from DDB file, mtyp:", itoa(mtyp)))

 ! MR: a new ddb is necessary for the longwave quantities due to incompability of it with authomatic reshapes
 ! that ddb%val and ddb%flg experience when passed as arguments of some routines

 ! Get Quadrupole tensor
 qdrp_cart=zero
 if (mtyp==33) then
   lwsym=1
   call ddb_lw_copy(ddb, ddb_lw, mpert, dtset%natom, dtset%ntypat)
   iblock_quadrupoles = ddb_lw%get_quadrupoles(lwsym, 33, qdrp_cart)
   call ddb_lw%free()
 end if
endif

 ! The default value is 1. Here we set the flags to zero if Q* is not available.
 if (iblock_quadrupoles == 0) then
   dtset%dipquad = 0
   dtset%quadquad = 0
 end if

 call ddb_hdr%free()

 if (my_rank == master) then
   if (iblock_quadrupoles == 0) then
     call wrtout(ab_out, sjoin("- Cannot find quadrupole tensor in DDB file:", ddb_filepath))
     call wrtout(ab_out, " Values initialized with zeros.")
   else
     call wrtout(ab_out, sjoin("- Found quadrupole tensor in DDB file:", ddb_filepath))
   end if
 end if

 call ifc_init(ifc, cryst, ddb, &
   dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, &
   dtset%ddb_ngqpt, ddb_nqshift, ddb_qshifts, dielt, zeff, &
   qdrp_cart, nsphere0, dtset%rifcsph, prtsrlr0, dtset%enunit, comm, &
   dipquad=dtset%dipquad, quadquad=dtset%quadquad)

 ABI_FREE(ddb_qshifts)
 if (my_rank == master) then
   call ifc%print(unit=std_out)
   !call ifc%print(unit=ab_out)
 end if

 ! Output phonon band structure (requires qpath)
 if (dtset%prtphbands /= 0) call ifc_mkphbs(ifc, cryst, dtset, dtfil%filnam_ds(4), comm)

 if (dtset%prtphdos == 1) then
   call wrtout(std_out, " Computing Phonon DOS. Use prtphdos 0 to disable this part.")
   wminmax = zero
   do
     call mkphdos(phdos, cryst, ifc, dtset%ph_intmeth, dtset%ph_wstep, dtset%ph_smear, dtset%ph_ngqpt, &
       dtset%ph_nqshift, dtset%ph_qshift, "", wminmax, count_wminmax, comm)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(msg, "(a, 2f8.5)") "Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
   end do

   if (my_rank == master) then
     if (dtset%prtvol > 0) then
       ! Disabled by default because it's slow and we use netcdf that is much better.
       path = strcat(dtfil%filnam_ds(4), "_PHDOS")
       call wrtout(ab_out, sjoin("- Writing phonon DOS to file:", path))
       call phdos%print(path)
     end if

#ifdef HAVE_NETCDF
     path = strcat(dtfil%filnam_ds(4), "_PHDOS.nc")
     call wrtout(ab_out, sjoin("- Writing phonon DOS to netcdf file:", path))
     ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, sjoin("Creating PHDOS.nc file:", path))
     NCF_CHECK(cryst%ncwrite(ncid))
     call phdos%ncwrite(ncid)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if
   call phdos%free()
 end if ! prtphdos

 if (dtset%prtbltztrp == 1 .and. my_rank == master) then
   call ifc%outphbtrap(cryst, dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, dtfil%filnam_ds(4))
   ! BoltzTraP output files in GENEric format
   call ebands_prtbltztrp(ebands, cryst, dtfil%filnam_ds(4))
 end if

 ! Output phonon isosurface in Xcrysden format.
 if (dtset%prtphsurf == 1) then
   path = strcat(dtfil%filnam_ds(4), "_PH.bxsf")
   call wrtout(ab_out, sjoin("- Writing phonon frequencies in Xcrysden format to file:", path))
   call ifc%printbxsf(cryst, dtset%ph_ngqpt, dtset%ph_nqshift, dtset%ph_qshift, path, comm)
 end if

 call cwtime_report(" eph%ifc:", cpu, wall, gflops)

 ! Initialize the object used to read DeltaVscf (required if eph_task /= 0)
 if (use_dvdb) then
   dvdb = dvdb_new(dvdb_filepath, comm)

   ! DVDB cryst comes from DPPT --> no time-reversal if q /= 0
   ! Change the value so that we use the same as the GS part.
   dvdb%cryst%timrev = cryst%timrev
   if (cryst%compare(dvdb%cryst, header=" Comparing WFK crystal with DVDB crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and DVDB do not agree! Check messages above!")
   end if
   if (dtset%prtvol > 10) dvdb%debug = .True.

   ! This to symmetrize the DFPT potentials.
   dvdb%symv1 = dtset%symv1scf

   ! Copy brav variable
   dvdb%brav = dtset%brav

   ! Select algorithm for generating the list of R-points and the weigths used to compute W(r,R)
   dvdb%rspace_cell = dtset%dvdb_rspace_cell

   !call dvdb%load_ddb(dtset%prtvol, comm, ddb=ddb)

   ! Set qdamp from frohl_params
   dvdb%qdamp = dtset%dvdb_qdamp

   ! Set quadrupoles
   dvdb%qstar = qdrp_cart; if (iblock_quadrupoles /= 0) dvdb%has_quadrupoles = .True.

   ! Set dielectric tensor, BECS and associated flags.
   ! This flag activates automatically the treatment of the long-range term in the Fourier interpolation
   ! of the DFPT potentials except when dvdb_add_lr == 0
   dvdb%add_lr = dtset%dvdb_add_lr
   if (iblock_dielt /= 0) then
     dvdb%has_dielt = .True.; dvdb%dielt = dielt
   end if
   if (iblock_dielt_zeff /= 0) then
     dvdb%has_zeff = .True.; dvdb%zeff = zeff; dvdb%zeff_raw = zeff_raw
   end if
   if (.not. dvdb%has_dielt .or. .not. (dvdb%has_zeff .or. dvdb%has_quadrupoles)) then
     if (dvdb%add_lr /= 0) then
       dvdb%add_lr = 0
       ABI_WARNING("Setting dvdb_add_lr to 0. Long-range term won't be substracted in Fourier interpolation.")
     end if
   end if

   if (dvdb%add_lr == 2) then
      if (dvdb%has_quadrupoles) then
        call wrtout(std_out, "dvdb_add_lr == 2 --> Quadrupoles are set to zero and won't be used in the interpolation")
      end if
      dvdb%has_quadrupoles = .False.
      dvdb%qstar = zero
   end if

   if (my_rank == master) then
     call dvdb%print()
     call dvdb%list_perts([-1, -1, -1], npert_miss)
     ABI_CHECK(npert_miss == 0, sjoin(itoa(npert_miss), "independent perturbation(s) are missing in the DVDB file!"))
   end if
 end if

 call pawfgr_init(pawfgr, dtset, mgfftf, nfftf, ecut_eff, ecutdg_eff, ngfftc, ngfftf, &
                  gsqcutc_eff=gsqcutc_eff, gsqcutf_eff=gsqcutf_eff, gmet=cryst%gmet, k0=k0)

 call print_ngfft(ngfftc, header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft(ngfftf, header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 call init_distribfft_seq(mpi_enreg%distribfft, 'f', ngfftf(2), ngfftf(3), 'all')

 ! I am not sure yet the EFMAS file will be needed as soon as eph_frohlichm/=0. To be decided later.
 if (dtset%eph_frohlichm /= 0) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(ncid, dtfil%fnameabi_efmas, xmpi_comm_self))
   call efmas_ncread(efmasdeg, efmasval, kpt_efmas, ncid)
   NCF_CHECK(nf90_close(ncid))
#else
   ABI_ERROR("netcdf support not enabled")
#endif
 end if

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset, dtfil, ecore, psp_gencond, gsqcutc_eff, gsqcutf_eff, pawrad, pawtab, psps, cryst%rprimd, comm_mpi=comm)

 ! TODO: Make sure that all subdrivers work with useylm == 1
 ABI_CHECK(dtset%useylm == 0, "useylm != 0 not implemented/tested")

 ! Relase nkpt-based arrays in dtset to decreased memory requirement if dense sampling.
 ! EPH routines should not access them after this point.
 if (all(dtset%eph_task /= [6, 10])) call dtset%free_nkpt_arrays()

 ! ====================================================
 ! === This is the real EPH stuff once all is ready ===
 ! ====================================================

 select case (dtset%eph_task)
 case (0)
   ! This is just to access the DDB post-processing tools
   continue

 case (1)
   ! Compute phonon linewidths in metals.
   call eph_phgamma(wfk0_path, dtfil, ngfftc, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                    pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

 case (2, -2)
   ! Compute e-ph matrix elements.
   call eph_gkk(wfk0_path, wfq_path, dtfil, ngfftc, ngfftf, dtset, cryst, ebands, ebands_kq, dvdb, ifc, &
                pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

 case (3)
   ! Compute phonon self-energy.
   call eph_phpi(wfk0_path, wfq_path, dtfil, ngfftc, ngfftf, dtset, cryst, ebands, ebands_kq, dvdb, ifc, &
                 pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

 case (4, -4)
   ! Compute electron self-energy (phonon contribution)
   call sigmaph(wfk0_path, dtfil, ngfftc, ngfftf, dtset, cryst, ebands, dvdb, ifc, wfk0_hdr, &
                pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

   ! Compute transport properties in the RTA/IBTE only if sigma_erange has been used
   if (dtset%eph_task == -4 .and. any(abs(dtset%sigma_erange) > zero)) then
     if (dtset%ibte_prep > 0) then
       call ibte_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm) ! Solve IBTE
     else
       call rta_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm)  ! Compute RTA
     end if
   end if

 case (5, -5)
   ! Interpolate the phonon potential.
   call dvdb%interpolate_and_write(dtset, dtfil%fnameabo_dvdb, ngfftc, ngfftf, cryst, ifc%ngqpt, ifc%nqshft, ifc%qshft, comm)

 case (6)
   ! Estimate zero-point renormalization and temperature-dependent electronic structure using the Frohlich model
   if (my_rank == master) then
     call frohlichmodel(cryst, dtset, efmasdeg, efmasval, ifc)
   end if

 case (7)
   ! Compute phonon-limited RTA from SIGEPH file.
   call rta_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm)

 case (8)
   ! Solve IBTE from SIGEPH file.
   call ibte_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm)

 case (10)
   ! Estimate polaron effective mass in the triply-degenerate VB or CB cubic case
   if (my_rank == master) then
     call polaronmass(cryst, dtset, efmasdeg, efmasval, ifc)
   end if

 case (15, -15)
   ! Write average of DFPT potentials to file.
   if (nprocs > 1) then
     ABI_WARNING("eph_task in [15, -15] does not support nprocs > 1. Running in sequential...")
   end if
   dvdb%comm = xmpi_comm_self
   if (my_rank == master) then
     call dvdb%open_read(ngfftf, xmpi_comm_self)
     !call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, comm)
     !call dvdb%set_pert_distrib(sigma%comm_pert, sigma%my_pinfo, sigma%pert_table)
     call dvdb%write_v1qavg(dtset, strcat(dtfil%filnam_ds(4), "_V1QAVG.nc"))
   end if

 case (-16, 16)
   if (nprocs > 1) then
     ABI_WARNING("eph_task in [16, -16] does not support nprocs > 1. Running in sequential...")
   end if

   call test_phrotation(ifc, cryst, dtset%ph_ngqpt, comm)

   dvdb%comm = xmpi_comm_self
   if (my_rank == master) then
     call dvdb%open_read(ngfftf, xmpi_comm_self)

     ! Compute \delta V_{q,nu)(r) and dump results to netcdf file.
     call ncwrite_v1qnu(dvdb, dtset, ifc, strcat(dtfil%filnam_ds(4), "_V1QNU.nc"))
   end if

 case default
   ABI_ERROR(sjoin("Unsupported value of eph_task:", itoa(dtset%eph_task)))
 end select

 !=====================
 !==== Free memory ====
 !=====================
 call cryst%free()
 call dvdb%free()
 call ddb%free()
 call ifc%free()
 call wfk0_hdr%free()
 call ebands_free(ebands)
 call ebands_free(ebands_kq)
 call pawfgr_destroy(pawfgr)
 call destroy_mpi_enreg(mpi_enreg)

 if (allocated(efmasdeg)) call efmasdeg_free_array(efmasdeg)
 if (allocated(efmasval)) call efmasval_free_array(efmasval)
 ABI_SFREE(kpt_efmas)

 ! Deallocation for PAW.
 if (dtset%usepaw == 1) then
   !call pawrhoij_free(pawrhoij)
   !ABI_FREE(pawrhoij)
   !call pawfgrtab_free(pawfgrtab)
   !ABI_FREE(pawfgrtab)
   !call paw_ij_free(paw_ij)
   !ABI_FREE(paw_ij)
   !call paw_an_free(paw_an)
   !ABI_FREE(paw_an)
 end if

 DBG_EXIT('COLL')

end subroutine eph
!!***

end module m_eph_driver
!!***
