MODULE diag_rb
!-------------------------------------------------------------------------------
!
!     Read binary files
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  use diag_header
  use netcdf
  implicit none

  public

  integer, parameter ::  &
      recl_cnt = nhead+DP+(2*nx+1)*(ny+1)*(2*nz)*(2*nv)*(nm+1)*(2*DP)+nhead,  &
      recl_fxv = nhead+DP+(2*nx+1)*(ny+1)*(2*nv)*(nm+1)*(2*DP)+nhead,         &
      recl_phi = nhead+DP+(2*nx+1)*(ny+1)*(2*nz)*(2*DP)+nhead,                &
      recl_Al  = nhead+DP+(2*nx+1)*(ny+1)*(2*nz)*(2*DP)+nhead,                &
      recl_mom = nhead+DP+(2*nx+1)*(ny+1)*(2*nz)*nmom*(2*DP)+nhead,           &
      recl_trn = nhead+DP+(2*nx+1)*(ny+1)*ntrn*DP+nhead,                      &
      recl_tri = nhead+DP+(2*nx+1)*(2*global_ny+1)*ntri*DP+nhead

  integer :: nloop_phi, nloop_Al, nloop_mom, nloop_trn, nloop_tri, &
             nloop_fxv, nloop_cnt  ! Total number of time step loop
  integer, dimension(1:enum) :: loop_phi_sta, loop_phi_end,  &
                                loop_Al_sta,  loop_Al_end,   &
                                loop_mom_sta, loop_mom_end,  &
                                loop_trn_sta, loop_trn_end,  &
                                loop_tri_sta, loop_tri_end,  &
                                loop_fxv_sta, loop_fxv_end,  &
                                loop_cnt_sta, loop_cnt_end

 !character(len=9), parameter :: flag_getsize = "stat"    ! GNU Fortran
  character(len=9), parameter :: flag_getsize = "inquire" ! Intel, Fujitsu Fortran

  character(len=9), save :: calc_type
  integer, save :: num_triad_diag
  integer, dimension(:), allocatable, save :: triad_diag_mxt, triad_diag_myt

  integer(kind=4), dimension(1:enum) :: phi_nc, &
                                Al_nc, &
                                mom_nc, &
                                trn_nc, &
                                tri_nc, &
                                fxv_nc, &
                                cnt_nc

 CONTAINS

SUBROUTINE rb_gmy2rankwmy( gmy, rankw, my )
!-------------------------------------------------------------------------------
!
!     Calculate rankw and my from gmy
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: gmy
  integer, intent(out) :: rankw, my

    rankw = gmy / ( ny+1 )
    !my = mod( gmy, ny+1 )
    my = gmy - ( ny+1 ) * rankw

END SUBROUTINE rb_gmy2rankwmy


SUBROUTINE rb_giz2rankziz( giz, rankz, iz )
!-------------------------------------------------------------------------------
!
!     Calculate rankz and iz from giz
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz
  integer, intent(out) :: rankz, iz

    rankz = ( giz + global_nz ) / ( 2*nz )
    !iz = mod( giz + global_nz, 2*nz ) - nz
    iz = giz + global_nz - 2*nz * rankz - nz

END SUBROUTINE rb_giz2rankziz


SUBROUTINE rb_giv2rankviv( giv, rankv, iv )
!-------------------------------------------------------------------------------
!
!     Calculate rankv and iv from giv
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giv
  integer, intent(out) :: rankv, iv

    rankv = ( giv - 1 ) / ( 2*nv )
    !iv = mod( giv - 1, 2*nv ) + 1
    iv = giv - 2*nv * rankv

END SUBROUTINE rb_giv2rankviv


SUBROUTINE rb_gim2rankmim( gim, rankm, im )
!-------------------------------------------------------------------------------
!
!     Calculate rankm and im from gim
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: gim
  integer, intent(out) :: rankm, im

    rankm = gim / ( nm+1 )
    !im = mod( gim, nm+1 )
    im = gim - (nm+1) * rankm

END SUBROUTINE rb_gim2rankmim


SUBROUTINE rb_phi_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *phi*
!                                                   (S. Maeyama, 28 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, rankz, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  character(len=3) :: cnum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do rankz = 0, nprocz-1
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*rankz
      !    write( crank, fmt="(i6.6)" ) ir
      !    open( unit=100000000+100000*inum+ir,                      &
      !          file="../phi/gkvp."//crank//".0.phi."//cnum,  &
      !          status="old", action="read",                        &
      !          form="unformatted", access="stream" )
                                        !--- for debug ---
                                        !write(*,*) 100000000+100000*inum+ir, &
                                        !           "../phi/gkvp."      &
                                        !           //crank//".0.phi."//cnum
                                        !-----------------
      !  end do
      !end do
      ierr_nf90 = nf90_open( path="../phi/gkvp.phi."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=phi_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_phi_fileopen


SUBROUTINE rb_phi_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *phi*
!                                                   (S. Maeyama, 28 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, rankz, rankw, ir
  integer :: inum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do rankz = 0, nprocz-1
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*rankz
      !    close( unit=100000000+100000*inum+ir )
      !  end do
      !end do

      ierr_nf90 = nf90_close( phi_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_phi_fileclose


SUBROUTINE rb_phi_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 27 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz
  integer(kind=4) :: len_kx, len_ky, len_zz

    do inum = snum, enum
      !ir = 0
      !read( unit=100000000+100000*inum+ir, pos=1 ) head
      !read( unit=100000000+100000*inum+ir, pos=recl_phi-4+1 ) foot

      !if ( head /= foot .or. head /= recl_phi - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *phi*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_phi-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if

      ierr_nf90 = nf90_inq_varid( phi_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( phi_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( phi_nc(inum), "zz", varid_zz )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( phi_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( phi_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( phi_nc(inum), varid_zz, len=len_zz )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           len_zz /= 2*global_nz ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.phi*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#              zz = ", len_zz
        write(*,*) "#     2*global_nz = ", 2*global_nz
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *phi* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_phi-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
    write(*,*) "# The netCDF files gkvp.phi*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#              zz = ", len_zz
    write(*,*) "#     2*global_nz = ", 2*global_nz
                                        !--- for debug ---
                                        !write(*,*) head, recl_phi-2*nhead, foot
                                        !-----------------

END SUBROUTINE rb_phi_check


SUBROUTINE rb_phi_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 27 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

    nloop_phi = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_phi_sta(inum) = nloop_phi + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../phi/gkvp.000000.0.phi."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../phi/gkvp.000000.0.phi."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_phi_end(inum) = nloop_phi + filesize/recl_phi
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_phi_sta = ",  &
                                        !                    loop_phi_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_phi_end = ",  &
                                        !                    loop_phi_end(inum)
                                        !-----------------
      !nloop_phi = loop_phi_end(inum)

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( phi_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( phi_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_phi_end(inum) = nloop_phi + t_count

      nloop_phi = loop_phi_end(inum)
    end do
    write(*,*) "# nloop_phi = ", nloop_phi

END SUBROUTINE rb_phi_setnloop


SUBROUTINE rb_phi_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 28 Feb. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_phi_sta(wknum) <= loop .and. loop <= loop_phi_end(wknum) ) then
        inum = wknum
        irec = loop - loop_phi_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_phi_loop2inumirec


SUBROUTINE rb_phi_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 27 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value


    call rb_phi_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_phi*(irec-1) + nhead
    !read( unit=100000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=phi_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)


END SUBROUTINE rb_phi_gettime


SUBROUTINE rb_phi_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 27 Feb. 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_phi_sta(snum):loop_phi_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir

    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_phi_sta(inum), loop_phi_end(inum)
        skipbyte = recl_phi*(irec-1) + nhead
        read( unit=100000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_phi_getalltime


SUBROUTINE rb_phi_mxmyizloop( mx, gmy, giz, loop, phi )
!-------------------------------------------------------------------------------
!
!     Get phi at mx, gmy, giz, loop
!                                                   (S. Maeyama, 3 March 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, loop
  complex(kind=DP), intent(out) :: phi

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: my, iz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_rephi, varid_imphi
  integer(kind=4) :: start_phi(1:4), count_phi(1:4)
  real(kind=DP), dimension(1) :: rephi
  real(kind=DP), dimension(1) :: imphi

    !call rb_gmy2rankwmy( gmy, rankw, my )
    !call rb_giz2rankziz( giz, rankz, iz )
    call rb_phi_loop2inumirec( loop, inum, irec )

    !ir = rankw + nprocw*rankz
    !skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
    !         + nhead + DP        &  ! header and time is also skipped.
    !         + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !read( unit=100000000+100000*inum+ir,  &
    !      pos=skipbyte+1 ) phi

    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "rephi", varid_rephi )
    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "imphi", varid_imphi )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_phi(:) = int((/ 1,1,1,1 /), kind=4)
    start_phi(:) = int((/ nx+mx+1,gmy+1,global_nz+giz+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_rephi, &
                              values=rephi, &
                              start=start_phi, count=count_phi )
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_imphi, &
                              values=imphi, &
                              start=start_phi, count=count_phi )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    phi = cmplx(rephi(1),imphi(1),kind=DP)

END SUBROUTINE rb_phi_mxmyizloop


SUBROUTINE rb_phi_mxmyloop( mx, gmy, loop, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(giz) at mx, gmy, loop
!                                                   (S. Maeyama, 3 March 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1) :: phi

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: my, iz, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_rephi, varid_imphi
  integer(kind=4) :: start_phi(1:4), count_phi(1:4)
  real(kind=DP), dimension(-global_nz:global_nz-1) :: rephi
  real(kind=DP), dimension(-global_nz:global_nz-1) :: imphi

    !call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_phi_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  ir = rankw + nprocw*rankz
    !  do iz = -nz, nz-1
    !    giz = - global_nz + 2*nz * rankz + iz + nz
    !    skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
    !             + nhead + DP        &  ! header and time is also skipped.
    !             + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !    read( unit=100000000+100000*inum+ir,  &
    !          pos=skipbyte+1 ) phi(giz)
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "rephi", varid_rephi )
    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "imphi", varid_imphi )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_phi(:) = int((/ 1,1,2*global_nz,1 /), kind=4)
    start_phi(:) = int((/ nx+mx+1,gmy+1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_rephi, &
                              values=rephi(:), &
                              start=start_phi, count=count_phi )
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_imphi, &
                              values=imphi(:), &
                              start=start_phi, count=count_phi )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    phi = cmplx(rephi,imphi,kind=DP)

END SUBROUTINE rb_phi_mxmyloop


SUBROUTINE rb_phi_izloop( giz, loop, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(mx,gmy) at giz, loop
!                                                   (S. Maeyama, 3 March 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz, loop
  complex(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny) :: phi

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: iz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_rephi, varid_imphi
  integer(kind=4) :: start_phi(1:4), count_phi(1:4)
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: rephi
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: imphi


    !call rb_giz2rankziz( giz, rankz, iz )
    call rb_phi_loop2inumirec( loop, inum, irec )

    !do rankw = 0, nprocw-1
    !  ir = rankw + nprocw*rankz
    !  skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
    !           + nhead + DP        &  ! header and time is also skipped.
    !           + ( (2*nx+1)*(ny+1)*(iz+nz) )*(2*DP)
    !  read( unit=100000000+100000*inum+ir, pos=skipbyte+1 )  &
    !      phi(-nx:nx,ist_y_g(rankw):iend_y_g(rankw))
    !end do

    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "rephi", varid_rephi )
    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "imphi", varid_imphi )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_phi(:) = int((/ 2*nx+1,global_ny+1,1,1 /), kind=4)
    start_phi(:) = int((/ 1,1,global_nz+giz+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_rephi, &
                              values=rephi(:,:), &
                              start=start_phi, count=count_phi )
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_imphi, &
                              values=imphi(:,:), &
                              start=start_phi, count=count_phi )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    phi = cmplx(rephi,imphi,kind=DP)

END SUBROUTINE rb_phi_izloop


SUBROUTINE rb_phi_loop( loop, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(mx,gmy,giz) at loop
!                                                   (S. Maeyama, 28 Feb. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wkphi
  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, gmy, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_rephi, varid_imphi
  integer(kind=4) :: start_phi(1:4), count_phi(1:4)
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: rephi
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: imphi


    call rb_phi_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*rankz
    !    skipbyte = recl_phi*(irec-1) + nhead + DP
    !    read( unit=100000000+100000*inum+ir, pos=skipbyte+1 ) wkphi
    !    do iz = -nz, nz-1
    !      giz = - global_nz + 2*nz * rankz + iz + nz
    !      do my = 0, ny
    !        gmy = ( ny+1 ) * rankw + my
    !        if ( gmy <= global_ny ) then
    !          do mx = -nx, nx
    !            phi(mx,gmy,giz) = wkphi(mx,my,iz)
    !          end do
    !        end if
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "rephi", varid_rephi )
    ierr_nf90 = nf90_inq_varid( phi_nc(inum), "imphi", varid_imphi )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_phi(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,1 /), kind=4)
    start_phi(:) = int((/ 1,1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_rephi, &
                              values=rephi(:,:,:), &
                              start=start_phi, count=count_phi )
    ierr_nf90 = nf90_get_var( phi_nc(inum), varid_imphi, &
                              values=imphi(:,:,:), &
                              start=start_phi, count=count_phi )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    phi = cmplx(rephi,imphi,kind=DP)


END SUBROUTINE rb_phi_loop


SUBROUTINE rb_phi_mxmyiz( mx, gmy, giz, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(loop) at mx, gmy, giz
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz
  complex(kind=DP), intent(out), &
    dimension(loop_phi_sta(snum):loop_phi_end(enum)) :: phi

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )

    ir = rankw + nprocw*rankz
    do inum = snum, enum
      irec = 1
      do loop = loop_phi_sta(inum), loop_phi_end(inum)
        skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
        read( unit=100000000+100000*inum+ir,  &
              pos=skipbyte+1 ) phi(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_phi_mxmyiz


SUBROUTINE rb_phi_mxmy( mx, gmy, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(giz,loop) at mx, gmy
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy
  complex(kind=DP), intent(out), &
    dimension(-global_nz:global_nz-1,loop_phi_sta(snum):loop_phi_end(enum)) :: phi

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )

    do inum = snum, enum
      irec = 1
      do loop = loop_phi_sta(inum), loop_phi_end(inum)
        do rankz = 0, nprocz-1
          ir = rankw + nprocw*rankz
          do iz = -nz, nz-1
            giz = - global_nz + 2*nz * rankz + iz + nz
            skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
                     + nhead + DP        &  ! header and time is also skipped.
                     + ( (2*nx+1)*(ny+1)*(iz+nz)+(2*nx+1)*my+(mx+nx) )*(2*DP)
            read( unit=100000000+100000*inum+ir,  &
                  pos=skipbyte+1 ) phi(giz,loop)
          end do
        end do
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_phi_mxmy


SUBROUTINE rb_phi_myloop( gmy, loop, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(mx,giz) at gmy, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: gmy, loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,-global_nz:global_nz-1) :: phi

  integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_phi_loop2inumirec( loop, inum, irec )

    do rankz = 0, nprocz-1
      ir = rankw + nprocw*rankz
      do iz = -nz, nz-1
        giz = - global_nz + 2*nz * rankz + iz + nz
        skipbyte = recl_phi*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(iz+nz)+(2*nx+1)*my )*(2*DP)
        read( unit=100000000+100000*inum+ir,  &
              pos=skipbyte+1 ) phi(-nx:nx,giz)
      end do
    end do

END SUBROUTINE rb_phi_myloop


SUBROUTINE rb_Al_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *Al*
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, rankz, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  character(len=3) :: cnum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do rankz = 0, nprocz-1
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*rankz
      !    write( crank, fmt="(i6.6)" ) ir
      !    open( unit=200000000+100000*inum+ir,                      &
      !          file="../phi/gkvp."//crank//".0.Al."//cnum,  &
      !          status="old", action="read",                        &
      !          form="unformatted", access="stream" )
                                        !--- for debug ---
                                        !write(*,*) 200000000+100000*inum+ir, &
                                        !           "../phi/gkvp."      &
                                        !           //crank//".0.Al."//cnum
                                        !-----------------
      !  end do
      !end do
      ierr_nf90 = nf90_open( path="../phi/gkvp.Al."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=Al_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_Al_fileopen


SUBROUTINE rb_Al_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *Al*
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, rankz, rankw, ir
  integer :: inum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do rankz = 0, nprocz-1
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*rankz
      !    close( unit=200000000+100000*inum+ir )
      !  end do
      !end do

      ierr_nf90 = nf90_close( Al_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_Al_fileclose


SUBROUTINE rb_Al_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz
  integer(kind=4) :: len_kx, len_ky, len_zz

    do inum = snum, enum
      !ir = 0
      !read( unit=200000000+100000*inum+ir, pos=1 ) head
      !read( unit=200000000+100000*inum+ir, pos=recl_Al-4+1 ) foot

      !if ( head /= foot .or. head /= recl_Al - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *Al*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_Al-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if

      ierr_nf90 = nf90_inq_varid( Al_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( Al_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( Al_nc(inum), "zz", varid_zz )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( Al_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( Al_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( Al_nc(inum), varid_zz, len=len_zz )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           len_zz /= 2*global_nz ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.Al*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#              zz = ", len_zz
        write(*,*) "#     2*global_nz = ", 2*global_nz
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *Al* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_Al-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
    write(*,*) "# The netCDF files gkvp.Al*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#              zz = ", len_zz
    write(*,*) "#     2*global_nz = ", 2*global_nz
                                        !--- for debug ---
                                        !write(*,*) head, recl_Al-2*nhead, foot
                                        !-----------------

END SUBROUTINE rb_Al_check


SUBROUTINE rb_Al_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

    nloop_Al = -1
    do inum = 1, enum
      write( cnum, fmt="(i3.3)" ) inum

      loop_Al_sta(inum) = nloop_Al + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../phi/gkvp.000000.0.Al."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../phi/gkvp.000000.0.Al."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_Al_end(inum) = nloop_Al + filesize/recl_Al
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_Al_sta = ",  &
                                        !                    loop_Al_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_Al_end = ",  &
                                        !                    loop_Al_end(inum)
                                        !-----------------
      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( Al_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( Al_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_Al_end(inum) = nloop_Al + t_count
      nloop_Al = loop_Al_end(inum)
    end do
    write(*,*) "# nloop_Al  = ", nloop_Al

END SUBROUTINE rb_Al_setnloop


SUBROUTINE rb_Al_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_Al_sta(wknum) <= loop .and. loop <= loop_Al_end(wknum) ) then
        inum = wknum
        irec = loop - loop_Al_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_Al_loop2inumirec


SUBROUTINE rb_Al_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value


    call rb_Al_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_Al*(irec-1) + nhead
    !read( unit=200000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=Al_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)



END SUBROUTINE rb_Al_gettime


SUBROUTINE rb_Al_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_Al_sta(snum):loop_Al_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir


    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_Al_sta(inum), loop_Al_end(inum)
        skipbyte = recl_Al*(irec-1) + nhead
        read( unit=200000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_Al_getalltime


SUBROUTINE rb_Al_mxmyizloop( mx, gmy, giz, loop, Al )
!-------------------------------------------------------------------------------
!
!     Get Al at mx, gmy, giz, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, loop
  complex(kind=DP), intent(out) :: Al

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: my, iz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_reAl, varid_imAl
  integer(kind=4) :: start_Al(1:4), count_Al(1:4)
  real(kind=DP), dimension(1) :: reAl
  real(kind=DP), dimension(1) :: imAl

    !call rb_gmy2rankwmy( gmy, rankw, my )
    !call rb_giz2rankziz( giz, rankz, iz )
    call rb_Al_loop2inumirec( loop, inum, irec )

    !ir = rankw + nprocw*rankz
    !skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
    !         + nhead + DP        &  ! header and time is also skipped.
    !         + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !read( unit=200000000+100000*inum+ir,  &
    !      pos=skipbyte+1 ) Al

    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "reAl", varid_reAl )
    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "imAl", varid_imAl )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_Al(:) = int((/ 1,1,1,1 /), kind=4)
    start_Al(:) = int((/ nx+mx+1,gmy+1,global_nz+giz+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_reAl, &
                              values=reAl, &
                              start=start_Al, count=count_Al )
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_imAl, &
                              values=imAl, &
                              start=start_Al, count=count_Al )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    Al = cmplx(reAl(1),imAl(1),kind=DP)

END SUBROUTINE rb_Al_mxmyizloop


SUBROUTINE rb_Al_mxmyloop( mx, gmy, loop, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(giz) at mx, gmy, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1) :: Al

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: my, iz, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_reAl, varid_imAl
  integer(kind=4) :: start_Al(1:4), count_Al(1:4)
  real(kind=DP), dimension(-global_nz:global_nz-1) :: reAl
  real(kind=DP), dimension(-global_nz:global_nz-1) :: imAl

    !call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_Al_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  ir = rankw + nprocw*rankz
    !  do iz = -nz, nz-1
    !    giz = - global_nz + 2*nz * rankz + iz + nz
    !    skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
    !             + nhead + DP        &  ! header and time is also skipped.
    !             + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !    read( unit=200000000+100000*inum+ir,  &
    !          pos=skipbyte+1 ) Al(giz)
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "reAl", varid_reAl )
    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "imAl", varid_imAl )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_Al(:) = int((/ 1,1,2*global_nz,1 /), kind=4)
    start_Al(:) = int((/ nx+mx+1,gmy+1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_reAl, &
                              values=reAl(:), &
                              start=start_Al, count=count_Al )
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_imAl, &
                              values=imAl(:), &
                              start=start_Al, count=count_Al )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    Al = cmplx(reAl,imAl,kind=DP)

END SUBROUTINE rb_Al_mxmyloop


SUBROUTINE rb_Al_izloop( giz, loop, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(mx,gmy) at giz, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz, loop
  complex(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny) :: Al

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: iz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_reAl, varid_imAl
  integer(kind=4) :: start_Al(1:4), count_Al(1:4)
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: reAl
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: imAl

    !call rb_giz2rankziz( giz, rankz, iz )
    call rb_Al_loop2inumirec( loop, inum, irec )

    !do rankw = 0, nprocw-1
    !  ir = rankw + nprocw*rankz
    !  skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
    !           + nhead + DP        &  ! header and time is also skipped.
    !           + ( (2*nx+1)*(ny+1)*(iz+nz) )*(2*DP)
    !  read( unit=200000000+100000*inum+ir, pos=skipbyte+1 )  &
    !      Al(-nx:nx,ist_y_g(rankw):iend_y_g(rankw))
    !end do

    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "reAl", varid_reAl )
    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "imAl", varid_imAl )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_Al(:) = int((/ 2*nx+1,global_ny+1,1,1 /), kind=4)
    start_Al(:) = int((/ 1,1,global_nz+giz+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_reAl, &
                              values=reAl(:,:), &
                              start=start_Al, count=count_Al )
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_imAl, &
                              values=imAl(:,:), &
                              start=start_Al, count=count_Al )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    Al = cmplx(reAl,imAl,kind=DP)

END SUBROUTINE rb_Al_izloop


SUBROUTINE rb_Al_loop( loop, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(mx,gmy,giz) at loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wkAl
  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, gmy, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_reAl, varid_imAl
  integer(kind=4) :: start_Al(1:4), count_Al(1:4)
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: reAl
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: imAl


    call rb_Al_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*rankz
    !    skipbyte = recl_Al*(irec-1) + nhead + DP
    !    read( unit=200000000+100000*inum+ir, pos=skipbyte+1 ) wkAl
    !    do iz = -nz, nz-1
    !      giz = - global_nz + 2*nz * rankz + iz + nz
    !      do my = 0, ny
    !        gmy = ( ny+1 ) * rankw + my
    !        if ( gmy <= global_ny ) then
    !          do mx = -nx, nx
    !            Al(mx,gmy,giz) = wkAl(mx,my,iz)
    !          end do
    !        end if
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "reAl", varid_reAl )
    ierr_nf90 = nf90_inq_varid( Al_nc(inum), "imAl", varid_imAl )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_Al(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,1 /), kind=4)
    start_Al(:) = int((/ 1,1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_reAl, &
                              values=reAl(:,:,:), &
                              start=start_Al, count=count_Al )
    ierr_nf90 = nf90_get_var( Al_nc(inum), varid_imAl, &
                              values=imAl(:,:,:), &
                              start=start_Al, count=count_Al )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    Al = cmplx(reAl,imAl,kind=DP)

END SUBROUTINE rb_Al_loop


SUBROUTINE rb_Al_mxmyiz( mx, gmy, giz, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(loop) at mx, gmy, giz
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz
  complex(kind=DP), intent(out), &
    dimension(loop_Al_sta(snum):loop_Al_end(enum)) :: Al

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )

    ir = rankw + nprocw*rankz
    do inum = snum, enum
      irec = 1
      do loop = loop_Al_sta(inum), loop_Al_end(inum)
        skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
        read( unit=200000000+100000*inum+ir,  &
              pos=skipbyte+1 ) Al(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_Al_mxmyiz


SUBROUTINE rb_Al_mxmy( mx, gmy, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(giz,loop) at mx, gmy
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy
  complex(kind=DP), intent(out), &
    dimension(-global_nz:global_nz-1,loop_Al_sta(snum):loop_Al_end(enum)) :: Al

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )

    do inum = snum, enum
      irec = 1
      do loop = loop_Al_sta(inum), loop_Al_end(inum)
        do rankz = 0, nprocz-1
          ir = rankw + nprocw*rankz
          do iz = -nz, nz-1
            giz = - global_nz + 2*nz * rankz + iz + nz
            skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
                     + nhead + DP        &  ! header and time is also skipped.
                     + ( (2*nx+1)*(ny+1)*(iz+nz)+(2*nx+1)*my+(mx+nx) )*(2*DP)
            read( unit=200000000+100000*inum+ir,  &
                  pos=skipbyte+1 ) Al(giz,loop)
          end do
        end do
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_Al_mxmy


SUBROUTINE rb_Al_myloop( gmy, loop, Al )
!-------------------------------------------------------------------------------
!
!     Get Al(mx,giz) at gmy, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: gmy, loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,-global_nz:global_nz-1) :: Al

  integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_Al_loop2inumirec( loop, inum, irec )

    do rankz = 0, nprocz-1
      ir = rankw + nprocw*rankz
      do iz = -nz, nz-1
        giz = - global_nz + 2*nz * rankz + iz + nz
        skipbyte = recl_Al*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(iz+nz)+(2*nx+1)*my )*(2*DP)
        read( unit=200000000+100000*inum+ir,  &
              pos=skipbyte+1 ) Al(-nx:nx,giz)
      end do
    end do

END SUBROUTINE rb_Al_myloop


SUBROUTINE rb_mom_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *mom*
!                                                   (S. Maeyama, 16 March 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankz, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  !character(len=1) :: srank
  character(len=3) :: cnum

  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do ranks = 0, nprocs-1
      !  write( srank, fmt="(i1.1)" ) ranks
      !  do rankz = 0, nprocz-1
      !    do rankw = 0, nprocw-1
      !      ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*ranks
      !      write( crank, fmt="(i6.6)" ) ir
      !      open( unit=400000000+100000*inum+ir,                      &
      !            file="../phi/gkvp."//crank//"."//srank//".mom."//cnum, &
      !            status="old", action="read",                        &
      !            form="unformatted", access="stream" )
                                        !--- for debug ---
                                        !write(*,*) 400000000+100000*inum+ir, &
                                        !           "../phi/gkvp."      &
                                        !   //crank//"."//srank//".mom."//cnum
                                        !-----------------
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_open( path="../phi/gkvp.mom."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=mom_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_mom_fileopen


SUBROUTINE rb_mom_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *mom*
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankz, rankw, ir
  integer :: inum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do ranks = 0, nprocs-1
      !  do rankz = 0, nprocz-1
      !    do rankw = 0, nprocw-1
      !      ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*ranks
      !      close( unit=400000000+100000*inum+ir )
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_close( mom_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_mom_fileclose


SUBROUTINE rb_mom_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_imom, varid_is
  integer(kind=4) :: len_kx, len_ky, len_zz, len_imom, len_is

    do inum = snum, enum
      !ir = 0
      !read( unit=400000000+100000*inum+ir, pos=1 ) head
      !read( unit=400000000+100000*inum+ir, pos=recl_mom-4+1 ) foot

      !if ( head /= foot .or. head /= recl_mom - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *mom*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_mom-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if
      ierr_nf90 = nf90_inq_varid( mom_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( mom_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( mom_nc(inum), "zz", varid_zz )
      ierr_nf90 = nf90_inq_varid( mom_nc(inum), "imom", varid_imom )
      ierr_nf90 = nf90_inq_varid( mom_nc(inum), "is", varid_is )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), varid_zz, len=len_zz )
      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), varid_imom, len=len_imom )
      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), varid_is, len=len_is )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           len_zz /= 2*global_nz .or. len_imom /= nmom .or. &
           len_is /= nprocs ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.mom*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#              zz = ", len_zz
        write(*,*) "#     2*global_nz = ", 2*global_nz
        write(*,*) "#            imom = ", len_imom
        write(*,*) "#            nmom = ", nmom
        write(*,*) "#              is = ", len_is
        write(*,*) "#          nprocs = ", nprocs
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *mom* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_mom-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
                                        !--- for debug ---
                                        !write(*,*) head, recl_mom-2*nhead, foot
                                        !-----------------
    write(*,*) "# The netCDF files gkvp.mom*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#              zz = ", len_zz
    write(*,*) "#     2*global_nz = ", 2*global_nz
    write(*,*) "#            imom = ", len_imom
    write(*,*) "#            nmom = ", nmom
    write(*,*) "#              is = ", len_is
    write(*,*) "#          nprocs = ", nprocs

END SUBROUTINE rb_mom_check


SUBROUTINE rb_mom_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

    nloop_mom = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_mom_sta(inum) = nloop_mom + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../phi/gkvp.000000.0.mom."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../phi/gkvp.000000.0.mom."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_mom_end(inum) = nloop_mom + filesize/recl_mom
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_mom_sta = ",  &
                                        !                    loop_mom_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_mom_end = ",  &
                                        !                    loop_mom_end(inum)
                                        !-----------------
      !nloop_mom = loop_mom_end(inum)

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( mom_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( mom_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_mom_end(inum) = nloop_mom + t_count

      nloop_mom = loop_mom_end(inum)
    end do
    write(*,*) "# nloop_mom = ", nloop_mom

END SUBROUTINE rb_mom_setnloop


SUBROUTINE rb_mom_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_mom_sta(wknum) <= loop .and. loop <= loop_mom_end(wknum) ) then
        inum = wknum
        irec = loop - loop_mom_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_mom_loop2inumirec


SUBROUTINE rb_mom_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value

    call rb_mom_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_mom*(irec-1) + nhead
    !read( unit=400000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=mom_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)

END SUBROUTINE rb_mom_gettime


SUBROUTINE rb_mom_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_mom_sta(snum):loop_mom_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir

    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_mom_sta(inum), loop_mom_end(inum)
        skipbyte = recl_mom*(irec-1) + nhead
        read( unit=400000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_mom_getalltime


SUBROUTINE rb_mom_mxmyizimomisloop( mx, gmy, giz, imom, is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom at mx, gmy, giz, imom, is, loop
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, imom, is, loop
  complex(kind=DP), intent(out) :: mom

  integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )
    call rb_mom_loop2inumirec( loop, inum, irec )

    ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
             + nhead + DP        &  ! header and time is also skipped.
             + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
               + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    read( unit=400000000+100000*inum+ir,  &
          pos=skipbyte+1 ) mom

END SUBROUTINE rb_mom_mxmyizimomisloop


SUBROUTINE rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(giz) at mx, gmy, imom, is, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, imom, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1) :: mom

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: my, iz, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_remom, varid_immom
  integer(kind=4) :: start_mom(1:6), count_mom(1:6)
  real(kind=DP), dimension(-global_nz:global_nz-1) :: remom
  real(kind=DP), dimension(-global_nz:global_nz-1) :: immom

    !call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_mom_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !  do iz = -nz, nz-1
    !    giz = - global_nz + 2*nz * rankz + iz + nz
    !    skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
    !             + nhead + DP        &  ! header and time is also skipped.
    !             + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
    !               + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !    read( unit=400000000+100000*inum+ir,  &
    !          pos=skipbyte+1 ) mom(giz)
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "remom", varid_remom )
    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "immom", varid_immom )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_mom(:) = int((/ 1,1,2*global_nz,1,1,1 /), kind=4)
    start_mom(:) = int((/ nx+mx+1,gmy+1,1,imom+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_remom, &
                              values=remom(:), &
                              start=start_mom, count=count_mom )
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_immom, &
                              values=immom(:), &
                              start=start_mom, count=count_mom )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    mom = cmplx(remom,immom,kind=DP)

END SUBROUTINE rb_mom_mxmyimomisloop


SUBROUTINE rb_mom_izimomisloop( giz, imom, is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(mx,gmy) at giz, imom, is, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz, imom, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny) :: mom

  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: iz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_remom, varid_immom
  integer(kind=4) :: start_mom(1:6), count_mom(1:6)
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: remom
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: immom

    !call rb_giz2rankziz( giz, rankz, iz )
    call rb_mom_loop2inumirec( loop, inum, irec )

    !do rankw = 0, nprocw-1
    !  ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !  skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
    !           + nhead + DP        &  ! header and time is also skipped.
    !           + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
    !             + (2*nx+1)*(ny+1)*(iz+nz) )*(2*DP)
    !  read( unit=400000000+100000*inum+ir, pos=skipbyte+1 )  &
    !      mom(-nx:nx,ist_y_g(rankw):iend_y_g(rankw))
    !end do

    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "remom", varid_remom )
    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "immom", varid_immom )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_mom(:) = int((/ 2*nx+1,global_ny+1,1,1,1,1 /), kind=4)
    start_mom(:) = int((/ 1,1,global_nz+giz+1,imom+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_remom, &
                              values=remom(:,:), &
                              start=start_mom, count=count_mom )
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_immom, &
                              values=immom(:,:), &
                              start=start_mom, count=count_mom )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    mom = cmplx(remom,immom,kind=DP)

END SUBROUTINE rb_mom_izimomisloop


SUBROUTINE rb_mom_imomisloop( imom, is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(mx,gmy,giz) at imom, is, loop
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: imom, is, loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: mom

  complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wkmom
  integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: mx, my, iz, gmy, giz

    call rb_mom_loop2inumirec( loop, inum, irec )

    do rankz = 0, nprocz-1
      do rankw = 0, nprocw-1
        ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
        skipbyte = recl_mom*(irec-1) + nhead + DP  &
                 + ( (2*nx+1)*(ny+1)*(2*nz)*imom )*(2*DP)
        read( unit=400000000+100000*inum+ir, pos=skipbyte+1 ) wkmom
        do iz = -nz, nz-1
          giz = - global_nz + 2*nz * rankz + iz + nz
          do my = 0, ny
            gmy = ( ny+1 ) * rankw + my
            if ( gmy <= global_ny ) then
              do mx = -nx, nx
                mom(mx,gmy,giz) = wkmom(mx,my,iz)
              end do
            end if
          end do
        end do
      end do
    end do

END SUBROUTINE rb_mom_imomisloop


SUBROUTINE rb_mom_isloop( is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(mx,gmy,giz,imom) at is, loop
!                                                   (S. Maeyama, 10 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: is, loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: mom

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,0:nmom-1) :: wkmom
  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, gmy, giz, imom

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_remom, varid_immom
  integer(kind=4) :: start_mom(1:6), count_mom(1:6)
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: remom
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: immom

    call rb_mom_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    skipbyte = recl_mom*(irec-1) + nhead + DP
    !    read( unit=400000000+100000*inum+ir, pos=skipbyte+1 ) wkmom
    !    do imom = 0, nmom-1
    !      do iz = -nz, nz-1
    !        giz = - global_nz + 2*nz * rankz + iz + nz
    !        do my = 0, ny
    !          gmy = ( ny+1 ) * rankw + my
    !          if ( gmy <= global_ny ) then
    !            do mx = -nx, nx
    !              mom(mx,gmy,giz,imom) = wkmom(mx,my,iz,imom)
    !            end do
    !          end if
    !        end do
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "remom", varid_remom )
    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "immom", varid_immom )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_mom(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,nmom,1,1 /), kind=4)
    start_mom(:) = int((/ 1,1,1,1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_remom, &
                              values=remom(:,:,:,:), &
                              start=start_mom, count=count_mom )
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_immom, &
                              values=immom(:,:,:,:), &
                              start=start_mom, count=count_mom )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    mom = cmplx(remom,immom,kind=DP)

END SUBROUTINE rb_mom_isloop


SUBROUTINE rb_mom_loop( loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(mx,gmy,giz,imom,is) at loop
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1,0:nmom-1) :: wkmom
  !integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, gmy, giz, imom, is

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_remom, varid_immom
  integer(kind=4) :: start_mom(1:6), count_mom(1:6)
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: remom
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: immom

    call rb_mom_loop2inumirec( loop, inum, irec )

    !do is = 0, ns-1
    !  do rankz = 0, nprocz-1
    !    do rankw = 0, nprocw-1
    !      ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !      skipbyte = recl_mom*(irec-1) + nhead + DP
    !      read( unit=400000000+100000*inum+ir, pos=skipbyte+1 ) wkmom
    !      do imom = 0, nmom-1
    !        do iz = -nz, nz-1
    !          giz = - global_nz + 2*nz * rankz + iz + nz
    !          do my = 0, ny
    !            gmy = ( ny+1 ) * rankw + my
    !            if ( gmy <= global_ny ) then
    !              do mx = -nx, nx
    !                mom(mx,gmy,giz,imom,is) = wkmom(mx,my,iz,imom)
    !              end do
    !            end if
    !          end do
    !        end do
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "remom", varid_remom )
    ierr_nf90 = nf90_inq_varid( mom_nc(inum), "immom", varid_immom )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_mom(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,nmom,ns,1 /), kind=4)
    start_mom(:) = int((/ 1,1,1,1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_remom, &
                              values=remom(:,:,:,:,:), &
                              start=start_mom, count=count_mom )
    ierr_nf90 = nf90_get_var( mom_nc(inum), varid_immom, &
                              values=immom(:,:,:,:,:), &
                              start=start_mom, count=count_mom )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    mom = cmplx(remom,immom,kind=DP)

END SUBROUTINE rb_mom_loop


SUBROUTINE rb_mom_mxmyizimomis( mx, gmy, giz, imom, is, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(loop) at mx, gmy, giz, imom, is
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, imom, is
  complex(kind=DP), intent(out), &
    dimension(loop_mom_sta(snum):loop_mom_end(enum)) :: mom

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )

    ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    do inum = snum, enum
      irec = 1
      do loop = loop_mom_sta(inum), loop_mom_end(inum)
        skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
                   + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
        read( unit=400000000+100000*inum+ir,  &
              pos=skipbyte+1 ) mom(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_mom_mxmyizimomis


SUBROUTINE rb_mom_mxmyimomis( mx, gmy, imom, is, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(giz,loop) at mx, gmy
!                                                   (S. Maeyama, 10 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, imom, is
  complex(kind=DP), intent(out), &
    dimension(-global_nz:global_nz-1,loop_mom_sta(snum):loop_mom_end(enum)) :: mom

  integer :: loop, inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )

    do inum = snum, enum
      irec = 1
      do loop = loop_mom_sta(inum), loop_mom_end(inum)
        do rankz = 0, nprocz-1
          ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
          do iz = -nz, nz-1
            giz = - global_nz + 2*nz * rankz + iz + nz
            skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
                     + nhead + DP        &  ! header and time is also skipped.
                     + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
                       + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
            read( unit=400000000+100000*inum+ir,  &
                  pos=skipbyte+1 ) mom(giz,loop)
          end do
        end do
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_mom_mxmyimomis


SUBROUTINE rb_mom_myimomisloop( gmy, imom, is, loop, mom )
!-------------------------------------------------------------------------------
!
!     Get mom(mx,giz) at gmy, imom, is, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: gmy, imom, is, loop
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,-global_nz:global_nz-1) :: mom

  integer :: inum, irec, skipbyte, ir, rankz, rankw
  integer :: my, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_mom_loop2inumirec( loop, inum, irec )

    do rankz = 0, nprocz-1
      ir = rankw + nprocw*rankz + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
      do iz = -nz, nz-1
        giz = - global_nz + 2*nz * rankz + iz + nz
        skipbyte = recl_mom*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(2*nz)*imom  &
                   + (2*nx+1)*(ny+1)*(iz+nz)+(2*nx+1)*my )*(2*DP)
        read( unit=400000000+100000*inum+ir,  &
              pos=skipbyte+1 ) mom(-nx:nx,giz)
      end do
    end do

END SUBROUTINE rb_mom_myimomisloop


SUBROUTINE rb_trn_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *trn*
!                                                   (S. Maeyama, 24 Aug 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  !character(len=1) :: srank
  character(len=3) :: cnum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do ranks = 0, nprocs-1
      !  write( srank, fmt="(i1.1)" ) ranks
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*nprocz*nprocv*nprocm*ranks
      !    write( crank, fmt="(i6.6)" ) ir
      !    open( unit=500000000+100000*inum+ir,                      &
      !          file="../phi/gkvp."//crank//"."//srank//".trn."//cnum, &
      !          status="old", action="read",                        &
      !          form="unformatted", access="stream" )
                                      !--- for debug ---
                                      !write(*,*) 500000000+100000*inum+ir, &
                                      !           "../phi/gkvp."      &
                                      !   //crank//"."//srank//".trn."//cnum
                                      !-----------------
      !  end do
      !end do

      ierr_nf90 = nf90_open( path="../phi/gkvp.trn."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=trn_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_trn_fileopen


SUBROUTINE rb_trn_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *trn*
!                                                   (S. Maeyama, 24 Aug 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankw, ir
  integer :: inum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do ranks = 0, nprocs-1
      !  do rankw = 0, nprocw-1
      !    ir = rankw + nprocw*nprocz*nprocv*nprocm*ranks
      !    close( unit=500000000+100000*inum+ir )
      !  end do
      !end do

      ierr_nf90 = nf90_close( trn_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_trn_fileclose


SUBROUTINE rb_trn_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 24 Aug 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_itrn, varid_is
  integer(kind=4) :: len_kx, len_ky, len_itrn, len_is

    do inum = snum, enum
      !ir = 0
      !read( unit=500000000+100000*inum+ir, pos=1 ) head
      !read( unit=500000000+100000*inum+ir, pos=recl_trn-4+1 ) foot

      !if ( head /= foot .or. head /= recl_trn - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *trn*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_trn-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if

      ierr_nf90 = nf90_inq_varid( trn_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( trn_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( trn_nc(inum), "itrn", varid_itrn )
      ierr_nf90 = nf90_inq_varid( trn_nc(inum), "is", varid_is )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( trn_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( trn_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( trn_nc(inum), varid_itrn, len=len_itrn )
      ierr_nf90 = nf90_inquire_dimension( trn_nc(inum), varid_is, len=len_is )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           len_itrn /= ntrn .or. len_is /= nprocs ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.trn*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#            itrn = ", len_itrn
        write(*,*) "#            ntrn = ", ntrn
        write(*,*) "#              is = ", len_is
        write(*,*) "#          nprocs = ", nprocs
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *trn* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_trn-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
                                        !--- for debug ---
                                        !write(*,*) head, recl_trn-2*nhead, foot
                                        !-----------------
    write(*,*) "# The netCDF files gkvp.trn*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#            itrn = ", len_itrn
    write(*,*) "#            ntrn = ", ntrn
    write(*,*) "#              is = ", len_is
    write(*,*) "#          nprocs = ", nprocs

END SUBROUTINE rb_trn_check


SUBROUTINE rb_trn_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 24 Aug 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

    nloop_trn = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_trn_sta(inum) = nloop_trn + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../phi/gkvp.000000.0.trn."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../phi/gkvp.000000.0.trn."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_trn_end(inum) = nloop_trn + filesize/recl_trn
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_trn_sta = ",  &
                                        !                    loop_trn_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_trn_end = ",  &
                                        !                    loop_trn_end(inum)
                                        !-----------------
      !nloop_trn = loop_trn_end(inum)

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( trn_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( trn_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_trn_end(inum) = nloop_trn + t_count

      nloop_trn = loop_trn_end(inum)
    end do
    write(*,*) "# nloop_trn = ", nloop_trn

END SUBROUTINE rb_trn_setnloop


SUBROUTINE rb_trn_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_trn_sta(wknum) <= loop .and. loop <= loop_trn_end(wknum) ) then
        inum = wknum
        irec = loop - loop_trn_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_trn_loop2inumirec


SUBROUTINE rb_trn_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value

    call rb_trn_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_trn*(irec-1) + nhead
    !read( unit=500000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( trn_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=trn_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)

END SUBROUTINE rb_trn_gettime


SUBROUTINE rb_trn_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_trn_sta(snum):loop_trn_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir

    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_trn_sta(inum), loop_trn_end(inum)
        skipbyte = recl_trn*(irec-1) + nhead
        read( unit=500000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_trn_getalltime


SUBROUTINE rb_trn_mxmyitrnisloop( mx, gmy, itrn, is, loop, trn )
!-------------------------------------------------------------------------------
!
!     Get trn at mx, gmy, itrn, is, loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, itrn, is, loop
  real(kind=DP), intent(out) :: trn

  integer :: inum, irec, skipbyte, ir, rankw
  integer :: my

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_trn_loop2inumirec( loop, inum, irec )

    ir = rankw + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    skipbyte = recl_trn*(irec-1) &  ! (irec-1) lines are skipped.
             + nhead + DP        &  ! header and time is also skipped.
             + ( (2*nx+1)*(ny+1)*itrn + (2*nx+1)*my + (mx+nx) )*DP
    read( unit=500000000+100000*inum+ir,  &
          pos=skipbyte+1 ) trn

END SUBROUTINE rb_trn_mxmyitrnisloop


SUBROUTINE rb_trn_itrnisloop( itrn, is, loop, trn )
!-------------------------------------------------------------------------------
!
!     Get trn(mx,gmy) at itrn, is, loop
!                                                   (S. Maeyama, 24 June 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: itrn, is, loop
  real(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny) :: trn

  integer :: inum, irec, skipbyte, ir, rankw

    call rb_trn_loop2inumirec( loop, inum, irec )

    do rankw = 0, nprocw-1
      ir = rankw + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
      skipbyte = recl_trn*(irec-1) &  ! (irec-1) lines are skipped.
               + nhead + DP        &  ! header and time is also skipped.
               + ( (2*nx+1)*(ny+1)*itrn )*DP
      read( unit=500000000+100000*inum+ir, pos=skipbyte+1 )  &
          trn(-nx:nx,ist_y_g(rankw):iend_y_g(rankw))
    end do

END SUBROUTINE rb_trn_itrnisloop


SUBROUTINE rb_trn_isloop( is, loop, trn )
!-------------------------------------------------------------------------------
!
!     Get trn(mx,gmy,itrn) at is, loop
!                                                   (S. Maeyama, 24 June 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: is, loop
  real(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,0:ntrn-1) :: trn

  !real(kind=DP), dimension(-nx:nx,0:ny,0:ntrn-1) :: wktrn
  !integer :: inum, irec, skipbyte, ir, rankw
  integer :: inum, irec
  !integer :: mx, my, gmy, itrn

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_trn
  integer(kind=4) :: start_trn(1:5), count_trn(1:5)

    call rb_trn_loop2inumirec( loop, inum, irec )

    !do rankw = 0, nprocw-1
    !  ir = rankw + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !  skipbyte = recl_trn*(irec-1) + nhead + DP
    !  read( unit=500000000+100000*inum+ir, pos=skipbyte+1 ) wktrn
    !  do itrn = 0, ntrn-1
    !    do my = 0, ny
    !      gmy = ( ny+1 ) * rankw + my
    !      if ( gmy <= global_ny ) then
    !        do mx = -nx, nx
    !          trn(mx,gmy,itrn) = wktrn(mx,my,itrn)
    !        end do
    !      end if
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( trn_nc(inum), "trn", varid_trn )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_trn(:) = int((/ 2*nx+1,global_ny+1,ntrn,1,1 /), kind=4)
    start_trn(:) = int((/ 1,1,1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( trn_nc(inum), varid_trn, &
                              values=trn(:,:,:), &
                              start=start_trn, count=count_trn )
    call check_nf90err(ierr_nf90, "nf90_get_var")

END SUBROUTINE rb_trn_isloop


SUBROUTINE rb_trn_loop( loop, trn )
!-------------------------------------------------------------------------------
!
!     Get trn(mx,gmy,itrn,is) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,0:ntrn-1,0:ns-1) :: trn

  !real(kind=DP), dimension(-nx:nx,0:ny,0:ntrn-1) :: wktrn
  !integer :: inum, irec, skipbyte, ir, rankw
  integer :: inum, irec
  !integer :: mx, my, gmy, itrn, is

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_trn
  integer(kind=4) :: start_trn(1:5), count_trn(1:5)

    call rb_trn_loop2inumirec( loop, inum, irec )

    !do is = 0, ns-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    skipbyte = recl_trn*(irec-1) + nhead + DP
    !    read( unit=500000000+100000*inum+ir, pos=skipbyte+1 ) wktrn
    !    do itrn = 0, ntrn-1
    !      do my = 0, ny
    !        gmy = ( ny+1 ) * rankw + my
    !        if ( gmy <= global_ny ) then
    !          do mx = -nx, nx
    !            trn(mx,gmy,itrn,is) = wktrn(mx,my,itrn)
    !          end do
    !        end if
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( trn_nc(inum), "trn", varid_trn )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_trn(:) = int((/ 2*nx+1,global_ny+1,ntrn,ns,1 /), kind=4)
    start_trn(:) = int((/ 1,1,1,1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( trn_nc(inum), varid_trn, &
                              values=trn(:,:,:,:), &
                              start=start_trn, count=count_trn )
    call check_nf90err(ierr_nf90, "nf90_get_var")

END SUBROUTINE rb_trn_loop


SUBROUTINE rb_trn_mxmyitrnis( mx, gmy, itrn, is, trn )
!-------------------------------------------------------------------------------
!
!     Get trn(loop) at mx, gmy, itrn, is
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, itrn, is
  real(kind=DP), intent(out), &
    dimension(loop_trn_sta(snum):loop_trn_end(enum)) :: trn

  integer :: loop, inum, irec, skipbyte, ir, rankw
  integer :: my

    call rb_gmy2rankwmy( gmy, rankw, my )

    ir = rankw + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    do inum = snum, enum
      irec = 1
      do loop = loop_trn_sta(inum), loop_trn_end(inum)
        skipbyte = recl_trn*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*itrn + (2*nx+1)*my + (mx+nx) )*DP
        read( unit=500000000+100000*inum+ir,  &
              pos=skipbyte+1 ) trn(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_trn_mxmyitrnis


!SUBROUTINE rb_trn_debug
!!-------------------------------------------------------------------------------
!!
!!     Debug rb_trn routines
!!                                                   (S. Maeyama, 24 Aug 2014)
!!
!!-------------------------------------------------------------------------------
!  real(kind=DP), dimension(:), allocatable :: wra
!  real(kind=DP) :: wr0, &
!                   wr2(-nx:nx,0:global_ny), &
!                   wr3(-nx:nx,0:global_ny,0:ntrn-1), &
!                   wr4(-nx:nx,0:global_ny,0:ntrn-1,0:ns-1)
!  integer :: mx, gmy, itrn, is, loop
!
!    allocate(wra(loop_trn_sta(snum):loop_trn_end(enum)))
!
!    write(*,*) "# test gettime and getalltime."
!    call rb_trn_getalltime( wra )
!    do loop = loop_trn_sta(snum), loop_trn_end(enum), nloop_trn/4
!      call rb_trn_gettime( loop, wr0 )
!      write(*,*) "loop = ", loop, ", time = ", wr0, wra(loop)
!    end do
!
!    !do loop = loop_trn_sta(snum), loop_trn_end(enum), nloop_trn/4
!    !  is = 1; itrn = 3; gmy = 0; mx = -nx
!    !  call rb_trn_mxmyitrnisloop( mx, gmy, itrn, is, loop, wr0 )
!    !  write(*,"(a,i7,a,g17.7e3)") "loop = ", loop, ", entropy_e = ", wr0
!    !end do
!
!    write(*,*) "# test mxmyitrnisloop, itrnisloop, isloop, and loop."
!    do loop = loop_trn_sta(snum), loop_trn_end(enum), nloop_trn/2
!      call rb_trn_loop( loop, wr4 )
!      do is = 0, ns-1
!        call rb_trn_isloop( is, loop, wr3 )
!        do itrn = 0, ntrn-1, ntrn/2
!          call rb_trn_itrnisloop( itrn, is, loop, wr2 )
!          do gmy = 0, global_ny, global_ny/2
!            do mx = -nx, nx, nx
!              call rb_trn_mxmyitrnisloop( mx, gmy, itrn, is, loop, wr0 )
!              write(*,"(4g17.7e3)") wr0, wr2(mx,gmy), wr3(mx,gmy,itrn), wr4(mx,gmy,itrn,is)
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    write(*,*) "# test mxmyitrnis."
!    do is = 0, ns-1
!      do itrn = 0, ntrn-1, ntrn/2
!        do gmy = 0, global_ny, global_ny/2
!          do mx = -nx, nx, nx
!            call rb_trn_mxmyitrnis( mx, gmy, itrn, is, wra )
!            do loop = loop_trn_sta(snum), loop_trn_end(enum), nloop_trn/2
!              call rb_trn_mxmyitrnisloop( mx, gmy, itrn, is, loop, wr0 )
!              write(*,"(4g17.7e3)") wr0, wra(loop)
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    deallocate(wra)
!
!    stop
!
!END SUBROUTINE rb_trn_debug


SUBROUTINE rb_tri_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *tri*
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, it, mxt, myt
  integer :: inum, ranks, it, mxt, myt
  character(len=4) :: cmx, cmy
  !character(len=1) :: srank
  character(len=3) :: cnum
  character(len=9) :: z_bound, z_filt, z_calc
  real(kind=DP) :: art_diff
  logical :: init_random

  integer(kind=4) :: ierr_nf90

  namelist /calct/ calc_type, z_bound, z_filt, z_calc, art_diff, &
                   init_random, num_triad_diag
  namelist /triad/ mxt, myt

    write( cnum, fmt="(i3.3)" ) snum
    open(inml, file="../gkvp_namelist."//cnum, status="old", action="read")
      
      calc_type="";z_bound="";z_filt="";z_calc="";art_diff=0._DP
      read(inml, nml=calct)
      if ( trim(calc_type) == "nonlinear" .and. num_triad_diag > 0 ) then
        allocate(triad_diag_mxt(0:num_triad_diag-1))
        allocate(triad_diag_myt(0:num_triad_diag-1))
        do it = 0, num_triad_diag-1
          read(inml, nml=triad)
          triad_diag_mxt(it) = mxt
          triad_diag_myt(it) = myt
        end do
      end if
    close(inml)
  if ( trim(calc_type) == "nonlinear" .and. num_triad_diag > 0 ) then
    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do ranks = 0, nprocs-1
      !  write( srank, fmt="(i1.1)" ) ranks
      !  do it = 0, num_triad_diag-1
      !    write( cmx, fmt="(i4.4)" ) triad_diag_mxt(it)
      !    write( cmy, fmt="(i4.4)" ) triad_diag_myt(it)
      !    open( unit=700000000+100000*inum+num_triad_diag*ranks+it,          &
      !          file="../phi/gkvp.s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnum, &
      !          status="old", action="read",                        &
      !          form="unformatted", access="stream" )
                                      !--- for debug ---
                                      !write(*,*) 700000000+100000*num_triad_diag*ranks+it, &
                                      !           "../phi/gkvp.s"     &
                                      !   //srank//"mx"//cmx//"my"//cmy//".tri."//cnum
                                      !-----------------
      !  end do
      !end do

      write( cmx, fmt="(i4.4)" ) triad_diag_mxt(0)
      write( cmy, fmt="(i4.4)" ) triad_diag_myt(0) 
      ierr_nf90 = nf90_open( path="../phi/gkvp.mx"//cmx//"my"//cmy//".tri."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=tri_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do
  end if

END SUBROUTINE rb_tri_fileopen


SUBROUTINE rb_tri_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *tri*
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, it
  integer :: inum
  integer(kind=4) :: ierr_nf90

  if ( trim(calc_type) == "nonlinear" .and. num_triad_diag > 0 ) then
    do inum = snum, enum
      !do ranks = 0, nprocs-1
      !  do it = 0, num_triad_diag-1
      !    close( unit=700000000+100000*inum+num_triad_diag*ranks+it )
      !  end do
      !end do

      ierr_nf90 = nf90_close( tri_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do
  end if

END SUBROUTINE rb_tri_fileclose


SUBROUTINE rb_tri_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_itri, varid_is
  integer(kind=4) :: len_kx, len_ky, len_itri, len_is

  if ( trim(calc_type) == "nonlinear" .and. num_triad_diag > 0 ) then
    do inum = snum, enum
      !read( unit=700000000+100000*inum+0, pos=1 ) head
      !read( unit=700000000+100000*inum+0, pos=recl_tri-4+1 ) foot

      !if ( head /= foot .or. head /= recl_tri - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *tri*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_tri-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if

      ierr_nf90 = nf90_inq_varid( tri_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( tri_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( tri_nc(inum), "itri", varid_itri )
      ierr_nf90 = nf90_inq_varid( tri_nc(inum), "is", varid_is )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( tri_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( tri_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( tri_nc(inum), varid_itri, len=len_itri )
      ierr_nf90 = nf90_inquire_dimension( tri_nc(inum), varid_is, len=len_is )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= 2*global_ny+1 .or. &
           len_itri /= ntri .or. len_is /= nprocs ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.tri*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#   2*global_ny+1 = ", 2*global_ny+1
        write(*,*) "#            itri = ", len_itri
        write(*,*) "#            ntri = ", ntri
        write(*,*) "#              is = ", len_is
        write(*,*) "#          nprocs = ", nprocs
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *tri* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_tri-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
                                        !--- for debug ---
                                        !write(*,*) head, recl_tri-2*nhead, foot
                                        !-----------------
    write(*,*) "# The netCDF files gkvp.tri*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#   2*global_ny+1 = ", 2*global_ny+1
    write(*,*) "#            itri = ", len_itri
    write(*,*) "#            ntri = ", ntri
    write(*,*) "#              is = ", len_is
    write(*,*) "#          nprocs = ", nprocs
  end if

END SUBROUTINE rb_tri_check


SUBROUTINE rb_tri_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !character(len=4) :: cmx, cmy
  !character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

  if ( trim(calc_type) == "nonlinear" .and. num_triad_diag > 0 ) then
    !write( cmx, fmt="(i4.4)" ) triad_diag_mxt(0)
    !write( cmy, fmt="(i4.4)" ) triad_diag_myt(0)

    nloop_tri = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_tri_sta(inum) = nloop_tri + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../phi/gkvp.s0mx"//cmx//"my"//cmy//".tri."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../phi/gkvp.s0mx"//cmx//"my"//cmy//".tri."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_tri_end(inum) = nloop_tri + filesize/recl_tri
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_tri_sta = ",  &
                                        !                    loop_tri_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_tri_end = ",  &
                                        !                    loop_tri_end(inum)
                                        !-----------------
      !nloop_tri = loop_tri_end(inum)

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( tri_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( tri_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_tri_end(inum) = nloop_tri + t_count

      nloop_tri = loop_tri_end(inum)
    end do
    write(*,*) "# nloop_tri = ", nloop_tri
  end if

END SUBROUTINE rb_tri_setnloop


SUBROUTINE rb_tri_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 6 Nov 2016)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_tri_sta(wknum) <= loop .and. loop <= loop_tri_end(wknum) ) then
        inum = wknum
        irec = loop - loop_tri_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_tri_loop2inumirec


SUBROUTINE rb_tri_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value

    call rb_tri_loop2inumirec( loop, inum, irec )

    !skipbyte = recl_tri*(irec-1) + nhead
    !read( unit=700000000+100000*inum+0, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( tri_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=tri_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)

END SUBROUTINE rb_tri_gettime


SUBROUTINE rb_tri_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 6 Nov 2016)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_tri_sta(snum):loop_tri_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte

    do inum = snum, enum
      irec = 1
      do loop = loop_tri_sta(inum), loop_tri_end(inum)
        skipbyte = recl_tri*(irec-1) + nhead
        read( unit=700000000+100000*inum+0,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_tri_getalltime


SUBROUTINE rb_tri_mxtmytisloop( mxt, myt, is, loop, tri )
!-------------------------------------------------------------------------------
!
!     Get tri(mx,my) of mxt, myt, is, at loop
!                                                   (S. Maeyama, 6 Nov 2016)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mxt, myt, is, loop
  real(kind=DP), intent(out), &
    dimension(-nx:nx,-global_ny:global_ny,0:ntri-1) :: tri

  !integer :: inum, irec, skipbyte, iw, it
  integer :: inum, irec, iw, it

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tri
  integer(kind=4) :: start_tri(1:5), count_tri(1:5)
  real(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:ntri-1) :: ff

    it = -1
    do iw = 0, num_triad_diag-1
      if (mxt == triad_diag_mxt(iw) .and. myt == triad_diag_myt(iw)) it = iw
    end do
    if (it == -1) then
      write(*,*) "# In rb_tri_mxtmytisloop, invalid mxt, myt.", mxt, myt
      stop
    end if
    call rb_tri_loop2inumirec( loop, inum, irec )

    !skipbyte = recl_tri*(irec-1) &  ! (irec-1) lines are skipped.
    !         + nhead + DP           ! header and time is also skipped.
    !read( unit=700000000+100000*inum+num_triad_diag*is+it,  &
    !      pos=skipbyte+1 ) tri

    ierr_nf90 = nf90_inq_varid( tri_nc(inum), "tri", varid_tri )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_tri(:) = int((/ 2*nx+1,2*global_ny+1,ntri,1,1 /), kind=4)
    start_tri(:) = int((/ 1,1,1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( tri_nc(inum), varid_tri, &
                              values=tri(:,:,:), &
                              start=start_tri, count=count_tri )
    call check_nf90err(ierr_nf90, "nf90_get_var")

END SUBROUTINE rb_tri_mxtmytisloop


SUBROUTINE rb_fxv_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *fxv*
!                                                   (S. Maeyama, 5 May 2020)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankm, rankv, rankz, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  !character(len=1) :: srank
  character(len=3) :: cnum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do ranks = 0, nprocs-1
      !  write( srank, fmt="(i1.1)" ) ranks
      !  do rankm = 0, nprocm-1
      !    do rankv = 0, nprocv-1
      !      do rankz = 0, nprocz-1
      !        do rankw = 0, nprocw-1
      !          ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
      !             + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*ranks
      !           write( crank, fmt="(i6.6)" ) ir
      !           open( unit=800000000+100000*inum+ir,                      &
      !                 file="../fxv/gkvp."//crank//"."//srank//".fxv."//cnum, &
      !                 status="old", action="read",                        &
      !                 form="unformatted", access="stream" )
                                        !--- for debug ---
                                        !write(*,*) 600000000+100000*inum+ir, &
                                        !           "../fxv/gkvp."      &
                                        !   //crank//"."//srank//".fxv."//cnum
                                        !-----------------
      !        end do
      !      end do
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_open( path="../fxv/gkvp.fxv."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=fxv_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_fxv_fileopen


SUBROUTINE rb_fxv_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *fxv*
!                                                   (S. Maeyama, 5 May 2020)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer :: inum, ranks, rankm, rankv, rankz, rankw, ir
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do ranks = 0, nprocs-1
      !  do rankm = 0, nprocm-1
      !    do rankv = 0, nprocv-1
      !      do rankz = 0, nprocz-1
      !        do rankw = 0, nprocw-1
      !          ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
      !             + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*ranks
      !           close( unit=800000000+100000*inum+ir )
      !        end do
      !      end do
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_close( fxv_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_fxv_fileclose


SUBROUTINE rb_fxv_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 5 May 2020)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir, ierr
  integer :: inum
  !character(len=3) :: cnum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is
  integer(kind=4) :: len_kx, len_ky, len_zz, len_vl, len_mu, len_is

    do inum = snum, enum
      !ir = 0

      !write( cnum, fmt="(i3.3)" ) inum
      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../fxv/gkvp.000000.0.fxv."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../fxv/gkvp.000000.0.fxv."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !if (filesize == 0) then
      !  write(*,*) "# WARNING: Filesize is zero for reading the binary files *fxv*."
      !  write(*,*) "# inum = ", inum
      !  write(*,*) "# filesize = ", filesize
      !else

        !read( unit=800000000+100000*inum+ir, pos=1 ) head
        !read( unit=800000000+100000*inum+ir, pos=recl_fxv-4+1 ) foot
      !  read( unit=800000000+100000*inum+ir, pos=1, iostat=ierr ) head
      !  read( unit=800000000+100000*inum+ir, pos=recl_fxv-4+1, iostat=ierr ) foot
      !  if (ierr /= 0 .and. ierr /= -1) then
      !    write(*,*) "# Error for reading the binary files *fxv*."
      !    write(*,*) "# inum = ", inum
      !    write(*,*) "# ierr = ", ierr
      !    stop
      !  end if
  
      !  if ( head /= foot .or. head /= recl_fxv - 2*nhead ) then
      !    write(*,*) "# Error for reading the binary files *fxv*."
      !    write(*,*) "#         head = ", head
      !    write(*,*) "# recl-2*nhead = ", recl_fxv-2*nhead, ", where nhead = ", nhead
      !    write(*,*) "#         foot = ", foot
      !    write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !    stop
      !  end if

      !end if

      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "zz", varid_zz )
      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "vl", varid_vl )
      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "mu", varid_mu )
      ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "is", varid_is )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_zz, len=len_zz )
      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_vl, len=len_vl )
      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_mu, len=len_mu )
      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), varid_is, len=len_is )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           !len_zz /= nprocz/2+1 .or. len_vl /= 2*global_nv .or. &
           len_zz /= nprocz .or. len_vl /= 2*global_nv .or. &
           len_mu /= global_nm+1 .or. len_is /= nprocs ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.fxv*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#              zz = ", len_zz
        !write(*,*) "#      nprocz/2+1 = ", nprocz/2+1
        write(*,*) "#          nprocz = ", nprocz
        write(*,*) "#              vl = ", len_vl
        write(*,*) "#     2*global_nv = ", 2*global_nv
        write(*,*) "#              mu = ", len_mu
        write(*,*) "#     global_nm+1 = ", global_nm+1
        write(*,*) "#              is = ", len_is
        write(*,*) "#          nprocs = ", nprocs
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *fxv* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_fxv-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
                                        !--- for debug ---
                                        !write(*,*) head, recl_fxv-2*nhead, foot
                                        !-----------------
    write(*,*) "# The netCDF files gkvp.fxv*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#              zz = ", len_zz
    !write(*,*) "#      nprocz/2+1 = ", nprocz/2+1
    write(*,*) "#          nprocz = ", nprocz
    write(*,*) "#              vl = ", len_vl
    write(*,*) "#     2*global_nv = ", 2*global_nv
    write(*,*) "#              mu = ", len_mu
    write(*,*) "#     global_nm+1 = ", global_nm+1
    write(*,*) "#              is = ", len_is
    write(*,*) "#          nprocs = ", nprocs

END SUBROUTINE rb_fxv_check


SUBROUTINE rb_fxv_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 5 May 2020)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid


    nloop_fxv = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_fxv_sta(inum) = nloop_fxv + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../fxv/gkvp.000000.0.fxv."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../fxv/gkvp.000000.0.fxv."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_fxv_end(inum) = nloop_fxv + filesize/recl_fxv
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_fxv_sta = ",  &
                                        !                    loop_fxv_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_fxv_end = ",  &
                                        !                    loop_fxv_end(inum)
                                        !-----------------
      !nloop_fxv = loop_fxv_end(inum)

      !if (filesize == 0) then
      !  write(*,*) "# WARNING: Filesize is zero for reading the binary files *fxv*."
      !  write(*,*) "# inum = ", inum
      !  write(*,*) "# filesize = ", filesize
      !  write(*,*) "# loop_fxv_sta(inum) = ", loop_fxv_sta(inum)
      !  write(*,*) "# loop_fxv_end(inum) = ", loop_fxv_end(inum)
      !  write(*,*) "# Subroutine rb_fxv_loop2inumirec will still work."
      !end if

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( fxv_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( fxv_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_fxv_end(inum) = nloop_fxv + t_count

      nloop_fxv = loop_fxv_end(inum)
    end do
    write(*,*) "# nloop_fxv = ", nloop_fxv
                                        !--- for debug ---
                                        !write(*,*) "# loop_fxv_sta = ", loop_fxv_sta(:)
                                        !write(*,*) "# loop_fxv_end = ", loop_fxv_end(:)
                                        !do loop = 0, nloop_fxv
                                        !  call rb_fxv_loop2inumirec(loop,inum,irec)
                                        !  write(*,*) loop, inum, irec
                                        !end do
                                        !-----------------

END SUBROUTINE rb_fxv_setnloop


SUBROUTINE rb_fxv_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_fxv_sta(wknum) <= loop .and. loop <= loop_fxv_end(wknum) ) then
        inum = wknum
        irec = loop - loop_fxv_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_fxv_loop2inumirec


SUBROUTINE rb_fxv_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value

    call rb_fxv_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_fxv*(irec-1) + nhead
    !read( unit=800000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=fxv_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "rb_fxv_gettime" )

    time = t_value(1)

END SUBROUTINE rb_fxv_gettime


SUBROUTINE rb_fxv_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_fxv_sta(snum):loop_fxv_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir

    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_fxv_sta(inum), loop_fxv_end(inum)
        skipbyte = recl_fxv*(irec-1) + nhead
        read( unit=800000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_fxv_getalltime


SUBROUTINE rb_fxv_mxmyrankzivimisloop( mx, gmy, rankz, giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, rankz, giv, gim, is, loop
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, rankz, giv, gim, is, loop
  complex(kind=DP), intent(out) :: ff

  integer :: inum, irec, skipbyte, ir, rankm, rankv, rankw
  integer :: my, iv, im

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giv2rankviv( giv, rankv, iv )
    call rb_gim2rankmim( gim, rankm, im )
    call rb_fxv_loop2inumirec( loop, inum, irec )

    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    skipbyte = recl_fxv*(irec-1) &  ! (irec-1) lines are skipped.
             + nhead + DP        &  ! header and time is also skipped.
             + ( (2*nx+1)*(ny+1)*(2*nv)*im + (2*nx+1)*(ny+1)*(iv-1) &
               + (2*nx+1)*my + (mx+nx) )*(2*DP)
    read( unit=800000000+100000*inum+ir,  &
          pos=skipbyte+1 ) ff

END SUBROUTINE rb_fxv_mxmyrankzivimisloop


SUBROUTINE rb_fxv_mxmyrankzisloop( mx, gmy, rankz, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, rankz, is, loop
!                                                   (S. Maeyama, 5 May 2020)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, rankz, is, loop
  complex(kind=DP), intent(out), dimension(1:2*global_nv,0:global_nm) :: ff

  !integer :: inum, irec, skipbyte, ir, rankm, rankv, rankw
  integer :: inum, irec
  !integer :: my, iv, im, giv, gim

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_refxv, varid_imfxv
  integer(kind=4) :: start_fxv(1:7), count_fxv(1:7)
  real(kind=DP), dimension(1:2*global_nv,0:global_nm) :: refxv
  real(kind=DP), dimension(1:2*global_nv,0:global_nm) :: imfxv

    !call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_fxv_loop2inumirec( loop, inum, irec )

    !do rankm = 0, nprocm-1
    !  do rankv = 0, nprocv-1
    !    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
    !       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    do im = 0, nm
    !      gim = (nm+1) * rankm + im
    !      do iv = 1, 2*nv
    !        giv = (2*nv) * rankv + iv
    !        skipbyte = recl_fxv*(irec-1) &  ! (irec-1) lines are skipped.
    !                 + nhead + DP        &  ! header and time is also skipped.
    !                 + ( (2*nx+1)*(ny+1)*(2*nv)*im + (2*nx+1)*(ny+1)*(iv-1)     &
    !                   + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !        read( unit=800000000+100000*inum+ir,  &
    !              pos=skipbyte+1 ) ff(giv,gim)
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "refxv", varid_refxv )
    ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "imfxv", varid_imfxv )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_fxv(:) = int((/ 1,1,1,2*global_nv,global_nm+1,1,1 /), kind=4)
    start_fxv(:) = int((/ nx+mx+1,gmy+1,rankz+1,1,1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( fxv_nc(inum), varid_refxv, &
                              values=refxv(:,:), &
                              start=start_fxv, count=count_fxv )
    ierr_nf90 = nf90_get_var( fxv_nc(inum), varid_imfxv, &
                              values=imfxv(:,:), &
                              start=start_fxv, count=count_fxv )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    ff = cmplx(refxv,imfxv,kind=DP)

END SUBROUTINE rb_fxv_mxmyrankzisloop


SUBROUTINE rb_fxv_rankzivimisloop( rankz, giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at rankz, giv, gim, is, loop
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: rankz, giv, gim, is, loop
  complex(kind=DP), intent(out), dimension(-nx:nx,0:global_ny) :: ff

  !complex(kind=DP), dimension(-nx:nx,0:ny) :: wkff
  !integer :: inum, irec, skipbyte, ir, rankm, rankv, rankw
  integer :: inum, irec
  !integer :: mx, my, gmy, iv, im

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_refxv, varid_imfxv
  integer(kind=4) :: start_fxv(1:7), count_fxv(1:7)
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: refxv
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: imfxv

    !call rb_giv2rankviv( giv, rankv, iv )
    !call rb_gim2rankmim( gim, rankm, im )
    call rb_fxv_loop2inumirec( loop, inum, irec )

    !do rankw = 0, nprocw-1
    !  ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
    !     + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !  skipbyte = recl_fxv*(irec-1) &  ! (irec-1) lines are skipped.
    !           + nhead + DP        &  ! header and time is also skipped.
    !           + ( (2*nx+1)*(ny+1)*(2*nv)*im + (2*nx+1)*(ny+1)*(iv-1) )*(2*DP)
    !  read( unit=800000000+100000*inum+ir,  &
    !        pos=skipbyte+1 ) wkff
    !  do my = 0, ny
    !    gmy = ( ny+1 ) * rankw + my
    !    if ( gmy <= global_ny ) then
    !      do mx = -nx, nx
    !        ff(mx,gmy) = wkff(mx,my)
    !      end do
    !    end if
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "refxv", varid_refxv )
    ierr_nf90 = nf90_inq_varid( fxv_nc(inum), "imfxv", varid_imfxv )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_fxv(:) = int((/ 2*nx+1,global_ny+1,1,1,1,1,1 /), kind=4)
    start_fxv(:) = int((/ 1,1,rankz+1,giv,gim+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( fxv_nc(inum), varid_refxv, &
                              values=refxv(:,:), &
                              start=start_fxv, count=count_fxv )
    ierr_nf90 = nf90_get_var( fxv_nc(inum), varid_imfxv, &
                              values=imfxv(:,:), &
                              start=start_fxv, count=count_fxv )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    ff = cmplx(refxv,imfxv,kind=DP)

END SUBROUTINE rb_fxv_rankzivimisloop


SUBROUTINE rb_cnt_fileopen
!-------------------------------------------------------------------------------
!
!     Open all files of *cnt*
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankm, rankv, rankz, rankw, ir
  integer :: inum
  !character(len=6) :: crank
  character(len=3) :: cnum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      !do ranks = 0, nprocs-1
      !  do rankm = 0, nprocm-1
      !    do rankv = 0, nprocv-1
      !      do rankz = 0, nprocz-1
      !        do rankw = 0, nprocw-1
      !          ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
      !             + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*ranks
      !           write( crank, fmt="(i6.6)" ) ir
      !           open( unit=600000000+100000*inum+ir,                      &
      !                 file="../cnt/gkvp."//crank//".cnt."//cnum, &
      !                 status="old", action="read",                        &
      !                 form="unformatted", access="stream" )
                                        !--- for debug ---
                                        !write(*,*) 600000000+100000*inum+ir, &
                                        !           "../cnt/gkvp."      &
                                        !   //crank//".cnt."//cnum
                                        !-----------------
      !        end do
      !      end do
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_open( path="../cnt/gkvp.cnt."//cnum//".nc",  &
                             mode=NF90_NOWRITE, ncid=cnt_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_open" )
    end do

END SUBROUTINE rb_cnt_fileopen


SUBROUTINE rb_cnt_fileclose
!-------------------------------------------------------------------------------
!
!     Close all files of *cnt*
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer :: inum, ranks, rankm, rankv, rankz, rankw, ir
  integer :: inum
  integer(kind=4) :: ierr_nf90

    do inum = snum, enum
      !do ranks = 0, nprocs-1
      !  do rankm = 0, nprocm-1
      !    do rankv = 0, nprocv-1
      !      do rankz = 0, nprocz-1
      !        do rankw = 0, nprocw-1
      !          ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
      !             + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*ranks
      !           close( unit=600000000+100000*inum+ir )
      !        end do
      !      end do
      !    end do
      !  end do
      !end do

      ierr_nf90 = nf90_close( cnt_nc(inum) )
      call check_nf90err( ierr_nf90, "nf90_close" )
    end do

END SUBROUTINE rb_cnt_fileclose


SUBROUTINE rb_cnt_check
!-------------------------------------------------------------------------------
!
!     Check reading a binary file (mainly endianness check)
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !integer(kind=nhead) :: head, foot
  !integer :: inum, ir
  integer :: inum

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is
  integer(kind=4) :: len_kx, len_ky, len_zz, len_vl, len_mu, len_is

    do inum = snum, enum
      !ir = 0
      !read( unit=600000000+100000*inum+ir, pos=1 ) head
      !read( unit=600000000+100000*inum+ir, pos=recl_cnt-4+1 ) foot

      !if ( head /= foot .or. head /= recl_cnt - 2*nhead ) then
      !  write(*,*) "# Error for reading the binary files *cnt*."
      !  write(*,*) "#         head = ", head
      !  write(*,*) "# recl-2*nhead = ", recl_cnt-2*nhead, ", where nhead = ", nhead
      !  write(*,*) "#         foot = ", foot
      !  write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
      !  stop
      !end if

      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "kx", varid_kx )
      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "ky", varid_ky )
      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "zz", varid_zz )
      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "vl", varid_vl )
      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "mu", varid_mu )
      ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "is", varid_is )
      call check_nf90err( ierr_nf90, "nf90_inq_varid" )

      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_kx, len=len_kx )
      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_ky, len=len_ky )
      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_zz, len=len_zz )
      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_vl, len=len_vl )
      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_mu, len=len_mu )
      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), varid_is, len=len_is )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      if ( len_kx /= 2*nx+1 .or. len_ky /= global_ny+1 .or. &
           len_zz /= 2*global_nz .or. len_vl /= 2*global_nv .or. &
           len_mu /= global_nm+1 .or. len_is /= nprocs ) then
        write(*,*) "# Error for reading the netCDF files *gkvp.cnt*.nc."
        write(*,*) "#   array size : "
        write(*,*) "#              kx = ", len_kx
        write(*,*) "#          2*nx+1 = ", 2*nx+1
        write(*,*) "#              ky = ", len_ky
        write(*,*) "#     global_ny+1 = ", global_ny+1
        write(*,*) "#              zz = ", len_zz
        write(*,*) "#     2*global_nz = ", 2*global_nz
        write(*,*) "#              vl = ", len_vl
        write(*,*) "#     2*global_nv = ", 2*global_nv
        write(*,*) "#              mu = ", len_mu
        write(*,*) "#     global_nm+1 = ", global_nm+1
        write(*,*) "#              is = ", len_is
        write(*,*) "#          nprocs = ", nprocs
        write(*,*) "# Check settings, e.g., grid size, nhead, Endian !"
        stop
      end if
    end do
    !write(*,*) "# The binary files *cnt* are correctly read."
    !write(*,*) "#         head = ", head
    !write(*,*) "# recl-2*nhead = ", recl_cnt-2*nhead, ", where nhead = ", nhead
    !write(*,*) "#         foot = ", foot
                                        !--- for debug ---
                                        !write(*,*) head, recl_cnt-2*nhead, foot
                                        !-----------------
    write(*,*) "# The netCDF files gkvp.cnt*nc are correctly read."
    write(*,*) "#   array size : "
    write(*,*) "#              kx = ", len_kx
    write(*,*) "#          2*nx+1 = ", 2*nx+1
    write(*,*) "#              ky = ", len_ky
    write(*,*) "#     global_ny+1 = ", global_ny+1
    write(*,*) "#              zz = ", len_zz
    write(*,*) "#     2*global_nz = ", 2*global_nz
    write(*,*) "#              vl = ", len_vl
    write(*,*) "#     2*global_nv = ", 2*global_nv
    write(*,*) "#              mu = ", len_mu
    write(*,*) "#     global_nm+1 = ", global_nm+1
    write(*,*) "#              is = ", len_is
    write(*,*) "#          nprocs = ", nprocs

END SUBROUTINE rb_cnt_check


SUBROUTINE rb_cnt_setnloop
!-------------------------------------------------------------------------------
!
!     Set total number of record, start and end for each inum
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  !character(len=3) :: cnum
  integer :: inum
  !integer :: filesize, statb(13)

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: t_count
  integer(kind=4) :: ndims, nvars, ngatts, unlimdimid

    nloop_cnt = -1
    do inum = 1, enum
      !write( cnum, fmt="(i3.3)" ) inum

      loop_cnt_sta(inum) = nloop_cnt + 1

      !if (trim(flag_getsize) == "stat") then
      !  call stat( "../cnt/gkvp.000000.cnt."//cnum, statb )
      !  filesize = statb(8)
      !else if (trim(flag_getsize) == "inquire") then
      !  inquire( file="../cnt/gkvp.000000.cnt."//cnum, size=filesize )
      !else
      !  write(*,*) "Wrong flag_getsize: ", flag_getsize
      !end if

      !loop_cnt_end(inum) = nloop_cnt + filesize/recl_cnt
                                        !--- for debug ---
                                        !write(*,"(a)") "# inum = "//cnum
                                        !write(*,"(a,i8.8)") "# loop_cnt_sta = ",  &
                                        !                    loop_cnt_sta(inum)
                                        !write(*,"(a,i8.8)") "# loop_cnt_end = ",  &
                                        !                    loop_cnt_end(inum)
                                        !-----------------
      !nloop_cnt = loop_cnt_end(inum)

      != Information about an open netCDF dataset
      ierr_nf90 = nf90_inquire( cnt_nc(inum), &
                                ndims, nvars, ngatts, unlimdimid )
      call check_nf90err( ierr_nf90, "nf90_inquire" )

      ierr_nf90 = nf90_inquire_dimension( cnt_nc(inum), unlimdimid, &
                                          len=t_count )
      call check_nf90err( ierr_nf90, "nf90_inquire_dimension" )

      loop_cnt_end(inum) = nloop_cnt + t_count

      nloop_cnt = loop_cnt_end(inum)
    end do
    write(*,*) "# nloop_cnt = ", nloop_cnt

END SUBROUTINE rb_cnt_setnloop


SUBROUTINE rb_cnt_loop2inumirec( loop, inum, irec )
!-------------------------------------------------------------------------------
!
!     Calculate inum and irec from loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  integer, intent(out) :: inum, irec

  integer :: wknum

    irec = -1
    inum = -1
    do wknum = snum, enum
      if ( loop_cnt_sta(wknum) <= loop .and. loop <= loop_cnt_end(wknum) ) then
        inum = wknum
        irec = loop - loop_cnt_sta(inum) + 1
        exit
      end if
    end do

END SUBROUTINE rb_cnt_loop2inumirec


SUBROUTINE rb_cnt_gettime( loop, time )
!-------------------------------------------------------------------------------
!
!     Get time at loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: loop
  real(kind=DP), intent(out) :: time

  !integer :: inum, irec, skipbyte, ir
  integer :: inum, irec

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_tt
  integer(kind=4) :: start_time(1),  count_time(1)
  real(kind=DP), dimension(1) :: t_value

    call rb_cnt_loop2inumirec( loop, inum, irec )

    !ir = 0
    !skipbyte = recl_cnt*(irec-1) + nhead
    !read( unit=600000000+100000*inum+ir, pos=skipbyte+1 ) time

    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "t", varid_tt )
    call check_nf90err( ierr_nf90, "nf90_inq_varid" )

    count_time(:) = 1
    start_time(:) = irec
    ierr_nf90 = nf90_get_var( ncid=cnt_nc(inum), varid=varid_tt, &
                              values=t_value, &
                              start=start_time, count=count_time )
    call check_nf90err( ierr_nf90, "nf90_get_var" )

    time = t_value(1)

END SUBROUTINE rb_cnt_gettime


SUBROUTINE rb_cnt_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(loop_cnt_sta(snum):loop_cnt_end(enum)), intent(out) :: time

  integer :: loop, inum, irec, skipbyte, ir

    ir = 0
    do inum = snum, enum
      irec = 1
      do loop = loop_cnt_sta(inum), loop_cnt_end(inum)
        skipbyte = recl_cnt*(irec-1) + nhead
        read( unit=600000000+100000*inum+ir,  &
              pos=skipbyte+1 ) time(loop)
        irec = irec + 1
      end do
    end do

END SUBROUTINE rb_cnt_getalltime


SUBROUTINE rb_cnt_mxmyizivimisloop( mx, gmy, giz, giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, giz, giv, gim, is, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, giv, gim, is, loop
  complex(kind=DP), intent(out) :: ff

  integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: my, iz, iv, im

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )
    call rb_giv2rankviv( giv, rankv, iv )
    call rb_gim2rankmim( gim, rankm, im )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
             + nhead + DP        &  ! header and time is also skipped.
             + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im + (2*nx+1)*(ny+1)*(2*nz)*(iv-1) &
               + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    read( unit=600000000+100000*inum+ir,  &
          pos=skipbyte+1 ) ff

END SUBROUTINE rb_cnt_mxmyizivimisloop


SUBROUTINE rb_cnt_mxmyizisloop( mx, gmy, giz, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, giz, is, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giz, is, loop
  complex(kind=DP), intent(out),  &
    dimension(1:2*global_nv,0:global_nm) :: ff

  integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: my, iz, iv, im, giv, gim

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giz2rankziz( giz, rankz, iz )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    do rankm = 0, nprocm-1
      do rankv = 0, nprocv-1
        ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
           + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
        do im = 0, nm
          gim = (nm+1) * rankm + im
          do iv = 1, 2*nv
            giv = (2*nv) * rankv + iv
            skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
                     + nhead + DP        &  ! header and time is also skipped.
                     + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
                       + (2*nx+1)*(ny+1)*(2*nz)*(iv-1)     &
                       + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
            read( unit=600000000+100000*inum+ir,  &
                  pos=skipbyte+1 ) ff(giv,gim)
          end do
        end do
      end do
    end do

END SUBROUTINE rb_cnt_mxmyizisloop


SUBROUTINE rb_cnt_mxmyimisloop( mx, gmy, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, gim, is, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!                                                   (Fujitsu,    22 Jan. 2021)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, gim, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv) :: ff

  !integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: inum, irec, rankm, rankw
  !integer :: my, iz, iv, im, giz, giv

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_recnt, varid_imcnt
  integer(kind=4) :: start_cnt(1:7), count_cnt(1:7)
  real(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv) :: recnt
  real(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv) :: imcnt

    !call rb_gmy2rankwmy( gmy, rankw, my )
    !call rb_gim2rankmim( gim, rankm, im )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    !do rankv = 0, nprocv-1
    !  do rankz = 0, nprocz-1
    !    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
    !       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    do iv = 1, 2*nv
    !      giv = (2*nv) * rankv + iv
    !      do iz = -nz, nz-1
    !        giz = - global_nz + 2*nz * rankz + iz + nz
    !        skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
    !                 + nhead + DP        &  ! header and time is also skipped.
    !                 + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
    !                   + (2*nx+1)*(ny+1)*(2*nz)*(iv-1)     &
    !                   + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
    !        read( unit=600000000+100000*inum+ir,  &
    !              pos=skipbyte+1 ) ff(giz,giv)
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "recnt", varid_recnt )
    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "imcnt", varid_imcnt )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_cnt(:) = int((/ 1,1,2*global_nz,2*global_nv,1,1,1 /), kind=4)
    start_cnt(:) = int((/ nx+mx+1,gmy+1,1,1,gim+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_recnt, &
                              values=recnt(:,:), &
                              start=start_cnt, count=count_cnt )
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_imcnt, &
                              values=imcnt(:,:), &
                              start=start_cnt, count=count_cnt )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    ff = cmplx(recnt,imcnt,kind=DP)

END SUBROUTINE rb_cnt_mxmyimisloop


SUBROUTINE rb_cnt_mxmyisloop( mx, gmy, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, is, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) :: ff

  integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: my, iz, iv, im, giz, giv, gim

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    do rankm = 0, nprocm-1
      do rankv = 0, nprocv-1
        do rankz = 0, nprocz-1
          ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
             + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
          do im = 0, nm
            gim = (nm+1) * rankm + im
            do iv = 1, 2*nv
              giv = (2*nv) * rankv + iv
              do iz = -nz, nz-1
                giz = - global_nz + 2*nz * rankz + iz + nz
                skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
                         + nhead + DP        &  ! header and time is also skipped.
                         + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
                           + (2*nx+1)*(ny+1)*(2*nz)*(iv-1)     &
                           + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
                read( unit=600000000+100000*inum+ir,  &
                      pos=skipbyte+1 ) ff(giz,giv,gim)
              end do
            end do
          end do
        end do
      end do
    end do

END SUBROUTINE rb_cnt_mxmyisloop


SUBROUTINE rb_cnt_ivimisloop( giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at giv, gim, is, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giv, gim, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: ff

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wkff
  !integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, iv, im, gmy, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_recnt, varid_imcnt
  integer(kind=4) :: start_cnt(1:7), count_cnt(1:7)
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: recnt
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: imcnt

    !call rb_giv2rankviv( giv, rankv, iv )
    !call rb_gim2rankmim( gim, rankm, im )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
    !       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
    !             + nhead + DP        &  ! header and time is also skipped.
    !             + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
    !               + (2*nx+1)*(ny+1)*(2*nz)*(iv-1) )*(2*DP)
    !    read( unit=600000000+100000*inum+ir, pos=skipbyte+1 ) wkff
    !    do iz = -nz, nz-1
    !      giz = - global_nz + 2*nz * rankz + iz + nz
    !      do my = 0, ny
    !        gmy = ( ny+1 ) * rankw + my
    !        if ( gmy <= global_ny ) then
    !          do mx = -nx, nx
    !            ff(mx,gmy,giz) = wkff(mx,my,iz)
    !          end do
    !        end if
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "recnt", varid_recnt )
    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "imcnt", varid_imcnt )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_cnt(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,1,1,1,1 /), kind=4)
    start_cnt(:) = int((/ 1,1,1,giv,gim+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_recnt, &
                              values=recnt(:,:,:), &
                              start=start_cnt, count=count_cnt )
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_imcnt, &
                              values=imcnt(:,:,:), &
                              start=start_cnt, count=count_cnt )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    ff = cmplx(recnt,imcnt,kind=DP)

END SUBROUTINE rb_cnt_ivimisloop


SUBROUTINE rb_cnt_izivimisloop( giz, giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at giz, giv, gim, is, loop
!                                                   (S. Maeyama, 12 March 2022)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz, giv, gim, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-nx:nx,0:global_ny) :: ff

  !complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wkff
  !integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: inum, irec
  !integer :: mx, my, iz, iv, im, gmy, giz

  integer(kind=4) :: ierr_nf90
  integer(kind=4) :: varid_recnt, varid_imcnt
  integer(kind=4) :: start_cnt(1:7), count_cnt(1:7)
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: recnt
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: imcnt

    !call rb_giv2rankviv( giv, rankv, iv )
    !call rb_gim2rankmim( gim, rankm, im )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    !do rankz = 0, nprocz-1
    !  do rankw = 0, nprocw-1
    !    ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
    !       + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
    !    skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
    !             + nhead + DP        &  ! header and time is also skipped.
    !             + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
    !               + (2*nx+1)*(ny+1)*(2*nz)*(iv-1) )*(2*DP)
    !    read( unit=600000000+100000*inum+ir, pos=skipbyte+1 ) wkff
    !    do iz = -nz, nz-1
    !      giz = - global_nz + 2*nz * rankz + iz + nz
    !      do my = 0, ny
    !        gmy = ( ny+1 ) * rankw + my
    !        if ( gmy <= global_ny ) then
    !          do mx = -nx, nx
    !            ff(mx,gmy,giz) = wkff(mx,my,iz)
    !          end do
    !        end if
    !      end do
    !    end do
    !  end do
    !end do

    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "recnt", varid_recnt )
    ierr_nf90 = nf90_inq_varid( cnt_nc(inum), "imcnt", varid_imcnt )
    call check_nf90err(ierr_nf90, "nf90_inq_varid" )

    count_cnt(:) = int((/ 2*nx+1,global_ny+1,1,1,1,1,1 /), kind=4)
    start_cnt(:) = int((/ 1,1,global_nz+giz+1,giv,gim+1,is+1,irec /), kind=4)
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_recnt, &
                              values=recnt(:,:), &
                              start=start_cnt, count=count_cnt )
    ierr_nf90 = nf90_get_var( cnt_nc(inum), varid_imcnt, &
                              values=imcnt(:,:), &
                              start=start_cnt, count=count_cnt )
    call check_nf90err(ierr_nf90, "nf90_get_var")

    ff = cmplx(recnt,imcnt,kind=DP)

END SUBROUTINE rb_cnt_izivimisloop


SUBROUTINE rb_cnt_mxmyivimisloop( mx, gmy, giv, gim, is, loop, ff )
!-------------------------------------------------------------------------------
!
!     Get ff at mx, gmy, giv, gim, is, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: mx, gmy, giv, gim, is, loop
  complex(kind=DP), intent(out),  &
    dimension(-global_nz:global_nz-1) :: ff

  integer :: inum, irec, skipbyte, ir, rankm, rankv, rankz, rankw
  integer :: my, iv, im, iz, giz

    call rb_gmy2rankwmy( gmy, rankw, my )
    call rb_giv2rankviv( giv, rankv, iv )
    call rb_gim2rankmim( gim, rankm, im )
    call rb_cnt_loop2inumirec( loop, inum, irec )

    do rankz = 0, nprocz-1
      ir = rankw + nprocw*rankz + nprocw*nprocz*rankv  &
         + nprocw*nprocz*nprocv*rankm + nprocw*nprocz*nprocv*nprocm*is  ! ranks=is
      do iz = -nz, nz-1
        giz = - global_nz + 2*nz * rankz + iz + nz
        skipbyte = recl_cnt*(irec-1) &  ! (irec-1) lines are skipped.
                 + nhead + DP        &  ! header and time is also skipped.
                 + ( (2*nx+1)*(ny+1)*(2*nz)*(2*nv)*im  &
                   + (2*nx+1)*(ny+1)*(2*nz)*(iv-1)     &
                   + (2*nx+1)*(ny+1)*(iz+nz) + (2*nx+1)*my + (mx+nx) )*(2*DP)
        read( unit=600000000+100000*inum+ir, pos=skipbyte+1 ) ff(giz)
      end do
    end do

END SUBROUTINE rb_cnt_mxmyivimisloop


!--------------------------------------
SUBROUTINE check_nf90err(werr, comment)
!--------------------------------------
!  Check error message of nf90

  integer(kind=4), intent(in) :: werr
  character(len=*), intent(in) :: comment

    if(werr /= nf90_noerr) then
      write(*,*) comment//" "//trim(nf90_strerror(werr))
      stop
    end if

END SUBROUTINE check_nf90err


END MODULE diag_rb
