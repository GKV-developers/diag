MODULE out_ffinzvm_vtk
!-------------------------------------------------------------------------------
!
!     Output ff in (z,v,m)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fkinzvm_vtk, fkinzvm_connect_vtk

  !character(len=9), parameter :: flag_normalize = "none" ! as-is
  character(len=9), parameter :: flag_normalize = "phi0" ! phi(z=0) normalize
  !character(len=9), parameter :: flag_normalize = "Al0"  ! Al(z=0) normalize

 CONTAINS


SUBROUTINE fkinzvm_vtk( mx, gmy, is, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (zz,vl,mu) in VTK binary format
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyisloop,               &
                      rb_phi_gettime, rb_phi_mxmyizloop, loop_phi_end, &
                      rb_Al_gettime, rb_Al_mxmyizloop, loop_Al_end
  use diag_geom, only : kx, gky, n_tht, vmax, gvp, beta
  use diag_functions, only : besselj0

  integer, intent(in) :: mx, gmy, is, loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) :: fk
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iz, iv, im, loop

    call rb_phi_gettime( loop_phi_end(enum), time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop_phi_end(enum), phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop_Al_end(enum), Al0 )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cis, fmt="(i1.1)" ) is

    do loop = loop_sta, loop_end, loop_skip
      write( cloop, '(i8.8)' ) loop
      call rb_cnt_gettime( loop, time )

!- |fk| -
      open( offinzvm, file="./data/absfkinzvm_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 
  
        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        call rb_cnt_mxmyisloop( mx, gmy, is, loop, fk )
        if ( trim(flag_normalize) == "phi0" ) fk(:,:,:) = fk(:,:,:) / phi0
        if ( trim(flag_normalize) == "Al0" )  fk(:,:,:) = fk(:,:,:) / (Al0/sqrt(beta))
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            do iz = -global_nz, global_nz-1
              !write( offinzvm, '(g17.7e3)' ) abs(fk(iz,iv,im))
              write( offinzvm, '(g17.7e3)' ) abs(2._DP*pi*gvp(iz,im)*fk(iz,iv,im))
            end do
          end do
        end do

      close( offinzvm )
!-

!- Re[fk] -
      open( offinzvm, file="./data/refkinzvm_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 

        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS Re[fk] float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        call rb_cnt_mxmyisloop( mx, gmy, is, loop, fk )
        if ( trim(flag_normalize) == "phi0" ) fk(:,:,:) = fk(:,:,:) / phi0
        if ( trim(flag_normalize) == "Al0" )  fk(:,:,:) = fk(:,:,:) / (Al0/sqrt(beta))
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            do iz = -global_nz, global_nz-1
              !write( offinzvm, '(g17.7e3)' ) real(fk(iz,iv,im))
              write( offinzvm, '(g17.7e3)' ) real(2._DP*pi*gvp(iz,im)*fk(iz,iv,im))
            end do
          end do
        end do

      close( offinzvm )
!-

!- Im[fk] -
      open( offinzvm, file="./data/imfkinzvm_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 

        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS Im[fk] float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        call rb_cnt_mxmyisloop( mx, gmy, is, loop, fk )
        if ( trim(flag_normalize) == "phi0" ) fk(:,:,:) = fk(:,:,:) / phi0
        if ( trim(flag_normalize) == "Al0" )  fk(:,:,:) = fk(:,:,:) / (Al0/sqrt(beta))
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            do iz = -global_nz, global_nz-1
              !write( offinzvm, '(g17.7e3)' ) aimag(fk(iz,iv,im))
              write( offinzvm, '(g17.7e3)' ) aimag(2._DP*pi*gvp(iz,im)*fk(iz,iv,im))
            end do
          end do
        end do

      close( offinzvm )
!-

    end do

END SUBROUTINE fkinzvm_vtk


SUBROUTINE fkinzvm_connect_vtk( mx, gmy, is, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (zz,vl,mu) in VTK binary format
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyivimisloop,           &
                      rb_phi_gettime, rb_phi_mxmyizloop, loop_phi_end, &
                      rb_Al_gettime, rb_Al_mxmyizloop, loop_Al_end
  use diag_geom, only : kx, gky, n_tht, vmax, beta, gvp, dj, ck

  integer, intent(in) :: mx, gmy, is, loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: fk
  complex(kind=DP) :: wf, phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, nconnect, mxw
  integer :: iz, iv, im, loop

    call rb_phi_gettime( loop_phi_end(enum), time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop_phi_end(enum), phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop_Al_end(enum), Al0 )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cis, fmt="(i1.1)" ) is

    do loop = loop_sta, loop_end, loop_skip
      write( cloop, '(i8.8)' ) loop
      call rb_cnt_gettime( loop, time )

!- |fk| -
      open( offinzvm, file="./data/absfkinzvm_connect_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

      if ( dj(gmy) == 0 ) then

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 
  
        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
            if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
            if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              wf = fk(iz)
              !write( offinzvm, '(g17.7e3)' ) abs(wf)
              write( offinzvm, '(g17.7e3)' ) abs(2._DP*pi*gvp(iz,im)*wf)
            end do
          end do
        end do

      else

        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
        nconnect = connect_max + connect_min + 1

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz*nconnect, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 

        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz*nconnect)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv

            if ( connect_min .ne. 0 ) then
              do iconnect = connect_min, 1, -1
                mxw = mx+iconnect*dj(gmy)
                call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
                if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
                if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
                do iz = -global_nz, global_nz-1
                  wf = ck(gmy)**iconnect * fk(iz)
                  !write( offinzvm, '(g17.7e3)' ) abs(wf)
                  write( offinzvm, '(g17.7e3)' ) abs(2._DP*pi*gvp(iz,im)*wf)
                end do
              end do
            end if
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(gmy)
              call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
              if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
              if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
              do iz = -global_nz, global_nz-1
                wf = conjg(ck(gmy)**iconnect) * fk(iz)
                !write( offinzvm, '(g17.7e3)' ) abs(wf)
                write( offinzvm, '(g17.7e3)' ) abs(2._DP*pi*gvp(iz,im)*wf)
              end do
            end do

          end do
        end do

      end if

      close( offinzvm )
!-

!- Re[fk] -
      open( offinzvm, file="./data/refkinzvm_connect_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

      if ( dj(gmy) == 0 ) then

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 
  
        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
            if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
            if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              wf = fk(iz)
              !write( offinzvm, '(g17.7e3)' ) real(wf)
              write( offinzvm, '(g17.7e3)' ) real(2._DP*pi*gvp(iz,im)*wf)
            end do
          end do
        end do

      else

        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
        nconnect = connect_max + connect_min + 1

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz*nconnect, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 

        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz*nconnect)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv

            if ( connect_min .ne. 0 ) then
              do iconnect = connect_min, 1, -1
                mxw = mx+iconnect*dj(gmy)
                call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
                if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
                if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
                do iz = -global_nz, global_nz-1
                  wf = ck(gmy)**iconnect * fk(iz)
                  !write( offinzvm, '(g17.7e3)' ) real(wf)
                  write( offinzvm, '(g17.7e3)' ) real(2._DP*pi*gvp(iz,im)*wf)
                end do
              end do
            end if
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(gmy)
              call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
              if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
              if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
              do iz = -global_nz, global_nz-1
                wf = conjg(ck(gmy)**iconnect) * fk(iz)
                !write( offinzvm, '(g17.7e3)' ) real(wf)
                write( offinzvm, '(g17.7e3)' ) real(2._DP*pi*gvp(iz,im)*wf)
              end do
            end do

          end do
        end do

      end if

      close( offinzvm )
!-

!- Im[fk] -
      open( offinzvm, file="./data/imfkinzvm_connect_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".vtk",  &
                       status="replace", action="write", form="formatted" )

      if ( dj(gmy) == 0 ) then

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 
  
        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv
            call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
            if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
            if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              wf = fk(iz)
              !write( offinzvm, '(g17.7e3)' ) aimag(wf)
              write( offinzvm, '(g17.7e3)' ) aimag(2._DP*pi*gvp(iz,im)*wf)
            end do
          end do
        end do

      else

        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
        nconnect = connect_max + connect_min + 1

        write( offinzvm, '(a)' ) "# vtk DataFile Version 2.0"
        write( offinzvm, '(a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3,a,g13.3e3)' )  &
                          "time= ", time, " kx= ", kx(mx), " ky= ", gky(gmy),  &
                          " lz= ", pi * real( n_tht, kind=DP ), " vmax= ", vmax
        write( offinzvm, '(a)' ) "ASCII"
        write( offinzvm, '(a)' ) "DATASET STRUCTURED_POINTS"
        write( offinzvm, '(a,3i17)' ) "DIMENSIONS ", 2*global_nz*nconnect, 2*global_nv, global_nm+1
        write( offinzvm, '(A)' ) "ORIGIN 0.0 0.0 0.0"
        write( offinzvm, '(A)' ) "SPACING 1.0 1.0 1.0"
        write( offinzvm, * ) 

        write( offinzvm, '(a,i17)' ) "POINT_DATA ", (2*global_nz*nconnect)*(2*global_nv)*(global_nm+1)
        write( offinzvm, '(a)' ) "SCALARS |fk| float 1"
        write( offinzvm, '(a)' ) "LOOKUP_TABLE default"
        do im = 0, global_nm
          do iv = 1, 2*global_nv

            if ( connect_min .ne. 0 ) then
              do iconnect = connect_min, 1, -1
                mxw = mx+iconnect*dj(gmy)
                call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
                if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
                if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
                do iz = -global_nz, global_nz-1
                  wf = ck(gmy)**iconnect * fk(iz)
                  !write( offinzvm, '(g17.7e3)' ) aimag(wf)
                  write( offinzvm, '(g17.7e3)' ) aimag(2._DP*pi*gvp(iz,im)*wf)
                end do
              end do
            end if
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(gmy)
              call rb_cnt_mxmyivimisloop( mxw, gmy, iv, im, is, loop, fk )
              if ( trim(flag_normalize) == "phi0" ) fk(:) = fk(:) / phi0
              if ( trim(flag_normalize) == "Al0" )  fk(:) = fk(:) / (Al0/sqrt(beta))
              do iz = -global_nz, global_nz-1
                wf = conjg(ck(gmy)**iconnect) * fk(iz)
                !write( offinzvm, '(g17.7e3)' ) aimag(wf)
                write( offinzvm, '(g17.7e3)' ) aimag(2._DP*pi*gvp(iz,im)*wf)
              end do
            end do

          end do
        end do

      end if

      close( offinzvm )
!-

    end do

END SUBROUTINE fkinzvm_connect_vtk


END MODULE out_ffinzvm_vtk
