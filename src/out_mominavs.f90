MODULE out_mominavs
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y,z)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinavs


 CONTAINS


SUBROUTINE phiinavs( flag, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in AVS binary format
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : lx, ly, dz, q_0, s_hat, r_minor, eps_r, ck, dj, n_alp
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: flag ! 0:coord(fluxtube), 1:var&header(fluxtube)
                              ! 2:coord(fulltorus), 3:var&header(fulltorus)
  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: wc2
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wr2
  real(kind=DP) :: wxx, wyy, wzz, wsr, wmr, wdx, wdy, wtheta, wzeta, wq_bar, wx_car, wy_car, wz_car
  integer :: mx, my, iz, loop, mwp, i_alp

    if ( flag == 0 ) then      ! output coordinates

      open( omominxyz, file="./data/phiinavs_tube_coord.dat", status="replace",  &
                       action="write", form="unformatted", access="stream",      &
                       convert="LITTLE_ENDIAN" )

      wdx = lx / real( nxw, kind=DP )
      wdy = ly / real( nyw, kind=DP )
      wq_bar = q_0 * sqrt( 1._DP - eps_r**2 )
      do iz = -global_nz, global_nz
        wzz = dz * real( iz, kind=DP )
        wtheta = 2._DP * atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) * tan(0.5_DP*wzz) )
        !wtheta = wzz
        do my = 0, 2*nyw
          wyy = -ly + wdy * real( my, kind=DP )
          do mx = 0, 2*nxw
            wxx = -lx + wdx * real( mx, kind=DP )
            wsr = r_minor + wxx
            wmr = r_minor / eps_r + wsr * cos( wtheta )
            !wzeta = wq_bar * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )
            wzeta = q_0 * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )
            wx_car = wmr * cos( wzeta )
            wy_car = wmr * sin( wzeta )
            wz_car = wsr * sin( wtheta )
            write( omominxyz ) real(wx_car, kind=4),  &
                               real(wy_car, kind=4),  &
                               real(wz_car, kind=4)
          end do
        end do
      end do

      close( omominxyz )

    else if ( flag == 1 ) then ! output variable

      open( omominxyz, file="./data/phiinavs_tube_var.dat", status="replace",  &
                       action="write", form="unformatted", access="stream",    &
                       convert="LITTLE_ENDIAN" )

      do loop = loop_sta, loop_end, loop_skip
        call rb_phi_gettime( loop, time )
        write( omominxyz ) real(time, kind=4)

        call rb_phi_loop( loop, phi )
        do iz = -global_nz, global_nz-1
          call fft_backward_xy( phi(:,:,iz), wr2 )
          do my = 0, 2*nyw-1
            do mx = 0, 2*nxw-1
              write( omominxyz ) real(wr2(mx,my), kind=4)
            end do
            ! mx = 2*nxw
            write( omominxyz ) real(wr2(0,my), kind=4)
          end do
          ! my = 2*nyw
            do mx = 0, 2*nxw-1
              write( omominxyz ) real(wr2(mx,0), kind=4)
            end do
            ! mx = 2*nxw
            write( omominxyz ) real(wr2(0,0), kind=4)
        end do
        ! iz = global_nz
          do my = 0, global_ny
            do mx = -nx, nx
              mwp = mx - dj(my)          ! --- mw = mx - dj for the positive-z 
                if( abs(mwp) > nx ) then
                  wc2(mx,my) = ( 0._DP, 0._DP )
                else
                  wc2(mx,my) = conjg( ck(my) ) * phi(mwp,my,-global_nz)
                end if
            end do
          end do
          call fft_backward_xy( wc2, wr2 )
          do my = 0, 2*nyw-1
            do mx = 0, 2*nxw-1
              write( omominxyz ) real(wr2(mx,my), kind=4)
            end do
            ! mx = 2*nxw
            write( omominxyz ) real(wr2(0,my), kind=4)
          end do
          ! my = 2*nyw
            do mx = 0, 2*nxw-1
              write( omominxyz ) real(wr2(mx,0), kind=4)
            end do
            ! mx = 2*nxw
            write( omominxyz ) real(wr2(0,0), kind=4)
      end do

      close( omominxyz )

      open( omominxyz, file="./data/phiinavs_tube_header.fld" )

        write( omominxyz, '(a)' ) "# AVS"
        write( omominxyz, '(a,i17)' ) "nstep=", int((loop_end-loop_sta)/loop_skip) + 1
        write( omominxyz, '(a)' ) "ndim=3"
        write( omominxyz, '(a,i17)' ) "dim1=", 2*nxw+1
        write( omominxyz, '(a,i17)' ) "dim2=", 2*nyw+1
        write( omominxyz, '(a,i17)' ) "dim3=", 2*global_nz+1
        write( omominxyz, '(a)' ) "nspace=3"
        write( omominxyz, '(a)' ) "veclen=1"
        write( omominxyz, '(a)' ) "data=float"
        write( omominxyz, '(a)' ) "field=irregular"
        write( omominxyz, '(a)' ) ""
        write( omominxyz, '(a)' ) "time file=phiinavs_tube_var.dat filetype=binary skip=0 close=0"
        write( omominxyz, '(a)' ) "coord 1 file=phiinavs_tube_coord.dat filetype=binary skip=0 stride=3 close=0"
        write( omominxyz, '(a)' ) "coord 2 file=phiinavs_tube_coord.dat filetype=binary skip=4 stride=3 close=0"
        write( omominxyz, '(a)' ) "coord 3 file=phiinavs_tube_coord.dat filetype=binary skip=8 stride=3 close=0"
        write( omominxyz, '(a)' ) "variable 1 file=phiinavs_tube_var.dat filetype=binary skip=4 close=0"
        write( omominxyz, '(a)' ) "EOT"
        write( omominxyz, '(a)' ) ""
        write( omominxyz, '(a)' ) "DO"
        write( omominxyz, '(a,i17,a)' ) "time file=phiinavs_tube_var.dat filetype=binary skip=", &
                                  (2*nxw+1)*(2*nyw+1)*(2*global_nz+1)*4," close=0"
        write( omominxyz, '(a)' ) "variable 1 file=phiinavs_tube_var.dat filetype=binary skip=4 close=0"
        write( omominxyz, '(a)' ) "EOT"
        write( omominxyz, '(a)' ) "ENDDO"

      close( omominxyz )

    else if ( flag == 2 ) then      ! output coordinates

      open( omominxyz, file="./data/phiinavs_full_coord.dat", status="replace",  &
                       action="write", form="unformatted", access="stream",      &
                       convert="LITTLE_ENDIAN" )

      wdx = lx / real( nxw, kind=DP )
      wdy = ly / real( nyw, kind=DP )
      wq_bar = q_0 * sqrt( 1._DP - eps_r**2 )
      do iz = -global_nz, global_nz
        wzz = dz * real( iz, kind=DP )
        wtheta = 2._DP * atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) * tan(0.5_DP*wzz) )
        !wtheta = wzz
        do i_alp = 0, n_alp-1
          do my = 0, 2*nyw-1
            wyy = -ly + wdy * real( my, kind=DP )
            do mx = 0, 2*nxw
              wxx = -lx + wdx * real( mx, kind=DP )
              wsr = r_minor + wxx
              wmr = r_minor / eps_r + wsr * cos( wtheta )
              !wzeta = wq_bar * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )
              !      - 2._DP * pi * real( i_alp, kind=DP ) / real( n_alp, kind=DP )
              wzeta = q_0 * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )  &
                    - 2._DP * pi * real( i_alp, kind=DP ) / real( n_alp, kind=DP )
              wx_car = wmr * cos( wzeta )
              wy_car = wmr * sin( wzeta )
              wz_car = wsr * sin( wtheta )
              write( omominxyz ) real(wx_car, kind=4),  &
                                 real(wy_car, kind=4),  &
                                 real(wz_car, kind=4)
            end do
          end do
        end do
        i_alp = n_alp-1
          my = 2*nyw
            wyy = -ly + wdy * real( my, kind=DP )
            do mx = 0, 2*nxw
              wxx = -lx + wdx * real( mx, kind=DP )
              wsr = r_minor + wxx
              wmr = r_minor / eps_r + wsr * cos( wtheta )
              !wzeta = wq_bar * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )
              !      - 2._DP * pi * real( i_alp, kind=DP ) / real( n_alp, kind=DP )
              wzeta = q_0 * ( wzz + ( s_hat * wxx * wzz - wyy ) / r_minor )  &
                    - 2._DP * pi * real( i_alp, kind=DP ) / real( n_alp, kind=DP )
              wx_car = wmr * cos( wzeta )
              wy_car = wmr * sin( wzeta )
              wz_car = wsr * sin( wtheta )
              write( omominxyz ) real(wx_car, kind=4),  &
                                 real(wy_car, kind=4),  &
                                 real(wz_car, kind=4)
            end do
          
      end do

      close( omominxyz )

    else if ( flag == 3 ) then ! output variable

      open( omominxyz, file="./data/phiinavs_full_var.dat", status="replace",  &
                       action="write", form="unformatted", access="stream",    &
                       convert="LITTLE_ENDIAN" )

      do loop = loop_sta, loop_end, loop_skip
        call rb_phi_gettime( loop, time )
        write( omominxyz ) real(time, kind=4)

        call rb_phi_loop( loop, phi )
        do iz = -global_nz, global_nz-1
          call fft_backward_xy( phi(:,:,iz), wr2 )
          do i_alp = 0, n_alp-1
            do my = 0, 2*nyw-1
              do mx = 0, 2*nxw-1
                write( omominxyz ) real(wr2(mx,my), kind=4)
              end do
              ! mx = 2*nxw
              write( omominxyz ) real(wr2(0,my), kind=4)
            end do
          end do
            ! my = 2*nyw
              do mx = 0, 2*nxw-1
                write( omominxyz ) real(wr2(mx,0), kind=4)
              end do
              ! mx = 2*nxw
              write( omominxyz ) real(wr2(0,0), kind=4)
        end do
        ! iz = global_nz
          do my = 0, global_ny
            do mx = -nx, nx
              mwp = mx - dj(my)          ! --- mw = mx - dj for the positive-z 
                if( abs(mwp) > nx ) then
                  wc2(mx,my) = ( 0._DP, 0._DP )
                else
                  wc2(mx,my) = conjg( ck(my) ) * phi(mwp,my,-global_nz)
                end if
            end do
          end do
          call fft_backward_xy( wc2, wr2 )
          do i_alp = 0, n_alp-1
            do my = 0, 2*nyw-1
              do mx = 0, 2*nxw-1
                write( omominxyz ) real(wr2(mx,my), kind=4)
              end do
              ! mx = 2*nxw
              write( omominxyz ) real(wr2(0,my), kind=4)
            end do
          end do
            ! my = 2*nyw
              do mx = 0, 2*nxw-1
                write( omominxyz ) real(wr2(mx,0), kind=4)
              end do
              ! mx = 2*nxw
              write( omominxyz ) real(wr2(0,0), kind=4)
      end do

      close( omominxyz )

      open( omominxyz, file="./data/phiinavs_full_header.fld" )

        write( omominxyz, '(a)' ) "# AVS"
        write( omominxyz, '(a,i17)' ) "nstep=", int((loop_end-loop_sta)/loop_skip) + 1
        write( omominxyz, '(a)' ) "ndim=3"
        write( omominxyz, '(a,i17)' ) "dim1=", 2*nxw+1
        write( omominxyz, '(a,i17)' ) "dim2=", 2*nyw*n_alp+1
        write( omominxyz, '(a,i17)' ) "dim3=", 2*global_nz+1
        write( omominxyz, '(a)' ) "nspace=3"
        write( omominxyz, '(a)' ) "veclen=1"
        write( omominxyz, '(a)' ) "data=float"
        write( omominxyz, '(a)' ) "field=irregular"
        write( omominxyz, '(a)' ) ""
        write( omominxyz, '(a)' ) "time file=phiinavs_full_var.dat filetype=binary skip=0 close=0"
        write( omominxyz, '(a)' ) "coord 1 file=phiinavs_full_coord.dat filetype=binary skip=0 stride=3 close=0"
        write( omominxyz, '(a)' ) "coord 2 file=phiinavs_full_coord.dat filetype=binary skip=4 stride=3 close=0"
        write( omominxyz, '(a)' ) "coord 3 file=phiinavs_full_coord.dat filetype=binary skip=8 stride=3 close=0"
        write( omominxyz, '(a)' ) "variable 1 file=phiinavs_full_var.dat filetype=binary skip=4 close=0"
        write( omominxyz, '(a)' ) "EOT"
        write( omominxyz, '(a)' ) ""
        write( omominxyz, '(a)' ) "DO"
        write( omominxyz, '(a,i17,a)' ) "time file=phiinavs_full_var.dat filetype=binary skip=", &
                                  (2*nxw+1)*(2*nyw*n_alp+1)*(2*global_nz+1)*4," close=0"
        write( omominxyz, '(a)' ) "variable 1 file=phiinavs_full_var.dat filetype=binary skip=4 close=0"
        write( omominxyz, '(a)' ) "EOT"
        write( omominxyz, '(a)' ) "ENDDO"

      close( omominxyz )

    else

      write(*,*) "flag in phiinavs is invalid !  flag =", flag
      stop

    end if

END SUBROUTINE phiinavs


END MODULE out_mominavs
