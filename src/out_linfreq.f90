MODULE out_linfreq
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public linfreqintime, linfreqinkxky


 CONTAINS


SUBROUTINE linfreqintime( mx, gmy, loop_sta, loop_end, loop_skp )
!-------------------------------------------------------------------------------
!
!     Output omega of a (mx,gmy) mode in time
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop
  use diag_geom, only : kx, gky
  real(kind=DP), parameter :: eps_omega = 1.d-3, &
                              eps_gamma = 1.d-3, &
                              eps_ineq  = 1.d-6

  integer, intent(in) :: mx, gmy, loop_sta, loop_end, loop_skp

  real(kind=DP) :: time, time0
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: phi, phi0
  complex(kind=DP) :: omega, omega0, phi0phi, diff
  real(kind=DP) :: phi_norm2, phi0_norm2, ineq
  logical :: freq_conv
  integer :: iz, loop
  character(len=4) :: cmx, cmy

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy

    open( olinfreq, file="./data/linfreqintime_mx"//cmx//"my"//cmy//".dat" )
      write( olinfreq, "(a,i17,a,g17.7e3)" ) "#   mx=",mx, ",  kx=",kx(mx)
      write( olinfreq, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,",  ky=",gky(gmy)
      write( olinfreq, "(99a17)" ) "#            time","growthrate","frequency",  &
                                   "diff(grow)","diff(freq)","1-Schwartz ineq","conv"

      call rb_phi_gettime( loop_sta, time )
      call rb_phi_mxmyloop( mx, gmy, loop_sta, phi )
      phi_norm2 = 0._DP
      do iz = -global_nz, global_nz-1
        phi_norm2 = phi_norm2 + abs(phi(iz))**2
      end do
      time0 = time
      phi0(:) = phi(:)
      phi0_norm2 = phi_norm2
      omega0 = (0._DP, 0._DP)
  
      do loop = loop_sta+loop_skp, loop_end, loop_skp
  
        call rb_phi_gettime( loop, time )
        call rb_phi_mxmyloop( mx, gmy, loop, phi )
  
        !- calculate interior products -
        phi0phi = (0._DP, 0._DP)
        do iz = -global_nz, global_nz-1
          phi0phi = phi0phi + conjg(phi0(iz)) * phi(iz)
        end do
        phi_norm2 = 0._DP
        do iz = -global_nz, global_nz-1
          phi_norm2 = phi_norm2 + abs(phi(iz))**2
        end do
  
        !- calculate frequency -
        omega = log(phi0phi / phi0_norm2) / (ui * (time0 - time))
  
        !- convergence check -
        diff = abs(real(omega - omega0, kind=DP) / real(omega, kind=DP)) &
             + ui * abs(aimag(omega - omega0) / aimag(omega))
        ineq = abs(phi0phi)**2 / (phi0_norm2 * phi_norm2)
        if ( real(diff, kind=DP) < eps_omega .and. &
                     aimag(diff) < eps_gamma .and. &
                  (1._DP - ineq) < eps_ineq ) then
          freq_conv = .true.
        else
          freq_conv = .false.
        end if
        write( olinfreq, "(6g17.7e3,L17)" ) &
            time, aimag(omega), real(omega, kind=DP), &
            aimag(diff), real(diff, kind=DP), 1._DP-ineq, freq_conv
  
        !- remember the values -
        time0 = time
        phi0(:) = phi(:)
        phi0_norm2 = phi_norm2
        omega0 = omega
  
      end do

    close( olinfreq )


END SUBROUTINE linfreqintime


SUBROUTINE linfreqinkxky( loop_end )
!-------------------------------------------------------------------------------
!
!     Output dispersion relation omega(kx,ky) at loop_end
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : kx, gky
  real(kind=DP), parameter :: eps_omega = 1.d-3, &
                              eps_gamma = 1.d-3, &
                              eps_ineq  = 1.d-6

  integer, intent(in) :: loop_end

  real(kind=DP) :: time, time0
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, phi0
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: omega, omega0, phi0phi, diff, omega_conv
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: phi_norm2, phi0_norm2, ineq
  logical, dimension(-nx:nx,0:global_ny) :: freq_conv
  integer :: mx, my, iz, loop, loop_sta, loop_skp

    loop_sta = loop_end-2
    loop_skp = 1

    call rb_phi_gettime( loop_sta, time )
    call rb_phi_loop( loop_sta, phi )
    phi_norm2(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      phi_norm2(:,:) = phi_norm2(:,:) + abs(phi(:,:,iz))**2
    end do
    time0 = time
    phi0(:,:,:) = phi(:,:,:)
    phi0_norm2(:,:) = phi_norm2(:,:)
    omega0(:,:) = (0._DP, 0._DP)
  
    do loop = loop_sta+loop_skp, loop_end, loop_skp

      call rb_phi_gettime( loop, time )
      call rb_phi_loop( loop, phi )

      !- calculate interior products -
      phi0phi(:,:) = (0._DP, 0._DP)
      do iz = -global_nz, global_nz-1
        phi0phi(:,:) = phi0phi(:,:) + conjg(phi0(:,:,iz)) * phi(:,:,iz)
      end do
      phi_norm2(:,:) = 0._DP
      do iz = -global_nz, global_nz-1
        phi_norm2(:,:) = phi_norm2(:,:) + abs(phi(:,:,iz))**2
      end do

      !- calculate frequency -
      omega(:,:) = log(phi0phi(:,:) / phi0_norm2(:,:)) / (ui * (time0 - time))

      !- convergence check -
      diff(:,:) = abs(real(omega(:,:) - omega0(:,:), kind=DP) / real(omega(:,:), kind=DP)) &
           + ui * abs(aimag(omega(:,:) - omega0(:,:)) / aimag(omega(:,:)))
      ineq(:,:) = abs(phi0phi(:,:))**2 / (phi0_norm2(:,:) * phi_norm2(:,:))
      do my = 0, global_ny
        do mx = -nx, nx
          if ( real(diff(mx,my), kind=DP) < eps_omega .and. &
                       aimag(diff(mx,my)) < eps_gamma .and. &
                    (1._DP - ineq(mx,my)) < eps_ineq ) then
            freq_conv(mx,my) = .true.
            omega_conv(mx,my) = omega(mx,my)
          else
            freq_conv(mx,my) = .false.
            omega_conv(mx,my) = (0._DP, 0._DP)
          end if
        end do
      end do

      !- remember the values -
      time0 = time
      phi0(:,:,:) = phi(:,:,:)
      phi0_norm2(:,:) = phi_norm2(:,:)
      omega0(:,:) = omega(:,:)

    end do

    omega(0,0) = (0._DP, 0._DP)
    diff(0,0) = (0._DP, 0._DP)
    ineq(0,0) = 1._DP

    open( olinfreq, file="./data/linfreqinkxky.dat" )
      write( olinfreq, "(99a17)" ) "#              kx","ky","growth_converged","freq_converged",  &
                                   "grow","freq","diff(grow)","diff(freq)","1-Schwartz ineq","conv"
      do my = 0, global_ny
        do mx = -nx, nx
          write( olinfreq, "(9g17.7e3,L17)" ) &
          kx(mx), gky(my), aimag(omega_conv(mx,my)), real(omega_conv(mx,my), kind=DP), &
          aimag(omega(mx,my)), real(omega(mx,my), kind=DP), &
          aimag(diff(mx,my)), real(diff(mx,my), kind=DP), 1._DP-ineq(mx,my), freq_conv(mx,my)
        end do
        write( olinfreq, * )
      end do
    close( olinfreq )


END SUBROUTINE linfreqinkxky


END MODULE out_linfreq
