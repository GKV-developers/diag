MODULE out_mominrz
!-------------------------------------------------------------------------------
!
!     Output moments in (R,Z)
!                                                   (S. Maeyama, 29 Oct. 2020)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinrz

  integer, parameter :: n_alp = 3
  ! Adapt a flux tube as 1/n_alp torus for visualization.
  ! Then, Larmor radius rho/r_major = pi*eps_r/(q_0*ly*n_alp)

  integer, parameter :: increase_nzw = 1
  ! Default poloidal grid number is, nzw_default = int(nyw*q_0*n_alp)
  ! Linear interpolation increases resolution by nzw = increase_nzw * nzw_default

 CONTAINS

SUBROUTINE phiinrz( loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 29 Oct. 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : xx, gky, gzz, ck, dj, eps_r, q_0, s_hat, lx, ly, &
                        lz, dz, n_tht
  use diag_fft, only : fft_backward_x1

  integer, intent(in) :: loop

  real(kind=DP) :: time
  integer :: nzw
  real(kind=DP) :: zeta
  complex(kind=DP), dimension(:,:,:), allocatable :: phikxkyz, phixkyz
  complex(kind=DP), dimension(:,:), allocatable :: phi_interp
  real(kind=DP), dimension(:,:), allocatable :: mr, z_car, phi_pol
  complex(kind=DP), dimension(0:2*nxw-1) :: phix
  real(kind=DP) :: alpha, wdx, wdz, wxx, wyy, wzz, rho, q_r, wtheta
  character(len=8) :: cloop
  integer :: ix, mx, my, iz, izw, mwp

    rho = pi * eps_r / (q_0 * ly * n_alp) ! = Larmor radius rho/r_major
    !write(*,*) "# Larmor radius rho/r_major = ", rho
    if (lx*rho > eps_r) then
      write(*,*) "# WARNING in out_mominvtk. lx*rho < eps_r is recommended. Set larger n_alp."
      write(*,*) "# lx=",lx,", rho=",rho,", eps_r=",eps_r,", n_alp=",n_alp 
    end if

    nzw = increase_nzw * int(nyw*n_alp*q_0)

    allocate(phikxkyz(-nx:nx,0:global_ny,-global_nz:global_nz))
    allocate(phixkyz(0:2*nxw,0:global_ny,-global_nz:global_nz))
    allocate(phi_interp(0:2*nxw,0:global_ny))
    allocate(mr(0:2*nxw,-nzw:nzw))
    allocate(z_car(0:2*nxw,-nzw:nzw))
    allocate(phi_pol(0:2*nxw,-nzw:nzw))

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phikxkyz(:,:,-global_nz:global_nz-1) )

    !- pseudo-periodic boundary in z -
    iz = global_nz
      do my = 0, global_ny
        do mx = -nx, nx
          mwp = mx - dj(my)          ! --- mw = mx - dj for the positive-z 
            if( abs(mwp) > nx ) then
              phikxkyz(mx,my,iz) = ( 0._DP, 0._DP )
            else
              phikxkyz(mx,my,iz) = conjg( ck(my) ) * phikxkyz(mwp,my,-global_nz)
            end if
        end do
      end do

    !- fft in x -
!$OMP parallel do default(none) &
!$OMP shared(phikxkyz,phixkyz) &
!$OMP private(my,iz,phix)
    do iz = -global_nz, global_nz
      do my = 0, global_ny
        call fft_backward_x1(phikxkyz(-nx:nx,my,iz), phix)
        phixkyz(0:2*nxw-1,my,iz) = phix(0:2*nxw-1)
        phixkyz(2*nxw,my,iz) = phix(0)
      end do
    end do

    zeta = 0._DP
    wdx = lx / real( nxw, kind=DP )
    !wdz = lz / real( nzw, kind=DP )
    wdz = lz / real( n_tht*nzw, kind=DP ) ! Modify for n_tht>1
!$OMP parallel do default(none) &
!$OMP shared(nzw,zeta,dz,gzz,gky,wdx,wdz,lx,ly,q_0,eps_r,s_hat,rho,phixkyz,phi_pol,n_tht) &
!$OMP private(ix,my,iz,izw,wzz,phi_interp,alpha,wtheta,wxx,q_r,wyy)
    do izw = -nzw, nzw

      !- interpolate along z -
      wzz = wdz * real(izw, kind=DP)
      if (izw == -nzw) then
        !phi_interp(:,:) = phixkyz(:,:,-global_nz)
        phi_interp(:,:) = phixkyz(:,:,-global_nz/n_tht) ! Modify for n_tht>1
      else if (izw == nzw) then
        !phi_interp(:,:) = phixkyz(:,:,global_nz)
        phi_interp(:,:) = phixkyz(:,:,global_nz/n_tht) ! Modify for n_tht>1
      else
        if (wzz > 0) then
          iz = int(wzz / dz) ! find position zz(iz) <= wzz < zz(iz+1)
        else
          iz = int(wzz / dz) -1 ! find position zz(iz) <= wzz < zz(iz+1)
        end if
        alpha = (wzz - gzz(iz)) / dz
        !write(*,*) gzz(iz), wzz, gzz(iz+1), alpha ! for debug
        phi_interp(:,:) = (1._DP - alpha) * phixkyz(:,:,iz) + alpha * phixkyz(:,:,iz+1)
      end if

      !- inverse fft in y at zeta=0 -
      wtheta = wzz
      do ix = 0, 2*nxw
        wxx = -lx + wdx * real( ix, kind=DP )
        q_r = q_0 * ( 1._DP + s_hat * wxx * rho / eps_r )
        wyy = eps_r * (q_r*wtheta - zeta) / (q_0 * rho)
        wyy = wyy + ly ! since -ly<=yy<ly in GKV, rather than 0<=yy<2*ly
        phi_pol(ix,izw) = real(phi_interp(ix,0), kind=DP)
        do my = 1, global_ny
          phi_pol(ix,izw) = phi_pol(ix,izw) + 2._DP &
                          * real(exp(ui * gky(my) * wyy) * phi_interp(ix,my), kind=DP)
        end do
      end do

    end do

    !- prepare structured grid in (R,Z) -
    call cartesian_coordinates(nzw, mr, z_car)

    write( cloop, fmt="(i8.8)" ) loop
    open( omominrz, file="./data/phiinrz_t"//cloop//".dat" )
      write( omominrz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominrz, "(a,i17,a,g17.7e3)" ) "#n_alp=",n_alp,", Larmor radius rho/r_major=",rho
      write( omominrz, "(a,i17,a,g17.7e3)" ) "#  nzw=",nzw,  ", zeta=",zeta
      write( omominrz, "(99a17)" ) "#         major R","height Z","phi"
      do iz = -nzw, nzw
        do ix = 0, 2*nxw
          write( omominrz, "(99g17.7e3)" ) mr(ix,iz), z_car(ix,iz), phi_pol(ix,iz)
        end do
        write( omominrz, * )
      end do
    close( omominrz )

    deallocate(phikxkyz)
    deallocate(phixkyz)
    deallocate(phi_interp)
    deallocate(mr)
    deallocate(z_car)
    deallocate(phi_pol)

END SUBROUTINE phiinrz


SUBROUTINE cartesian_coordinates( nzw, mr, z_car )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : q_0, s_hat, eps_r
  implicit none
  integer, intent(in) :: nzw
  real(kind=DP), intent(out), &
    dimension(0:2*nxw,-nzw:nzw) :: mr, z_car

    call cartesian_coordinates_salpha(mr,z_car,nzw,q_0,s_hat,eps_r)
    !call cartesian_coordinates_miller(mr,z_car,nzw,q_0,s_hat,eps_r,&
    !       dRmildr=-0.1d0,dZmildr=0.d0,kappa=1.5d0,s_kappa=0.7d0,delta=0.4d0,s_delta=1.3d0,zetasq=0.d0,s_zetasq=0.d0)

END SUBROUTINE cartesian_coordinates


SUBROUTINE cartesian_coordinates_salpha( mr, z_car, nzw, q_0, s_hat, eps_r )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz, n_tht
  implicit none
  real(kind=DP), intent(out), &
    dimension(0:2*nxw,-nzw:nzw) :: mr, z_car
  integer, intent(in) :: nzw
  real(kind=DP), intent(in) :: q_0, s_hat, eps_r
  real(kind=DP) :: rho, q_r
  real(kind=DP) :: wxx, wzz, wsr, wmr, wdx, wdz, wtheta, wz_car
  integer :: mx, iz

    rho = pi * eps_r / (q_0 * ly * n_alp) ! = Larmor radius rho/r_major

    wdx = lx / real( nxw, kind=DP )
    !wdz = lz / real( nzw, kind=DP )
    wdz = lz / real( n_tht*nzw, kind=DP ) ! Modify for n_tht>1
    do iz = -nzw, nzw
      wzz = wdz * real( iz, kind=DP )
      wtheta = wzz
      do mx = 0, 2*nxw
        wxx = -lx + wdx * real( mx, kind=DP )
        q_r = q_0 * ( 1._DP + s_hat * wxx * rho / eps_r )

        wsr = eps_r + rho * wxx
        wmr = 1._DP + wsr * cos( wtheta )
        wz_car = wsr * sin( wtheta )

        mr(mx,iz) = wmr
        z_car(mx,iz) = wz_car
      end do
    end do

END SUBROUTINE cartesian_coordinates_salpha


SUBROUTINE cartesian_coordinates_miller( mr, z_car, nzw, q_0, s_hat, eps_r, &
             dRmildr, dZmildr, kappa, s_kappa, delta, s_delta, zetasq, s_zetasq )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz, n_tht
  implicit none
  real(kind=DP), intent(out), &
    dimension(0:2*nxw,-nzw:nzw) :: mr, z_car
  integer, intent(in) :: nzw
  real(kind=DP), intent(in) :: q_0, s_hat, eps_r, dRmildr, dZmildr, kappa, s_kappa, delta, s_delta, zetasq, s_zetasq
  real(kind=DP) :: rho, q_r, kappa_r, delta_r, zetasq_r, Rmil_r, Zmil_r
  real(kind=DP) :: wxx, wzz, wsr, wmr, wdx, wdz, wtheta, wz_car
  integer :: mx, iz

    rho = pi * eps_r / (q_0 * ly * n_alp) ! = Larmor radius rho/r_major

    wdx = lx / real( nxw, kind=DP )
    !wdz = lz / real( nzw, kind=DP )
    wdz = lz / real( n_tht*nzw, kind=DP ) ! Modify for n_tht>1
    do iz = -nzw, nzw
      wzz = wdz * real( iz, kind=DP )
      wtheta = wzz
      do mx = 0, 2*nxw
        wxx = -lx + wdx * real( mx, kind=DP )
        q_r = q_0 * ( 1._DP + s_hat * wxx * rho / eps_r )
        kappa_r = kappa * ( 1._DP + s_kappa * wxx * rho / eps_r )
        delta_r = delta + sqrt(1._DP - delta**2) * s_delta * wxx * rho / eps_r
        zetasq_r = zetasq + s_zetasq * wxx * rho / eps_r
        Rmil_r = 1.d0 + dRmildr * wxx * rho
        Zmil_r = 0.d0 + dZmildr * wxx * rho

        wsr = eps_r + rho * wxx
        wmr = Rmil_r + wsr * cos( wtheta + asin(delta_r) * sin(wtheta) )
        wz_car = Zmil_r + wsr * kappa_r * sin( wtheta + zetasq_r * sin(2._DP*wtheta) )

        mr(mx,iz) = wmr
        z_car(mx,iz) = wz_car
      end do
    end do

END SUBROUTINE cartesian_coordinates_miller


END MODULE out_mominrz
