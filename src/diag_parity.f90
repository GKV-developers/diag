MODULE diag_parity
!-------------------------------------------------------------------------------
!
!     Module for integral in phase space
!                                                   (S. Maeyama, 16 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public parity_mom3, parity_mom2


 CONTAINS


SUBROUTINE parity_mom3(ww, weven, wodd)
!-------------------------------------------------------------------------------
!
!     weven(kx,ky,zz) = (ww(kx,ky,zz) + ww(-kx,ky,-zz) * exp(-2*i*kx*xs)) / 2
!      wodd(kx,ky,zz) = (ww(kx,ky,zz) - ww(-kx,ky,-zz) * exp(-2*i*kx*xs)) / 2
!
!     where xs = cy*q0/s_hat = 2*pi*del_c/m_j is the mode rational surface
!     in a flux-tube geometry.
!                                                   (S. Maeyama, 16 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : dj, ck, del_c, s_hat, kymin, kx
  implicit none
  complex(kind=DP), intent(in), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) ::  ww
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) ::  weven, wodd
  complex(kind=DP), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz) ::  wk
  real(kind=DP) :: xs
  complex(kind=DP) :: wp, wn, cp
  integer :: mx, my, iz, mwp

  !%%% Copy to wk %%%
    wk(-nx:nx,0:global_ny,-global_nz:global_nz-1) = ww(-nx:nx,0:global_ny,-global_nz:global_nz-1)
    do my = 0, global_ny
      if ( dj(my) == 0 ) then
        wk(-nx:nx,my,global_nz) = wk(-nx:nx,my,-global_nz)
      else
        do mx = -nx, nx
          mwp = mx - dj(my) ! --- mw = mx - dj for the positive-z 

          if( abs(mwp) > nx ) then
            !wk(mx,my,global_nz) = 0._DP ! zero-fixed boundary for diagnostics
            wk(mx,my,global_nz) = wk(mx,my,global_nz-1) ! Extrapolation for diagnostics
          else
            wk(mx,my,global_nz) = conjg(ck(my)) * wk(mwp,my,-global_nz)
          end if
        end do
      end if
    end do

    xs = - del_c / (s_hat * kymin)

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wp = wk(mx,my,iz)
          wn = wk(-mx,my,-iz)
          cp = exp(- ui * 2._DP * kx(mx) * xs)
          weven(mx,my,iz) = 0.5_DP * (wp + wn * cp)
           wodd(mx,my,iz) = 0.5_DP * (wp - wn * cp)
        end do
      end do
    end do
  

END SUBROUTINE parity_mom3


SUBROUTINE parity_mom2(gmy, ww, weven, wodd)
!-------------------------------------------------------------------------------
!
!     For a given ky,
!
!     weven(kx,zz) = (ww(kx,zz) + ww(-kx,-zz) * exp(-2*i*kx*xs)) / 2
!      wodd(kx,zz) = (ww(kx,zz) - ww(-kx,-zz) * exp(-2*i*kx*xs)) / 2
!
!     where xs = cy*q0/s_hat = 2*pi*del_c/m_j is the mode rational surface
!     in a flux-tube geometry.
!                                                   (S. Maeyama, 16 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : dj, ck, del_c, s_hat, kymin, kx
  implicit none
  integer :: gmy
  complex(kind=DP), intent(in), &
    dimension(-nx:nx,-global_nz:global_nz-1) ::  ww
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,-global_nz:global_nz-1) ::  weven, wodd
  complex(kind=DP), &
    dimension(-nx:nx,-global_nz:global_nz) ::  wk
  real(kind=DP) :: xs
  complex(kind=DP) :: wp, wn, cp
  integer :: mx, iz, mwp

  !%%% Copy to wk %%%
    wk(-nx:nx,-global_nz:global_nz-1) = ww(-nx:nx,-global_nz:global_nz-1)
    if ( dj(gmy) == 0 ) then
      wk(-nx:nx,global_nz) = wk(-nx:nx,-global_nz)
    else
      do mx = -nx, nx
        mwp = mx - dj(gmy) ! --- mw = mx - dj for the positive-z 

        if( abs(mwp) > nx ) then
          !wk(mx,global_nz) = 0._DP ! zero-fixed boundary for diagnostics
          wk(mx,global_nz) = wk(mx,global_nz-1) ! Extrapolation for diagnostics
        else
          wk(mx,global_nz) = conjg(ck(gmy)) * wk(mwp,-global_nz)
        end if
      end do
    end if

    xs = - del_c / (s_hat * kymin)

    do iz = -global_nz, global_nz-1
      do mx = -nx, nx
        wp = wk(mx,iz)
        wn = wk(-mx,-iz)
        cp = exp(- ui * 2._DP * kx(mx) * xs)
        weven(mx,iz) = 0.5_DP * (wp + wn * cp)
         wodd(mx,iz) = 0.5_DP * (wp - wn * cp)
      end do
    end do
  

END SUBROUTINE parity_mom2


END MODULE diag_parity
