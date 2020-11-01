MODULE diag_intgrl
!-------------------------------------------------------------------------------
!
!     Module for integral in phase space
!                                                   (S. Maeyama, 8 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public intgrl_thet

  INTERFACE intgrl_thet
    module procedure intgrl_thet_r, intgrl_thet_z
  END INTERFACE intgrl_thet


 CONTAINS


SUBROUTINE intgrl_thet_z(wc3, wc2)
!-------------------------------------------------------------------------------
!
!     Field-alinged (zz) average of a complex wc3(kx,ky,zz) -> wc2(kx,ky)
!                                                   (S. Maeyama, 8 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : grootg, cfsrf
  implicit none
  complex(kind=DP), intent(in), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wc3
  complex(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny) :: wc2
  real(kind=DP) :: fct
  integer :: mx, my, iz
 
!$OMP parallel workshare
    wc2(:,:) = (0._DP, 0._DP)
!$OMP end parallel workshare

!$OMP parallel default(none)       &
!$OMP shared(wc2,wc3,grootg,cfsrf) &
!$OMP private(mx,my,iz,fct)
!$OMP do reduction(+:wc2)
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf != sqrt(g) dz / (int sqrt(g) dz)
      do my = 0, global_ny
        do mx = -nx, nx
          wc2(mx,my) = wc2(mx,my) + fct * wc3(mx,my,iz)
        end do
      end do
    end do
!$OMP end do
!$OMP end parallel

END SUBROUTINE intgrl_thet_z


SUBROUTINE intgrl_thet_r(wr3, wr2)
!-------------------------------------------------------------------------------
!
!     Field-alinged (zz) average of a real wr3(kx,ky,zz) -> wr2(kx,ky)
!                                                   (S. Maeyama, 8 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : grootg, cfsrf
  implicit none
  real(kind=DP), intent(in), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), intent(out), &
    dimension(-nx:nx,0:global_ny) :: wr2
  real(kind=DP) :: fct
  integer :: mx, my, iz
 
!$OMP parallel workshare
    wr2(:,:) = 0._DP
!$OMP end parallel workshare

!$OMP parallel default(none)       &
!$OMP shared(wr2,wr3,grootg,cfsrf) &
!$OMP private(mx,my,iz,fct)
!$OMP do reduction(+:wr2)
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf != sqrt(g) dz / (int sqrt(g) dz)
      do my = 0, global_ny
        do mx = -nx, nx
          wr2(mx,my) = wr2(mx,my) + fct * wr3(mx,my,iz)
        end do
      end do
    end do
!$OMP end do
!$OMP end parallel

END SUBROUTINE intgrl_thet_r

END MODULE diag_intgrl
