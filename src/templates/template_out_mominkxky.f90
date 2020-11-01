MODULE out_mominkxky
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinkxky
         

 CONTAINS


SUBROUTINE phiinkxky( loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (kx,ky) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : kx, gky
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: phi2
  character(len=8) :: cloop
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = 0.5_DP * abs(phi(mx,my,iz))**2
        end do
      end do
    end do
    call intgrl_thet(wr3, phi2)

    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/phiinkxky_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(99a17)" ) "#              kx","ky","<|phi|^2>"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my), phi2(mx,my)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE phiinkxky


END MODULE out_mominkxky
