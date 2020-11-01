MODULE out_textfile
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiintext


 CONTAINS


SUBROUTINE phiintext
!-------------------------------------------------------------------------------
!
!     Output phi(kx,ky,zz,time) in textfile
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, loop_phi_sta, loop_phi_end
  use diag_geom, only : kx, gky, gzz

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  character(len=3) :: cnum
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: mx, my, iz, loop, inum

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      open( otextfile, file="./data/phi."//cnum )
      do loop = loop_phi_sta(inum), loop_phi_end(inum)
        call rb_phi_gettime( loop, time )
        call rb_phi_loop( loop, phi )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              write( otextfile, "(99e17.7e3)" ) kx(mx), gky(my), gzz(iz), time, phi(mx,my,iz)
            end do
          end do
        end do
      end do
      close( otextfile )
    end do

END SUBROUTINE phiintext


END MODULE out_textfile
