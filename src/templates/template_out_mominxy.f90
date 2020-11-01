MODULE out_mominxy
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinxy


 CONTAINS


SUBROUTINE phiinxy( giz, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: phikxky
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)   :: phixy
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: mx, my

    call rb_phi_gettime( loop, time )
    call rb_phi_izloop( giz, loop, phikxky )
    call fft_backward_xy( phikxky, phixy )

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/phiinxy_z"//ciz//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy","phi"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my), phixy(mx,my)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE phiinxy


END MODULE out_mominxy
