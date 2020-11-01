MODULE out_ffinkxky
!-------------------------------------------------------------------------------
!
!     Output ff in (v,m)
!                                                   (S. Maeyama, 5 May 2022)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fkinkxky_fxv, fkinkxky_cnt


 CONTAINS


SUBROUTINE fkinkxky_fxv( rankz, giv, gim, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk in (kx,ky) at a given (zz,vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_fxv_gettime, rb_fxv_rankzivimisloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu

  integer, intent(in) :: rankz, giv, gim, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: fk
  character(len=4) :: crz, civ, cim
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: giz, mx, my

    call rb_fxv_gettime( loop, time )
    call rb_fxv_rankzivimisloop( rankz, giv, gim, is, loop, fk )

    write( crz, fmt="(i4.4)" ) rankz
    write( civ, fmt="(i4.4)" ) giv
    write( cim, fmt="(i4.4)" ) gim
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    giz = -global_nz + rankz*(2*nz)
    open( offinkxky, file="./data/fkinkxky_rankz"//crz//"v"//civ//"m"//cim//"s"//cis//"_t"//cloop//".dat" )
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  giv=",giv,  ",  gvl=",gvl(giv)
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  gim=",gim,  ",  gmu=",gmu(gim)
      write( offinkxky, "(99a17)" ) "#              kx","ky","|ff|^2"
      do my = 0, global_ny
        do mx = -nx, nx
          write( offinkxky, '(99g17.7e3)' ) kx(mx), gky(my), abs(fk(mx,my))**2
        end do
        write( offinkxky, * )
      end do
    close( offinkxky )

END SUBROUTINE fkinkxky_fxv


SUBROUTINE fkinkxky_cnt( giz, giv, gim, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk in (kx,ky) at a given (zz,vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_ivimisloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu

  integer, intent(in) :: giz, giv, gim, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: fk
  character(len=4) :: ciz, civ, cim
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my

    call rb_cnt_gettime( loop, time )
    call rb_cnt_ivimisloop( giv, gim, is, loop, fk )

    write( ciz, fmt="(i4.4)" ) giz
    write( civ, fmt="(i4.4)" ) giv
    write( cim, fmt="(i4.4)" ) gim
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinkxky, file="./data/fkinkxky_z"//ciz//"v"//civ//"m"//cim//"s"//cis//"_t"//cloop//".dat" )
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  giv=",giv,  ",  gvl=",gvl(giv)
      write( offinkxky, "(a,i17,a,g17.7e3)" ) "#  gim=",gim,  ",  gmu=",gmu(gim)
      write( offinkxky, "(99a17)" ) "#              kx","ky","|ff|^2"
      do my = 0, global_ny
        do mx = -nx, nx
          write( offinkxky, '(99g17.7e3)' ) kx(mx), gky(my), abs(fk(mx,my,giz))**2
        end do
        write( offinkxky, * )
      end do
    close( offinkxky )

END SUBROUTINE fkinkxky_cnt


END MODULE out_ffinkxky
