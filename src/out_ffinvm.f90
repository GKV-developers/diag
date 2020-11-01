MODULE out_ffinvm
!-------------------------------------------------------------------------------
!
!     Output ff in (v,m)
!                                                   (S. Maeyama, 5 May 2022)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fkinvm_fxv, fkinvm_cnt


 CONTAINS


SUBROUTINE fkinvm_fxv( mx, gmy, rankz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_fxv_gettime, rb_fxv_mxmyrankzisloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu, gvp

  integer, intent(in) :: mx, gmy, rankz, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(1:2*global_nv,0:global_nm) :: fk
  character(len=4) :: cmx, cmy, crz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: giz, iv, im

    call rb_fxv_gettime( loop, time )
    call rb_fxv_mxmyrankzisloop( mx, gmy, rankz, is, loop, fk )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( crz, fmt="(i4.4)" ) rankz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    giz = -global_nz + rankz*(2*nz)
    open( offinvm, file="./data/fkinvm_mx"//cmx//"my"//cmy//"rankz"//crz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(99a17)" ) "#              vl","mu","Re[ff]","Im[ff]"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), &
                                          real(fk(iv,im), kind=DP), aimag(fk(iv,im))
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fkinvm_fxv


SUBROUTINE fkinvm_cnt( mx, gmy, giz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyizisloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu, gvp

  integer, intent(in) :: mx, gmy, giz, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(1:2*global_nv,0:global_nm) :: fk
  character(len=4) :: cmx, cmy, ciz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iv, im

    call rb_cnt_gettime( loop, time)
    call rb_cnt_mxmyizisloop( mx, gmy, giz, is, loop, fk )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( ciz, fmt="(i4.4)" ) giz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinvm, file="./data/fkinvm_mx"//cmx//"my"//cmy//"z"//ciz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(99a17)" ) "#              vl","mu","Re[ff]","Im[ff]"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), &
                                          real(fk(iv,im), kind=DP), aimag(fk(iv,im))
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fkinvm_cnt


END MODULE out_ffinvm
