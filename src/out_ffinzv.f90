MODULE out_ffinzv
!-------------------------------------------------------------------------------
!
!     Output ff in (z,v)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fkinzv


 CONTAINS


SUBROUTINE fkinzv( mx, gmy, gim, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (zz,vl)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyimisloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu

  integer, intent(in) :: mx, gmy, gim, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv) :: fk
  character(len=4) :: cmx, cmy, cim
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iz, iv

    call rb_cnt_gettime( loop, time )
    call rb_cnt_mxmyimisloop( mx, gmy, gim, is, loop, fk )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cim, fmt="(i4.4)" ) gim
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinzv, file="./data/fkinzv_mx"//cmx//"my"//cmy//"m"//cim//"s"//cis//"_t"//cloop//".dat" )
      write( offinzv, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinzv, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( offinzv, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( offinzv, "(a,i17,a,g17.7e3)" ) "#  gim=",gim,  ",  gmu=",gmu(gim)
      write( offinzv, "(99a17)" ) "#              zz","vl","Re[ff]","Im[ff]"
      do iv = 1, 2*global_nv
        do iz = -global_nz, global_nz-1
          write( offinzv, '(99g17.7e3)' ) gzz(iz), gvl(iv),  &
                                          real(fk(iz,iv), kind=DP), aimag(fk(iz,iv))
        end do
        write( offinzv, * )
      end do
    close( offinzv )

END SUBROUTINE fkinzv


END MODULE out_ffinzv
