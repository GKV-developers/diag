MODULE out_triinkxky
!-------------------------------------------------------------------------------
!
!     Output transfer diagnostics in (kx,ky)
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public triinkxky


 CONTAINS


SUBROUTINE triinkxky( mxt, myt, is, loop )
!-------------------------------------------------------------------------------
!
!     Output transfer diagnostics in (kx,ky) of mxt, myt, is at loop
!                                                   (S. Maeyama, 6 Nov 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_tri_gettime, rb_tri_mxtmytisloop
  use diag_geom, only : kx, gky

  integer, intent(in) :: mxt, myt, is, loop

  real(kind=DP) :: time
  real(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:ntri-1) :: tri
  character(len=1) :: cis
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: px, py

    call rb_tri_gettime( loop, time )
    call rb_tri_mxtmytisloop( mxt, myt, is, loop, tri )

    write( cis, fmt="(i1.1)" ) is
    write( cmx, fmt="(i4.4)" ) mxt
    write( cmy, fmt="(i4.4)" ) myt
    write( cloop, fmt="(i8.8)" ) loop
    open( otriinkxky, file="./data/triinkxky_s"//cis//"mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( otriinkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( otriinkxky, "(a,i17)" )           "#   is=",is
      write( otriinkxky, "(a,i17,a,g17.7e3)" ) "#  mxt=",mxt,  ",  kx=",kx(mxt)
      write( otriinkxky, "(a,i17,a,g17.7e3)" ) "#  myt=",myt,  ",  ky=",gky(myt)
      write( otriinkxky, "(99a17)" ) "#              px","py","jkpq_es","jpqk_es","jqkp_es",&
                                                              "jkpq_em","jpqk_em","jqkp_em"
      do py = -global_ny, -1
        do px = -nx, nx
          write( otriinkxky, "(99g17.7e3)" ) kx(px), -gky(-py), tri(px,py,:)
        end do
        write( otriinkxky, * )
      end do
      do py = 0, global_ny
        do px = -nx, nx
          write( otriinkxky, "(99g17.7e3)" ) kx(px), gky(py), tri(px,py,:)
        end do
        write( otriinkxky, * )
      end do
    close( otriinkxky )

                                     !%%% For debug %%%
                                     !write(9000+100*mxt+10*myt+is,*) &
                                     !time, sum(tri(:,:,0)), sum(tri(:,:,3))
                                     !%%%%%%%%%%%%%%%%%

END SUBROUTINE triinkxky


END MODULE out_triinkxky
