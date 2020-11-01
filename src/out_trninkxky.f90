MODULE out_trninkxky
!-------------------------------------------------------------------------------
!
!     Output transfer diagnostics in (kx,ky)

!       trn(:,:, 0) :: Entropy S_s
!       trn(:,:, 1) :: Electrostatic field energy W_E
!       trn(:,:, 2) :: Magnetic field energy W_M
!       trn(:,:, 3) :: W_E to S_s interaction R_sE
!       trn(:,:, 4) :: W_M to S_s interaction R_sM
!       trn(:,:, 5) :: Entropy transfer via ExB nonlinearity I_sE
!       trn(:,:, 6) :: Entropy transfer via magnetic nonlinearity I_sM
!       trn(:,:, 7) :: Collisional dissipation D_s
!       trn(:,:, 8) :: Particle flux by ExB flows G_sE
!       trn(:,:, 9) :: Particle flux by magnetic flutters G_sM
!       trn(:,:,10) :: Energy flux by ExB flows Q_sE
!       trn(:,:,11) :: Energy flux by magnetic flutters Q_sM
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public trninkxky


 CONTAINS


SUBROUTINE trninkxky( is, loop )
!-------------------------------------------------------------------------------
!
!     Output transfer diagnostics in (kx,ky) at is, loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_trn_gettime, rb_trn_isloop
  use diag_geom, only : kx, gky

  integer, intent(in) :: is, loop

  real(kind=DP) :: time
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:ntrn-1) :: trn
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my

    call rb_trn_gettime( loop, time )
    call rb_trn_isloop( is, loop, trn )

    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( otrninkxky, file="./data/trninkxky_s"//cis//"_t"//cloop//".dat" )
      write( otrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( otrninkxky, "(a,i17)" ) "#   is=",is
      write( otrninkxky, "(99a17)" ) "#              kx","ky","S_s",         &
                                     "W_E","W_M","R_sE","R_sM","I_sE","I_sM",&
                                     "D_s","G_sE","G_sM","Q_sE","Q_sM"
      do my = 0, global_ny
        do mx = -nx, nx
          write( otrninkxky, "(99g17.7e3)" ) kx(mx), gky(my), trn(mx,my,:)
        end do
        write( otrninkxky, * )
      end do
    close( otrninkxky )

END SUBROUTINE trninkxky


SUBROUTINE trninkxky_all( loop )
!-------------------------------------------------------------------------------
!
!     Output transfer diagnostics in (kx,ky) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_trn_gettime, rb_trn_loop
  use diag_geom, only : kx, gky

  integer, intent(in) :: loop

  real(kind=DP) :: time
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:ntrn-1,0:ns-1) :: trn
  character(len=8) :: cloop
  integer :: mx, my

    call rb_trn_gettime( loop, time )
    call rb_trn_loop( loop, trn )

    write( cloop, fmt="(i8.8)" ) loop
    open( otrninkxky, file="./data/trninkxky_t"//cloop//".dat" )
      write( otrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( otrninkxky, "(99a17)" ) "#              kx","ky","S_s",         &
                                     "W_E","W_M","R_sE","R_sM","I_sE","I_sM",&
                                     "D_s","G_sE","G_sM","Q_sE","Q_sM"
      do my = 0, global_ny
        do mx = -nx, nx
          write( otrninkxky, "(99g17.7e3)" ) kx(mx), gky(my), trn(mx,my,:,:)
        end do
        write( otrninkxky, * )
      end do
    close( otrninkxky )


                                     !%%% For debug %%%
                                     !mx=0
                                     !my=2
                                     !is=0
                                     !write(999,*) time, &
                                     !trn(mx,my,3,is), trn(mx,my,4,is)
                                     !%%%%%%%%%%%%%%%%%

END SUBROUTINE trninkxky_all


END MODULE out_trninkxky
