MODULE out_mominxz
!-------------------------------------------------------------------------------
!
!     Output moments in (x,z)
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinxz, Alinxz, mominxz


 CONTAINS


SUBROUTINE phiinxz( gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,zz) at gmy, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_myloop, rb_Al_mxmyizloop
  use diag_geom, only : xx, gky, gzz, kx, r_minor, s_hat, lx
  use diag_fft, only : fft_backward_x1

  integer, intent(in) :: gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) :: phikxz
  complex(kind=DP), dimension(0:2*nxw-1,-global_nz:global_nz-1) :: phixz
  complex(kind=DP) :: Al0
  character(len=4) :: cmy
  character(len=8) :: cloop
  integer :: mx, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_myloop( gmy, loop, phikxz )
    call rb_Al_mxmyizloop( 0, gmy, 0, loop, Al0 )
    do iz = -global_nz, global_nz-1 ! Rewind phase with ck
      do mx = -nx, nx
        phikxz(mx,iz) = exp(ui*kx(mx)*lx) * exp(ui*kx(mx)*(-r_minor/s_hat))  &
                      * phikxz(mx,iz)
                      !* phikxz(mx,iz) / Al0
      end do
    end do
    do iz = -global_nz, global_nz-1
      call fft_backward_x1( phikxz(-nx:nx,iz), phixz(0:2*nxw-1,iz) )
    end do

    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxz, file="./data/phiinxz_my"//cmy//"_t"//cloop//".dat" )
      write( omominxz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominxz, "(99a17)" ) "#              xx","zz","phi"
      do iz = -global_nz, global_nz-1
        do mx = 0, 2*nxw-1
          write( omominxz, "(99g17.7e3)" ) xx(mx), gzz(iz),    &
            real(phixz(mx,iz), kind=DP), aimag(phixz(mx,iz))
        end do
        write( omominxz, * )
      end do
    close( omominxz )

END SUBROUTINE phiinxz


SUBROUTINE Alinxz( gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (xx,zz) at gmy, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_myloop, rb_Al_mxmyizloop
  use diag_geom, only : xx, gky, gzz, kx, r_minor, s_hat, lx
  use diag_fft, only : fft_backward_x1

  integer, intent(in) :: gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) :: Alkxz, bxkxz, bykxz
  complex(kind=DP), dimension(0:2*nxw-1,-global_nz:global_nz-1) :: Alxz, bxxz, byxz
  complex(kind=DP) :: Al0
  character(len=4) :: cmy
  character(len=8) :: cloop
  integer :: mx, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_myloop( gmy, loop, Alkxz )
    call rb_Al_mxmyizloop( 0, gmy, 0, loop, Al0 )
    do iz = -global_nz, global_nz-1 ! Rewind phase with ck
      do mx = -nx, nx
        Alkxz(mx,iz) = exp(ui*kx(mx)*lx) * exp(ui*kx(mx)*(-r_minor/s_hat))  &
                     * Alkxz(mx,iz)
                     !* Alkxz(mx,iz) / Al0
        bxkxz(mx,iz) = ui * gky(gmy) * Alkxz(mx,iz)
        bykxz(mx,iz) = - ui * kx(mx) * Alkxz(mx,iz)
      end do
    end do
    do iz = -global_nz, global_nz-1
      call fft_backward_x1( Alkxz(-nx:nx,iz), Alxz(0:2*nxw-1,iz) )
      call fft_backward_x1( bxkxz(-nx:nx,iz), bxxz(0:2*nxw-1,iz) )
      call fft_backward_x1( bykxz(-nx:nx,iz), byxz(0:2*nxw-1,iz) )
    end do

    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxz, file="./data/Alinxz_my"//cmy//"_t"//cloop//".dat" )
      write( omominxz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominxz, "(99a17)" ) "#              xx","zz","Al"
      do iz = -global_nz, global_nz-1
        do mx = 0, 2*nxw-1
          write( omominxz, "(99g17.7e3)" ) xx(mx), gzz(iz),    &
            real(Alxz(mx,iz), kind=DP), aimag(Alxz(mx,iz)),    &
            real(bxxz(mx,iz), kind=DP), aimag(bxxz(mx,iz)),    &
            real(byxz(mx,iz), kind=DP), aimag(byxz(mx,iz))
        end do
        write( omominxz, * )
      end do
    close( omominxz )

END SUBROUTINE Alinxz


SUBROUTINE mominxz( gmy, imom, is, loop )
!-------------------------------------------------------------------------------
!
!     Output mom in (xx,zz) at gmy, imom, is, loop
!                                                   (S. Maeyama, 11 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_myimomisloop, rb_Al_mxmyizloop
  use diag_geom, only : xx, gky, gzz, kx, r_minor, s_hat, lx
  use diag_fft, only : fft_backward_x1

  integer, intent(in) :: gmy, imom, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) :: momkxz
  complex(kind=DP), dimension(0:2*nxw-1,-global_nz:global_nz-1) :: momxz
  complex(kind=DP) :: Al0
  character(len=4) :: cmy
  character(len=1) :: cimom, cis
  character(len=8) :: cloop
  integer :: mx, iz

    call rb_mom_gettime( loop, time )
    call rb_mom_myimomisloop( gmy, imom, is, loop, momkxz )
    call rb_Al_mxmyizloop( 0, gmy, 0, loop, Al0 )
    do iz = -global_nz, global_nz-1 ! Rewind phase with ck
      do mx = -nx, nx
        momkxz(mx,iz) = exp(ui*kx(mx)*lx) * exp(ui*kx(mx)*(-r_minor/s_hat))  &
                      * momkxz(mx,iz)
                      !* momkxz(mx,iz) / Al0
      end do
    end do
    do iz = -global_nz, global_nz-1
      call fft_backward_x1( momkxz(-nx:nx,iz), momxz(0:2*nxw-1,iz) )
    end do

    write( cmy, fmt="(i4.4)" ) gmy
    write( cimom, fmt="(i1.1)" ) imom
    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxz, file="./data/mominxz_my"//cmy//"mom"//cimom//"s"//cis//"_t"//cloop//".dat" )
      write( omominxz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxz, "(a,i17,a,g17.7e3)" ) "#   is=",is,   ", imom=",imom
      write( omominxz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominxz, "(99a17)" ) "#              xx","zz","mom"
      do iz = -global_nz, global_nz-1
        do mx = 0, 2*nxw-1
          write( omominxz, "(99g17.7e3)" ) xx(mx), gzz(iz),    &
            real(momxz(mx,iz), kind=DP), aimag(momxz(mx,iz)),  &
            atan2(aimag(momxz(mx,iz)),dble(momxz(mx,iz)))
        end do
        write( omominxz, * )
      end do
    close( omominxz )

END SUBROUTINE mominxz


END MODULE out_mominxz
