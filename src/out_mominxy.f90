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

  public phiinxy, Alinxy, mominxy,  &
         phiinxy_parity, Alinxy_parity, mominxy_parity


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


SUBROUTINE phiinxy_parity( giz, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, loop

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi,phi_even,phi_odd
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)   :: phixy_even, phixy_odd
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: mx, my

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call parity_mom3(phi, phi_even, phi_odd)
          ! phi_even(-kx,ky,-zz) = phi_even(kx,ky,zz) <Ballooning parity of phi>
          ! phi_odd(-kx,ky,-zz) = - phi_odd(kx,ky,zz) <Tearing parity of phi>
    call fft_backward_xy( phi_even(:,:,giz), phixy_even )
    call fft_backward_xy( phi_odd(:,:,giz), phixy_odd )

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/phiinxy_parity_z"//ciz//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy","phi_even","phi_odd"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my),  &
                                           phixy_even(mx,my), phixy_odd(mx,my)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE phiinxy_parity


SUBROUTINE Alinxy( giz, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_izloop
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: Alkxky
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)   :: Alxy
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: mx, my

    call rb_Al_gettime( loop, time )
    call rb_Al_izloop( giz, loop, Alkxky )
    call fft_backward_xy( Alkxky, Alxy )

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/Alinxy_z"//ciz//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy","Al"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my), Alxy(mx,my)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE Alinxy


SUBROUTINE Alinxy_parity( giz, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, loop

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al,Al_even,Al_odd
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)   :: Alxy_even, Alxy_odd
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: mx, my

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )
    call parity_mom3(Al, Al_even, Al_odd)
          ! Al_even(-kx,ky,-zz) = Al_even(kx,ky,zz) <Tearing parity of Al>
          ! Al_odd(-kx,ky,-zz) = - Al_odd(kx,ky,zz) <Ballooning parity of Al>
    call fft_backward_xy( Al_even(:,:,giz), Alxy_even )
    call fft_backward_xy( Al_odd(:,:,giz), Alxy_odd )

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/Alinxy_parity_z"//ciz//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy","Al_even","Al_odd"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my),  &
                                           Alxy_even(mx,my), Alxy_odd(mx,my)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE Alinxy_parity


SUBROUTINE mominxy( giz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output moments in (xx,yy) at giz, is, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_izimomisloop
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1) :: momkxky
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nmom-1)   :: momxy
  character(len=4) :: ciz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, imom

    call rb_mom_gettime( loop, time )
    do imom = 0, nmom-1
      call rb_mom_izimomisloop( giz, imom, is, loop, momkxky(:,:,imom) )
      call fft_backward_xy( momkxky(:,:,imom), momxy(:,:,imom) )
    end do

    write( ciz, fmt="(i4.4)" ) giz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/mominxy_z"//ciz//"s"//cis//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#   is=",is
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy",&
                                   "dens","upara","ppara",  &
                                   "pperp","qlpara","qlperp"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my), momxy(mx,my,:)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE mominxy


SUBROUTINE mominxy_parity( giz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output moments in (xx,yy) at giz, is, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_isloop
  use diag_parity, only : parity_mom3
  use diag_geom, only : xx, yy, gzz
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, is, loop

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: mom,mom_even,mom_odd
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nmom-1) :: momxy_even, momxy_odd
  character(len=4) :: ciz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, imom

    call rb_mom_gettime( loop, time )
    call rb_mom_isloop( is, loop, mom )
    do imom = 0, nmom-1
      call parity_mom3(mom(:,:,:,imom), mom_even(:,:,:,imom), mom_odd(:,:,:,imom))
      call fft_backward_xy( mom_even(:,:,giz,imom), momxy_even(:,:,imom) )
      call fft_backward_xy( mom_odd(:,:,giz,imom), momxy_odd(:,:,imom) )
    end do
    ! mom_even(-kx,ky,-zz) = mom_even(kx,ky,zz)
    ! mom_odd(-kx,ky,-zz) = - mom_odd(kx,ky,zz)
    ! <Ballooning parity> dens_even, upara_odd, ppara_even, pperp_even, qlpara_odd, qlperp_odd
    !    <Tearing parity> dens_odd, upara_even, ppara_odd, pperp_odd, qlpara_even, qlperp_even

    write( ciz, fmt="(i4.4)" ) giz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/mominxy_parity_z"//ciz//"s"//cis//"_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#   is=",is
      write( omominxy, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( omominxy, "(99a17)" ) "#              xx","yy",&
                                   "dens_even","upara_even","ppara_even",  &
                                   "pperp_even","qlpara_even","qlperp_even", &
                                   "dens_odd","upara_odd","ppara_odd",  &
                                   "pperp_odd","qlpara_odd","qlperp_odd"
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(mx), yy(my), momxy_even(mx,my,:), momxy_odd(mx,my,:)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

END SUBROUTINE mominxy_parity


END MODULE out_mominxy
