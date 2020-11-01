MODULE out_fluidtotaltrans
!-------------------------------------------------------------------------------
!
!     Total triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fluidtotaltrans_izloop, fluidtotaltrans_loop


 CONTAINS


SUBROUTINE fluidtotaltrans_izloop( giz, loop )
!-------------------------------------------------------------------------------
!
!     Total triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop, rb_Al_izloop, rb_mom_izimomisloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gomg, gksq, g0, g1, kx, gky, gzz
  use diag_fft, only : fft_backward_xy, fft_forward_xy

  integer, intent(in) :: giz, loop 

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: phi, Al
  complex(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1,0:ns-1) :: mom
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: ikxf, ikyf, ikxp, ikyp, nf
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: dpdx, dpdy, wkxy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nmom-1,0:ns-1) :: dfdx, dfdy
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1,0:ns-1) :: tk_es, tk_em
  real(kind=DP) :: bb
  integer :: ix, iy, mx, my, imom, is
  character(len=4) :: ciz
  character(len=8) :: cloop, crow

!--- initialize ---
    call rb_phi_gettime( loop, time )
    call rb_phi_izloop( giz, loop, phi )
    call rb_Al_izloop( giz, loop, Al )
    do is = 0, ns-1
      do imom = 0, nmom-1
        call rb_mom_izimomisloop( giz, imom, is, loop, mom(:,:,imom,is) )
      end do
    end do
    !- moments transform: gyrokinetic distribution -> non-adiabatic part
    do is = 0, ns-1
      do my = 0, global_ny 
        do mx = -nx, nx
          bb = gksq(mx,my,giz)*tau(is)*Anum(is)/(Znum(is)**2*gomg(giz)**2)
          mom(mx,my,0,is) = mom(mx,my,0,is)  &
                          + sgn(is) * fcs(is) * g0(mx,my,giz,is) * phi(mx,my) / tau(is)
          mom(mx,my,2,is) = mom(mx,my,2,is)  &
                          + 0.5_DP * sgn(is) * fcs(is) * g0(mx,my,giz,is) * phi(mx,my)
          mom(mx,my,3,is) = mom(mx,my,3,is)  &
                          + sgn(is) * fcs(is) * phi(mx,my)  &
                          * ((1._DP - bb) * g0(mx,my,giz,is) + bb * g1(mx,my,giz,is))
        end do
      end do
    end do
    !--- moments transform: non-adiabatic part -> Hermite-Laguerre coefficients
    do is = 0, ns-1
      do my = 0, global_ny 
        do mx = -nx, nx
          mom(mx,my,0,is) = Znum(is) * mom(mx,my,0,is) / fcs(is)
          mom(mx,my,1,is) = sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,1,is) / fcs(is)
          mom(mx,my,2,is) = 2._DP * Znum(is) * mom(mx,my,2,is) / (fcs(is) * tau(is)) - mom(mx,my,0,is)
          mom(mx,my,3,is) = - Znum(is) * mom(mx,my,3,is) / (fcs(is) * tau(is)) + mom(mx,my,0,is)
          mom(mx,my,4,is) = 2._DP * sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,4,is)  &
                          / (fcs(is) * tau(is)) - 3._DP * mom(mx,my,1,is)
          mom(mx,my,5,is) = - sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,5,is)  &
                          / (fcs(is) * tau(is)) + mom(mx,my,1,is)
        end do
      end do
    end do
!------------------

!--- calc. total transfer ---
    do is = 0, ns-1
      do imom = 0, nmom-1
        do my = 0, global_ny
          do mx = -nx, nx
            ikxf(mx,my) = ui * kx(mx) * mom(mx,my,imom,is)
            ikyf(mx,my) = ui * gky(my) * mom(mx,my,imom,is)
          end do
        end do
        call fft_backward_xy(ikxf, dfdx(:,:,imom,is))
        call fft_backward_xy(ikyf, dfdy(:,:,imom,is))
      end do
    end do

    do my = 0, global_ny
      do mx = -nx, nx
        ikxp(mx,my) = ui * kx(mx) * phi(mx,my)
        ikyp(mx,my) = ui * gky(my) * phi(mx,my)
      end do
    end do
    call fft_backward_xy(ikxp, dpdx)
    call fft_backward_xy(ikyp, dpdy)
    do is = 0, ns-1
      do imom = 0, nmom-1
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,imom,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,imom,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
          do my = 0, global_ny
            do mx = -nx, nx
              tk_es(mx,my,imom,is) = (fcs(is) * tau(is) / Znum(is))  &
                                   * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP )
            end do
          end do
        else if ( imom == 2 ) then
          do my = 0, global_ny
            do mx = -nx, nx
              tk_es(mx,my,imom,is) = (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                   * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP )
            end do
          end do
        else if ( imom == 4 ) then
          do my = 0, global_ny
            do mx = -nx, nx
              tk_es(mx,my,imom,is) = (fcs(is) * tau(is) / Znum(is)) * 0.166666666666666_DP  &
                                   * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP )
            end do
          end do
        else
          write( *, * ) "nmom is wrong.", nmom, imom
          stop
        end if
      end do
    end do

    do my = 0, global_ny
      do mx = -nx, nx
        ikxp(mx,my) = ui * kx(mx) * Al(mx,my)
        ikyp(mx,my) = ui * gky(my) * Al(mx,my)
      end do
    end do
    call fft_backward_xy(ikxp, dpdx)
    call fft_backward_xy(ikyp, dpdy)
    tk_em(:,:,:,:) = 0._DP
    do is = 0, ns-1
      imom = 0
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,0,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,0,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
          end do
        end do
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,1,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,1,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,0,is)) * nf(mx,my), kind=DP )
          end do
        end do
      imom = 1
        !do iy = 0, 2*nyw-1
        !  do ix = 0, 2*nxw-1
        !    wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,1,is)  &
        !                  + dpdy(ix,iy) * dfdx(ix,iy,1,is)
        !  end do
        !end do
        !call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP )
          end do
        end do
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,2,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,2,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
          end do
        end do
      imom = 2
        !do iy = 0, 2*nyw-1
        !  do ix = 0, 2*nxw-1
        !    wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,2,is)  &
        !                  + dpdy(ix,iy) * dfdx(ix,iy,2,is)
        !  end do
        !end do
        !call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,4,is)) * nf(mx,my), kind=DP )
          end do
        end do
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,4,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,4,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP )
          end do
        end do
      imom = 3
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,3,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,3,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,5,is)) * nf(mx,my), kind=DP )
          end do
        end do
        do iy = 0, 2*nyw-1
          do ix = 0, 2*nxw-1
            wkxy(ix,iy) = - dpdx(ix,iy) * dfdy(ix,iy,5,is)  &
                          + dpdy(ix,iy) * dfdx(ix,iy,5,is)
          end do
        end do
        call fft_forward_xy(wkxy, nf)
        do my = 0, global_ny
          do mx = -nx, nx
            tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is)            &
                                 + (fcs(is) * tau(is) / Znum(is))  &
                                 * (- sqrt(tau(is) / Anum(is)))    &
                                 * real( conjg(mom(mx,my,3,is)) * nf(mx,my), kind=DP )
          end do
        end do
    end do
!----------------------------

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    write( crow, '(i0)' ) 2*nx+1
!--- ascii output for gnuplot ---
    open( oflttrninkxky, file="./data/fluidtotaltransinkxky_z"//ciz//"_t"//cloop//".dat" )
      write( oflttrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( oflttrninkxky, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( oflttrninkxky, "(99a17)" ) "#              kx","ky","tk_es","tk_em"
      do my = 0, global_ny
        do mx = -nx, nx
          write( oflttrninkxky, "(99g17.7e3)" ) kx(mx), gky(my), tk_es(mx,my,:,:), tk_em(mx,my,:,:)
        end do
        write( oflttrninkxky, * )
      end do
    close( oflttrninkxky )
!!--- binary output for gnuplot ---
!    !splot './data/fluidtotaltransinkxky_z"//ciz//"_t"//cloop//".dat' "//&
!    !      "binary record=("//trim(crow)//",-1) format='%double%double"//&
!    !      "%double%double%double%double%double%double%double%double%double%double%double%double"//&
!    !      "%double%double%double%double%double%double%double%double%double%double%double%double"//&
!    !      "' u 1:2:3 ti '' w pm3d"
!    open( unit=oflttrninkxky, file="./data/fluidtotaltransinkxky_z"//ciz//"_t"//cloop//".dat", &
!          status="replace", action="write", form="unformatted", access="stream" )
!      do my = 0, global_ny
!        do mx = -nx, nx
!          write( oflttrninkxky ) kx(mx), ky(my), tk_es(mx,my,:,:), tk_em(mx,my,:,:)
!        end do
!      end do
!    close( oflttrninkxky )
!!-----------------------------------

END SUBROUTINE fluidtotaltrans_izloop


SUBROUTINE fluidtotaltrans_loop( loop )
!-------------------------------------------------------------------------------
!
!     Total triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_isloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gomg, gksq, g0, g1, kx, gky, grootg, cfsrf
  use diag_fft, only : fft_backward_xy, fft_forward_xy

  integer, intent(in) :: loop 

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: ikxf, ikyf, ikxp, ikyp, nf
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wkxy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,-global_nz:global_nz-1) :: dpdx, dpdy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: dfdx, dfdy
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1,0:ns-1) :: tk_es, tk_em
  real(kind=DP) :: bb, fct
  integer :: ix, iy, mx, my, iz, imom, is
  character(len=8) :: cloop, crow

!--- initialize ---
    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call rb_Al_loop( loop, Al )
    do is = 0, ns-1
      call rb_mom_isloop( is, loop, mom(:,:,:,:,is) )
    end do
    !- moments transform: gyrokinetic distribution -> non-adiabatic part
    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny 
          do mx = -nx, nx
            bb = gksq(mx,my,iz)*tau(is)*Anum(is)/(Znum(is)**2*gomg(iz)**2)
            mom(mx,my,iz,0,is) = mom(mx,my,iz,0,is)  &
                               + sgn(is) * fcs(is) * g0(mx,my,iz,is) * phi(mx,my,iz) / tau(is)
            mom(mx,my,iz,2,is) = mom(mx,my,iz,2,is)  &
                               + 0.5_DP * sgn(is) * fcs(is) * g0(mx,my,iz,is) * phi(mx,my,iz)
            mom(mx,my,iz,3,is) = mom(mx,my,iz,3,is)  &
                               + sgn(is) * fcs(is) * phi(mx,my,iz)  &
                               * ((1._DP - bb) * g0(mx,my,iz,is) + bb * g1(mx,my,iz,is))
          end do
        end do
      end do
    end do
    !--- moments transform: non-adiabatic part -> Hermite-Laguerre coefficients
    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny 
          do mx = -nx, nx
            mom(mx,my,iz,0,is) = Znum(is) * mom(mx,my,iz,0,is) / fcs(is)
            mom(mx,my,iz,1,is) = sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,iz,1,is) / fcs(is)
            mom(mx,my,iz,2,is) = 2._DP * Znum(is) * mom(mx,my,iz,2,is) / (fcs(is) * tau(is)) - mom(mx,my,iz,0,is)
            mom(mx,my,iz,3,is) = - Znum(is) * mom(mx,my,iz,3,is) / (fcs(is) * tau(is)) + mom(mx,my,iz,0,is)
            mom(mx,my,iz,4,is) = 2._DP * sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,iz,4,is)  &
                               / (fcs(is) * tau(is)) - 3._DP * mom(mx,my,iz,1,is)
            mom(mx,my,iz,5,is) = - sqrt(Anum(is) / tau(is)) * Znum(is) * mom(mx,my,iz,5,is)  &
                               / (fcs(is) * tau(is)) + mom(mx,my,iz,1,is)
          end do
        end do
      end do
    end do
!------------------

!--- calc. total transfer ---
    do is = 0, ns-1
      do imom = 0, nmom-1
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              ikxf(mx,my) = ui * kx(mx) * mom(mx,my,iz,imom,is)
              ikyf(mx,my) = ui * gky(my) * mom(mx,my,iz,imom,is)
            end do
          end do
          call fft_backward_xy(ikxf, dfdx(:,:,iz,imom,is))
          call fft_backward_xy(ikyf, dfdy(:,:,iz,imom,is))
        end do
      end do
    end do

    tk_es(:,:,:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          ikxp(mx,my) = ui * kx(mx) * phi(mx,my,iz)
          ikyp(mx,my) = ui * gky(my) * phi(mx,my,iz)
        end do
      end do
      call fft_backward_xy(ikxp, dpdx(:,:,iz))
      call fft_backward_xy(ikyp, dpdy(:,:,iz))
    end do
    do is = 0, ns-1
      do imom = 0, nmom-1

        if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
          do iz = -global_nz, global_nz-1
            fct = grootg(iz) / cfsrf
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,imom,is)  &
                              + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,imom,is)
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do my = 0, global_ny
              do mx = -nx, nx
                tk_es(mx,my,imom,is) = tk_es(mx,my,imom,is) + fct * (  &! flux-surface average
                                       (fcs(is) * tau(is) / Znum(is))  &
                                     * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) )
              end do
            end do
          end do
        else if ( imom == 2 ) then
          do iz = -global_nz, global_nz-1
            fct = grootg(iz) / cfsrf
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,imom,is)  &
                              + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,imom,is)
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do my = 0, global_ny
              do mx = -nx, nx
                tk_es(mx,my,imom,is) = tk_es(mx,my,imom,is) + fct * (  &! flux-surface average
                                       (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                     * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) )
              end do
            end do
          end do
        else if ( imom == 4 ) then
          do iz = -global_nz, global_nz-1
            fct = grootg(iz) / cfsrf
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,imom,is)  &
                              + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,imom,is)
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do my = 0, global_ny
              do mx = -nx, nx
                tk_es(mx,my,imom,is) = tk_es(mx,my,imom,is) + fct * (  &! flux-surface average
                                       (fcs(is) * tau(is) / Znum(is)) * 0.166666666666666_DP  &
                                     * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) )
              end do
            end do
          end do
        else
          write( *, * ) "nmom is wrong.", nmom, imom
          stop
        end if

      end do
    end do

    tk_em(:,:,:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          ikxp(mx,my) = ui * kx(mx) * Al(mx,my,iz)
          ikyp(mx,my) = ui * gky(my) * Al(mx,my,iz)
        end do
      end do
      call fft_backward_xy(ikxp, dpdx(:,:,iz))
      call fft_backward_xy(ikyp, dpdy(:,:,iz))
    end do
    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / cfsrf

        imom = 0
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,0,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,0,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,1,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,1,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,0,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
        imom = 1
          !do iy = 0, 2*nyw-1
          !  do ix = 0, 2*nxw-1
          !    wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,1,is)  &
          !                  + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,1,is)
          !  end do
          !end do
          !call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,2,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,2,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
        imom = 2
          !do iy = 0, 2*nyw-1
          !  do ix = 0, 2*nxw-1
          !    wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,2,is)  &
          !                  + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,2,is)
          !  end do
          !end do
          !call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,4,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,4,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,4,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is)) * 0.5_DP  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
        imom = 3
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,3,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,3,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,5,is)) * nf(mx,my), kind=DP ) )
            end do
          end do
          do iy = 0, 2*nyw-1
            do ix = 0, 2*nxw-1
              wkxy(ix,iy) = - dpdx(ix,iy,iz) * dfdy(ix,iy,iz,5,is)  &
                            + dpdy(ix,iy,iz) * dfdx(ix,iy,iz,5,is)
            end do
          end do
          call fft_forward_xy(wkxy, nf)
          do my = 0, global_ny
            do mx = -nx, nx
              tk_em(mx,my,imom,is) = tk_em(mx,my,imom,is) + fct * (  &! flux-surface average
                                     (fcs(is) * tau(is) / Znum(is))  &
                                   * (- sqrt(tau(is) / Anum(is)))    &
                                   * real( conjg(mom(mx,my,iz,3,is)) * nf(mx,my), kind=DP ) )
            end do
          end do

      end do
    end do
!----------------------------

    write( cloop, fmt="(i8.8)" ) loop
    write( crow, '(i0)' ) 2*nx+1
!--- ascii output for gnuplot ---
    open( oflttrninkxky, file="./data/fluidtotaltransinkxky_t"//cloop//".dat" )
      write( oflttrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( oflttrninkxky, "(99a17)" ) "#              kx","ky","tk_es","tk_em"
      do my = 0, global_ny
        do mx = -nx, nx
          write( oflttrninkxky, "(99g17.7e3)" ) kx(mx), gky(my), tk_es(mx,my,:,:), tk_em(mx,my,:,:)
        end do
        write( oflttrninkxky, * )
      end do
    close( oflttrninkxky )
!!--- binary output for gnuplot ---
!    !splot './data/fluidtotaltransinkxky_t"//cloop//".dat' "//&
!    !      "binary record=("//trim(crow)//",-1) format='%double%double"//&
!    !      "%double%double%double%double%double%double%double%double%double%double%double%double"//&
!    !      "%double%double%double%double%double%double%double%double%double%double%double%double"//&
!    !      "' u 1:2:3 ti '' w pm3d"
!    open( unit=oflttrninkxky, file="./data/fluidtotaltransinkxky_t"//cloop//".dat", &
!          status="replace", action="write", form="unformatted", access="stream" )
!      do my = 0, global_ny
!        do mx = -nx, nx
!          write( oflttrninkxky ) kx(mx), ky(my), tk_es(mx,my,:,:), tk_em(mx,my,:,:)
!        end do
!      end do
!    close( oflttrninkxky )
!!-----------------------------------

END SUBROUTINE fluidtotaltrans_loop


END MODULE out_fluidtotaltrans
