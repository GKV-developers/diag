MODULE out_fluidsubsptrans
!-------------------------------------------------------------------------------
!
!     Sub-space transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fluidsubsptrans_izloop, fluidsubsptrans_izloop_open,  &
         fluidsubsptrans_izloop_close, fluidsubsptrans_loop,  &
         fluidsubsptrans_loop_open, fluidsubsptrans_loop_close

  integer, parameter :: nfil = 14  ! number of filter

 CONTAINS


SUBROUTINE fluidsubsptrans_izloop( giz, loop )
!-------------------------------------------------------------------------------
!
!     Sub-space transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop, rb_Al_izloop, rb_mom_izimomisloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gomg, gksq, g0, g1, kx, gky
  use diag_fft, only : fft_backward_xy, fft_forward_xy

  integer, intent(in) :: giz, loop 

  real(kind=DP), dimension(-nx:nx,0:global_ny,0:nfil-1) :: filter
  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: phi, Al
  complex(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1,0:ns-1) :: mom
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: ikxf, ikyf, ikxp, ikyp, nf
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wkxy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nfil-1) :: dpdx, dpdy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nfil-1,0:nmom-1,0:ns-1) :: dfdx, dfdy
  real(kind=DP), dimension(0:nfil-1,0:nfil-1,0:nfil-1,0:nmom-1,0:ns-1) :: subtrans_es, subtrans_em
  real(kind=DP), dimension(0:nfil-1,0:nfil-1,0:nfil-1,0:ns-1) :: subtrans
  real(kind=DP) :: bb
  integer :: ix, iy, mx, my, imom, is, kfil, ifil, jfil

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
    !--- set filter
    filter(:,:,:) = 0._DP
    do my = 0, global_ny
      do mx = -nx, nx
      !- medium resolution -
      !  if ( gksq(mx,my,giz) < 1._DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,0) = 1._DP    ! Ion-scale ZF
      !    else
      !      filter(mx,my,1) = 1._DP    ! Ion-scale Turb
      !    end if
      !  else if ( gksq(mx,my,giz) < 4.0_DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,2) = 1._DP    ! Med-scale ZF
      !    else
      !      filter(mx,my,3) = 1._DP    ! Med-scale Turb
      !    end if
      !  else if ( gksq(mx,my,giz) < 16.0_DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,4) = 1._DP    ! Electron-scale ZF
      !    else
      !      filter(mx,my,5) = 1._DP    ! Electron-scale Turb
      !    end if
      !  else
      !    if ( my == 0 ) then
      !      filter(mx,my,6) = 1._DP    ! Ambient-scale ZF
      !    else
      !      filter(mx,my,7) = 1._DP    ! Ambient-scale Turb
      !    end if
      !  end if
      !- fine resolution -
        if ( gksq(mx,my,giz) < 0.5_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,0) = 1._DP    ! Ion-scale ZF
          else
            filter(mx,my,1) = 1._DP    ! Ion-scale Turb
          end if
        else if ( gksq(mx,my,giz) < 1._DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,2) = 1._DP    ! Ion-scale ZF
          else
            filter(mx,my,3) = 1._DP    ! Ion-scale Turb
          end if
        else if ( gksq(mx,my,giz) < 2._DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,4) = 1._DP    ! Med-scale ZF
          else
            filter(mx,my,5) = 1._DP    ! Med-scale Turb
          end if
        else if ( gksq(mx,my,giz) < 4.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,6) = 1._DP    ! Med-scale ZF
          else
            filter(mx,my,7) = 1._DP    ! Med-scale Turb
          end if
        else if ( gksq(mx,my,giz) < 8.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,8) = 1._DP    ! Electron-scale ZF
          else
            filter(mx,my,9) = 1._DP    ! Electron-scale Turb
          end if
        else if ( gksq(mx,my,giz) < 16.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,10) = 1._DP   ! Electron-scale ZF
          else
            filter(mx,my,11) = 1._DP   ! Electron-scale Turb
          end if
        else
          if ( my == 0 ) then
            filter(mx,my,12) = 1._DP   ! Ambient-scale ZF
          else
            filter(mx,my,13) = 1._DP   ! Ambient-scale Turb
          end if
        end if
      end do
    end do
!------------------

!--- calc. subsp transfer ---
    subtrans_es(:,:,:,:,:) = 0._DP
    do ifil = 0, nfil-1
      do my = 0, global_ny
        do mx = -nx, nx
          ikxp(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * phi(mx,my)
          ikyp(mx,my) = filter(mx,my,ifil) * ui * gky(my) * phi(mx,my)
        end do
      end do
      call fft_backward_xy(ikxp, dpdx(:,:,ifil))
      call fft_backward_xy(ikyp, dpdy(:,:,ifil))
    end do
    do is = 0, ns-1
      do imom = 0, nmom-1
        do ifil = 0, nfil-1
          do my = 0, global_ny
            do mx = -nx, nx
              ikxf(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * mom(mx,my,imom,is)
              ikyf(mx,my) = filter(mx,my,ifil) * ui * gky(my) * mom(mx,my,imom,is)
            end do
          end do
          call fft_backward_xy(ikxf, dfdx(:,:,ifil,imom,is))
          call fft_backward_xy(ikyf, dfdy(:,:,ifil,imom,is))
        end do
      end do
    end do
    do is = 0, ns-1
      do imom = 0, nmom-1

        if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP )
                  end do
              end do
            end do
          end do
        else if ( imom == 2 ) then
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                  end do
              end do
            end do
          end do
        else if ( imom == 4 ) then
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP ) * 0.166666666666666_DP
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_es(kfil,ifil,jfil,imom,is) = subtrans_es(kfil,ifil,jfil,imom,is)  &
                                + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                                * real( conjg(mom(mx,my,imom,is)) * nf(mx,my), kind=DP ) * 0.166666666666666_DP
                  end do
              end do
            end do
          end do
        else
          write( *, * ) "nmom is wrong.", nmom, imom
          stop
        end if

      end do
    end do

    subtrans_em(:,:,:,:,:) = 0._DP
    do ifil = 0, nfil-1
      do my = 0, global_ny
        do mx = -nx, nx
          ikxp(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * Al(mx,my)
          ikyp(mx,my) = filter(mx,my,ifil) * ui * gky(my) * Al(mx,my)
        end do
      end do
      call fft_backward_xy(ikxp, dpdx(:,:,ifil))
      call fft_backward_xy(ikyp, dpdy(:,:,ifil))
    end do
    do is = 0, ns-1

      imom = 0
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,0,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,0,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,0,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,0,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
                end do
            end do
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,1,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,1,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,1,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,1,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,0,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,0,is)) * nf(mx,my), kind=DP )
                end do
            end do
          end do
        end do
      imom = 1
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,1,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,1,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,1,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,1,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP )
                end do
            end do
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,2,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,2,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,2,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,2,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,1,is)) * nf(mx,my), kind=DP )
                end do
            end do
          end do
        end do
      imom = 2
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,2,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,2,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,2,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,2,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,4,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,4,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                end do
            end do
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,4,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,4,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,4,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,4,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,2,is)) * nf(mx,my), kind=DP ) * 0.5_DP
                end do
            end do
          end do
        end do
      imom = 3
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,3,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,3,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,3,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,3,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,5,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,5,is)) * nf(mx,my), kind=DP )
                end do
            end do
            do iy = 0, 2*nyw-1
              do ix = 0, 2*nxw-1
                wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,5,is)  &
                                         + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,5,is)  &
                                         - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,5,is)  &
                                         + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,5,is) )
              end do
            end do
            call fft_forward_xy(wkxy, nf)
            do kfil = 0, nfil-1
              do my = 1, global_ny
                do mx = -nx, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,3,is)) * nf(mx,my), kind=DP )
                end do
              end do
              my = 0
                do mx = 1, nx
                  subtrans_em(kfil,ifil,jfil,imom,is) = subtrans_em(kfil,ifil,jfil,imom,is)  &
                              + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                              * (- sqrt(tau(is) / Anum(is)))                                 &
                              * real( conjg(mom(mx,my,3,is)) * nf(mx,my), kind=DP )
                end do
            end do
          end do
        end do

    end do

    subtrans(:,:,:,:) = 0._DP
    do is = 0, ns-1
      do imom = 0, nmom-1
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do kfil = 0, nfil-1
              subtrans(kfil,ifil,jfil,is) = subtrans(kfil,ifil,jfil,is)          &
                                          + subtrans_es(kfil,ifil,jfil,imom,is)  &
                                          + subtrans_em(kfil,ifil,jfil,imom,is)
            end do
          end do
        end do
      end do
    end do
!----------------------------

!--- ascii output for gnuplot ---
    do is = 0, ns-1
      do kfil = 0, nfil-1
        write( oflstrninkxky+kfil+is*5000000,  &
               "(g17.7e3)", advance="NO" ) time
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            if ( ifil == jfil ) then
              write( oflstrninkxky+kfil+is*5000000,  &
                     "(g17.7e3)", advance="NO" ) subtrans(kfil,ifil,jfil,is)
            else
              write( oflstrninkxky+kfil+is*5000000,  &
                     "(g17.7e3)", advance="NO" ) 2._DP * subtrans(kfil,ifil,jfil,is)
            end if
          end do
        end do
        write( oflstrninkxky+kfil+is*5000000, * )
      end do
    end do
!--------------------------------

    is = 0
    write(99,*) time, subtrans(1,2,3,is), subtrans(2,3,1,is), subtrans(3,1,2,is),  &
                      subtrans(2,1,3,is), subtrans(1,3,2,is), subtrans(3,2,1,is)
    write(98,*) time, subtrans(7,10,13,is), subtrans(10,13,7,is), subtrans(13,7,10,is),  &
                      subtrans(10,7,13,is), subtrans(7,13,10,is), subtrans(13,10,7,is)
    is = 1
    write(97,*) time, subtrans(1,2,3,is), subtrans(2,3,1,is), subtrans(3,1,2,is),  &
                      subtrans(2,1,3,is), subtrans(1,3,2,is), subtrans(3,2,1,is)
    write(96,*) time, subtrans(7,10,13,is), subtrans(10,13,7,is), subtrans(13,7,10,is),  &
                      subtrans(10,7,13,is), subtrans(7,13,10,is), subtrans(13,10,7,is)

END SUBROUTINE fluidsubsptrans_izloop


SUBROUTINE fluidsubsptrans_izloop_open( giz )
!-------------------------------------------------------------------------------
!
!     file open
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz

  character(len=4) :: ciz, cfil
  character(len=1) :: cis
  integer :: is, kfil

    write( ciz, "(i4.4)" ) giz
    do is = 0, ns-1
      write( cis, "(i1.1)" ) is
      do kfil = 0, nfil-1
        write( cfil, "(i4.4)" ) kfil 
        open( oflstrninkxky+kfil+is*5000000,  &
              file="./data/fluidsubsptrans_z"//ciz//"kfil"//cfil//"s"//cis//".dat" )
      end do
    end do

END SUBROUTINE fluidsubsptrans_izloop_open


SUBROUTINE fluidsubsptrans_izloop_close
!-------------------------------------------------------------------------------
!
!     file close
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer :: is, kfil

    do is = 0, ns-1
      do kfil = 0, nfil-1
        close( oflstrninkxky+kfil+is*5000000 )
      end do
    end do

END SUBROUTINE fluidsubsptrans_izloop_close


SUBROUTINE fluidsubsptrans_loop( loop )
!-------------------------------------------------------------------------------
!
!     Sub-space transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_isloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gomg, gksq, g0, g1, kx, gky, grootg, cfsrf
  use diag_fft, only : fft_backward_xy, fft_forward_xy

  integer, intent(in) :: loop 

  real(kind=DP), dimension(-nx:nx,0:global_ny,0:nfil-1) :: filter
  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: ikxf, ikyf, ikxp, ikyp, nf
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wkxy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nfil-1) :: dpdx, dpdy
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:nfil-1,0:nmom-1,0:ns-1) :: dfdx, dfdy
  real(kind=DP), dimension(0:nfil-1,0:nfil-1,0:nfil-1,0:nmom-1,0:ns-1) :: subtrans_es, subtrans_em
  real(kind=DP), dimension(0:nfil-1,0:nfil-1,0:nfil-1,0:ns-1) :: subtrans
  real(kind=DP) :: bb, fct
  integer :: ix, iy, mx, my, iz, imom, is, kfil, ifil, jfil

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
    !--- set filter
    filter(:,:,:) = 0._DP
    iz = 0
    do my = 0, global_ny
      do mx = -nx, nx
      !- medium resolution -
      !  if ( gksq(mx,my,iz) < 1._DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,0) = 1._DP    ! Ion-scale ZF
      !    else
      !      filter(mx,my,1) = 1._DP    ! Ion-scale Turb
      !    end if
      !  else if ( gksq(mx,my,iz) < 4.0_DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,2) = 1._DP    ! Med-scale ZF
      !    else
      !      filter(mx,my,3) = 1._DP    ! Med-scale Turb
      !    end if
      !  else if ( gksq(mx,my,iz) < 16.0_DP**2 ) then
      !    if ( my == 0 ) then
      !      filter(mx,my,4) = 1._DP    ! Electron-scale ZF
      !    else
      !      filter(mx,my,5) = 1._DP    ! Electron-scale Turb
      !    end if
      !  else
      !    if ( my == 0 ) then
      !      filter(mx,my,6) = 1._DP    ! Ambient-scale ZF
      !    else
      !      filter(mx,my,7) = 1._DP    ! Ambient-scale Turb
      !    end if
      !  end if
      !- fine resolution -
        if ( gksq(mx,my,iz) < 0.5_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,0) = 1._DP    ! Ion-scale ZF
          else
            filter(mx,my,1) = 1._DP    ! Ion-scale Turb
          end if
        else if ( gksq(mx,my,iz) < 1._DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,2) = 1._DP    ! Ion-scale ZF
          else
            filter(mx,my,3) = 1._DP    ! Ion-scale Turb
          end if
        else if ( gksq(mx,my,iz) < 2._DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,4) = 1._DP    ! Med-scale ZF
          else
            filter(mx,my,5) = 1._DP    ! Med-scale Turb
          end if
        else if ( gksq(mx,my,iz) < 4.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,6) = 1._DP    ! Med-scale ZF
          else
            filter(mx,my,7) = 1._DP    ! Med-scale Turb
          end if
        else if ( gksq(mx,my,iz) < 8.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,8) = 1._DP    ! Electron-scale ZF
          else
            filter(mx,my,9) = 1._DP    ! Electron-scale Turb
          end if
        else if ( gksq(mx,my,iz) < 16.0_DP**2 ) then
          if ( my == 0 ) then
            filter(mx,my,10) = 1._DP   ! Electron-scale ZF
          else
            filter(mx,my,11) = 1._DP   ! Electron-scale Turb
          end if
        else
          if ( my == 0 ) then
            filter(mx,my,12) = 1._DP   ! Ambient-scale ZF
          else
            filter(mx,my,13) = 1._DP   ! Ambient-scale Turb
          end if
        end if
      end do
    end do
!------------------

!--- calc. subsp transfer ---
    subtrans_es(:,:,:,:,:) = 0._DP
    subtrans_em(:,:,:,:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf

      do ifil = 0, nfil-1
        do my = 0, global_ny
          do mx = -nx, nx
            ikxp(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * phi(mx,my,iz)
            ikyp(mx,my) = filter(mx,my,ifil) * ui * gky(my) * phi(mx,my,iz)
          end do
        end do
        call fft_backward_xy(ikxp, dpdx(:,:,ifil))
        call fft_backward_xy(ikyp, dpdy(:,:,ifil))
      end do
      do is = 0, ns-1
        do imom = 0, nmom-1
          do ifil = 0, nfil-1
            do my = 0, global_ny
              do mx = -nx, nx
                ikxf(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * mom(mx,my,iz,imom,is)
                ikyf(mx,my) = filter(mx,my,ifil) * ui * gky(my) * mom(mx,my,iz,imom,is)
              end do
            end do
            call fft_backward_xy(ikxf, dfdx(:,:,ifil,imom,is))
            call fft_backward_xy(ikyf, dfdy(:,:,ifil,imom,is))
          end do
        end do
      end do
      do is = 0, ns-1
        do imom = 0, nmom-1
  
          if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
            do jfil = 0, nfil-1
              do ifil = 0, nfil-1
                do iy = 0, 2*nyw-1
                  do ix = 0, 2*nxw-1
                    wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                             + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                             - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                             + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                  end do
                end do
                call fft_forward_xy(wkxy, nf)
                do kfil = 0, nfil-1
                  do my = 1, global_ny
                    do mx = -nx, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) )
                    end do
                  end do
                  my = 0
                    do mx = 1, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) )
                    end do
                end do
              end do
            end do
          else if ( imom == 2 ) then
            do jfil = 0, nfil-1
              do ifil = 0, nfil-1
                do iy = 0, 2*nyw-1
                  do ix = 0, 2*nxw-1
                    wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                             + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                             - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                             + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                  end do
                end do
                call fft_forward_xy(wkxy, nf)
                do kfil = 0, nfil-1
                  do my = 1, global_ny
                    do mx = -nx, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                    end do
                  end do
                  my = 0
                    do mx = 1, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                    end do
                end do
              end do
            end do
          else if ( imom == 4 ) then
            do jfil = 0, nfil-1
              do ifil = 0, nfil-1
                do iy = 0, 2*nyw-1
                  do ix = 0, 2*nxw-1
                    wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,imom,is)  &
                                             + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,imom,is)  &
                                             - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,imom,is)  &
                                             + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,imom,is) )
                  end do
                end do
                call fft_forward_xy(wkxy, nf)
                do kfil = 0, nfil-1
                  do my = 1, global_ny
                    do mx = -nx, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) * 0.166666666666666_DP )
                    end do
                  end do
                  my = 0
                    do mx = 1, nx
                      subtrans_es(kfil,ifil,jfil,imom,is) =                              &
                          subtrans_es(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                          + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                          * real( conjg(mom(mx,my,iz,imom,is)) * nf(mx,my), kind=DP ) * 0.166666666666666_DP )
                    end do
                end do
              end do
            end do
          else
            write( *, * ) "nmom is wrong.", nmom, imom
            stop
          end if
  
        end do
      end do


      do ifil = 0, nfil-1
        do my = 0, global_ny
          do mx = -nx, nx
            ikxp(mx,my) = filter(mx,my,ifil) * ui * kx(mx) * Al(mx,my,iz)
            ikyp(mx,my) = filter(mx,my,ifil) * ui * gky(my) * Al(mx,my,iz)
          end do
        end do
        call fft_backward_xy(ikxp, dpdx(:,:,ifil))
        call fft_backward_xy(ikyp, dpdy(:,:,ifil))
      end do
      do is = 0, ns-1
        imom = 0
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,0,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,0,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,0,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,0,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,1,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,1,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,1,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,1,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,0,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,0,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
            end do
          end do
        imom = 1
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,1,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,1,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,1,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,1,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,2,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,2,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,2,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,2,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,1,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
            end do
          end do
        imom = 2
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,2,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,2,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,2,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,2,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,4,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,4,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                  end do
              end do
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,4,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,4,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,4,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,4,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,2,is)) * nf(mx,my), kind=DP ) * 0.5_DP )
                  end do
              end do
            end do
          end do
        imom = 3
          do jfil = 0, nfil-1
            do ifil = 0, nfil-1
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,3,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,3,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,3,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,3,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,5,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,5,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
              do iy = 0, 2*nyw-1
                do ix = 0, 2*nxw-1
                  wkxy(ix,iy) = 0.5_DP * ( - dpdx(ix,iy,ifil) * dfdy(ix,iy,jfil,5,is)  &
                                           + dpdy(ix,iy,ifil) * dfdx(ix,iy,jfil,5,is)  &
                                           - dpdx(ix,iy,jfil) * dfdy(ix,iy,ifil,5,is)  &
                                           + dpdy(ix,iy,jfil) * dfdx(ix,iy,ifil,5,is) )
                end do
              end do
              call fft_forward_xy(wkxy, nf)
              do kfil = 0, nfil-1
                do my = 1, global_ny
                  do mx = -nx, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,3,is)) * nf(mx,my), kind=DP ) )
                  end do
                end do
                my = 0
                  do mx = 1, nx
                    subtrans_em(kfil,ifil,jfil,imom,is) =                              &
                        subtrans_em(kfil,ifil,jfil,imom,is) + fct * ( & ! flux-surface average
                        + 2._DP * filter(mx,my,kfil) * (fcs(is) * tau(is) / Znum(is))  &
                        * (- sqrt(tau(is) / Anum(is)))                                 &
                        * real( conjg(mom(mx,my,iz,3,is)) * nf(mx,my), kind=DP ) )
                  end do
              end do
            end do
          end do
  
      end do

    end do

    subtrans(:,:,:,:) = 0._DP
    do is = 0, ns-1
      do imom = 0, nmom-1
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            do kfil = 0, nfil-1
              subtrans(kfil,ifil,jfil,is) = subtrans(kfil,ifil,jfil,is)          &
                                          + subtrans_es(kfil,ifil,jfil,imom,is)  &
                                          + subtrans_em(kfil,ifil,jfil,imom,is)
            end do
          end do
        end do
      end do
    end do
!----------------------------

!--- ascii output for gnuplot ---
    do is = 0, ns-1
      do kfil = 0, nfil-1
        write( oflstrninkxky+kfil+is*5000000,  &
               "(g17.7e3)", advance="NO" ) time
        do jfil = 0, nfil-1
          do ifil = 0, nfil-1
            if ( ifil == jfil ) then
              write( oflstrninkxky+kfil+is*5000000,  &
                     "(g17.7e3)", advance="NO" ) subtrans(kfil,ifil,jfil,is)
            else
              write( oflstrninkxky+kfil+is*5000000,  &
                     "(g17.7e3)", advance="NO" ) 2._DP * subtrans(kfil,ifil,jfil,is)
            end if
          end do
        end do
        write( oflstrninkxky+kfil+is*5000000, * )
      end do
    end do
!--------------------------------

    !is = 0
    !write(99,*) time, subtrans(1,2,3,is), subtrans(2,3,1,is), subtrans(3,1,2,is),  &
    !                  subtrans(2,1,3,is), subtrans(1,3,2,is), subtrans(3,2,1,is)
    !write(98,*) time, subtrans(7,10,13,is), subtrans(10,13,7,is), subtrans(13,7,10,is),  &
    !                  subtrans(10,7,13,is), subtrans(7,13,10,is), subtrans(13,10,7,is)
    !is = 1
    !write(97,*) time, subtrans(1,2,3,is), subtrans(2,3,1,is), subtrans(3,1,2,is),  &
    !                  subtrans(2,1,3,is), subtrans(1,3,2,is), subtrans(3,2,1,is)
    !write(96,*) time, subtrans(7,10,13,is), subtrans(10,13,7,is), subtrans(13,7,10,is),  &
    !                  subtrans(10,7,13,is), subtrans(7,13,10,is), subtrans(13,10,7,is)

END SUBROUTINE fluidsubsptrans_loop


SUBROUTINE fluidsubsptrans_loop_open
!-------------------------------------------------------------------------------
!
!     file open
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------

  character(len=4) :: cfil
  character(len=1) :: cis
  integer :: is, kfil

    do is = 0, ns-1
      write( cis, "(i1.1)" ) is
      do kfil = 0, nfil-1
        write( cfil, "(i4.4)" ) kfil 
        open( oflstrninkxky+kfil+is*5000000,  &
              file="./data/fluidsubsptrans_kfil"//cfil//"s"//cis//".dat" )
      end do
    end do

END SUBROUTINE fluidsubsptrans_loop_open


SUBROUTINE fluidsubsptrans_loop_close
!-------------------------------------------------------------------------------
!
!     file close
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer :: is, kfil

    do is = 0, ns-1
      do kfil = 0, nfil-1
        close( oflstrninkxky+kfil+is*5000000 )
      end do
    end do

END SUBROUTINE fluidsubsptrans_loop_close


END MODULE out_fluidsubsptrans
