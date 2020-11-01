MODULE out_fluiddetailtrans
!-------------------------------------------------------------------------------
!
!     Detailed triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fluiddetailtrans_izloop, fluiddetailtrans_izloop_open,  &
         fluiddetailtrans_izloop_close, fluiddetailtrans_loop,   &
         fluiddetailtrans_loop_open, fluiddetailtrans_loop_close

  real(kind=DP), dimension(-global_ny:global_ny) :: wky

 CONTAINS


SUBROUTINE fluiddetailtrans_izloop( giz, loop )
!-------------------------------------------------------------------------------
!
!     Detailed triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop, rb_Al_izloop, rb_mom_izimomisloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gomg, gksq, g0, g1, gky

  integer, intent(in) :: giz, loop 

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_ny:global_ny) :: phi, Al
  complex(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:nmom-1,0:ns-1) :: mom
  real(kind=DP) :: bb
  integer :: mx, my, imom, is

!--- initialize ---
    call rb_phi_gettime( loop, time )
    call rb_phi_izloop( giz, loop, phi(:,0:global_ny) )
    call rb_Al_izloop( giz, loop, Al(:,0:global_ny) )
    do is = 0, ns-1
      do imom = 0, nmom-1
        call rb_mom_izimomisloop( giz, imom, is, loop, mom(:,0:global_ny,imom,is) )
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
    do my = 1, global_ny
      do mx = -nx, nx
        phi(-mx,-my) = conjg( phi(mx,my) )
        Al(-mx,-my) = conjg( Al(mx,my) )
      end do
    end do
    do is = 0, ns-1
      do imom = 0, nmom-1
        do my = 1, global_ny
          do mx = -nx, nx
            mom(-mx,-my,imom,is) = conjg( mom(mx,my,imom,is) )
          end do
        end do
      end do
    end do
    do my = 0, global_ny 
      wky(-my) = - gky(my)
      wky(my) = gky(my)
    end do
!------------------

    call fluiddetailtrans_izloop_calc( loop,  0, 150, giz, time, phi, Al, mom ) ! ETG
    call fluiddetailtrans_izloop_calc( loop,  0,  66, giz, time, phi, Al, mom ) ! Streamer
    call fluiddetailtrans_izloop_calc( loop,  0,   5, giz, time, phi, Al, mom ) ! ITG
    call fluiddetailtrans_izloop_calc( loop,  0,   3, giz, time, phi, Al, mom ) ! Turb
    call fluiddetailtrans_izloop_calc( loop, 66,   0, giz, time, phi, Al, mom ) ! High ZF
    call fluiddetailtrans_izloop_calc( loop, 24,   0, giz, time, phi, Al, mom ) ! Med ZF
    call fluiddetailtrans_izloop_calc( loop,  1,   0, giz, time, phi, Al, mom ) ! Low ZF

END SUBROUTINE fluiddetailtrans_izloop


SUBROUTINE fluiddetailtrans_izloop_calc( loop, diag_mx, diag_my, giz, &
                                         time, phi, Al, mom )
!-------------------------------------------------------------------------------
!
!     Detailed triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, kx, gky, gzz

  integer, intent(in) :: loop, diag_mx, diag_my, giz

  real(kind=DP), intent(in) :: time
  complex(kind=DP), intent(in), dimension(-nx:nx,-global_ny:global_ny) :: phi, Al
  complex(kind=DP), intent(in), dimension(-nx:nx,-global_ny:global_ny,0:nmom-1,0:ns-1) :: mom

  real(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:nmom-1,0:ns-1) ::  &
                             jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em
  real(kind=DP), dimension(0:nmom-1,0:ns-1) :: tk_es, tk_es_pos, tk_es_neg,  &
                                               tk_em, tk_em_pos, tk_em_neg
  real(kind=DP) :: ceff
  character(len=4) :: cmx, cmy, ciz
  character(len=8) :: cloop, crow
  integer :: mx, my, px, py, qx, qy, imom, is

!--- calc. detail transfer ---
    jkpq_es(:,:,:,:) = 0._DP
    jpqk_es(:,:,:,:) = 0._DP
    jqkp_es(:,:,:,:) = 0._DP
    mx = diag_mx
    my = diag_my
    do is = 0, ns-1
      do imom = 0, nmom-1

        if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is)
        else if ( imom == 2 ) then
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * 0.5_DP
        else if ( imom == 4 ) then
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * 0.166666666666666_DP
        else
          write( *, * ) "nmom is wrong.", nmom, imom
          stop
        end if
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          qy = - py - my
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            qx = - px - mx
            jkpq_es(px,py,imom,is) = - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                                       * real( (   phi(px,py) * mom(qx,qy,imom,is)     &
                                                 - phi(qx,qy) * mom(px,py,imom,is) )   &
                                               * mom(mx,my,imom,is), kind=DP )
            jpqk_es(px,py,imom,is) = - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                                       * real( (   phi(qx,qy) * mom(mx,my,imom,is)     &
                                                 - phi(mx,my) * mom(qx,qy,imom,is) )   &
                                               * mom(px,py,imom,is), kind=DP )
            jqkp_es(px,py,imom,is) = - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                                       * real( (   phi(mx,my) * mom(px,py,imom,is)     &
                                                 - phi(px,py) * mom(mx,my,imom,is) )   &
                                               * mom(qx,qy,imom,is), kind=DP )
          end do
        end do

      end do
    end do

    jkpq_em(:,:,:,:) = 0._DP
    jpqk_em(:,:,:,:) = 0._DP
    jqkp_em(:,:,:,:) = 0._DP
    mx = diag_mx
    my = diag_my
    do is = 0, ns-1

      imom = 0
        ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          qy = - py - my
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            qx = - px - mx
            jkpq_em(px,py,imom,is) = - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                              * real(   Al(px,py) * mom(qx,qy,0,is) * mom(mx,my,1,is)  &
                                      - Al(qx,qy) * mom(px,py,0,is) * mom(mx,my,1,is)  &
                                      + Al(px,py) * mom(qx,qy,1,is) * mom(mx,my,0,is)  &
                                      - Al(qx,qy) * mom(px,py,1,is) * mom(mx,my,0,is), kind=DP )
            jpqk_em(px,py,imom,is) = - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                              * real(   Al(qx,qy) * mom(mx,my,0,is) * mom(px,py,1,is)  &
                                      - Al(mx,my) * mom(qx,qy,0,is) * mom(px,py,1,is)  &
                                      + Al(qx,qy) * mom(mx,my,1,is) * mom(px,py,0,is)  &
                                      - Al(mx,my) * mom(qx,qy,1,is) * mom(px,py,0,is), kind=DP )
            jqkp_em(px,py,imom,is) = - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                              * real(   Al(mx,my) * mom(px,py,0,is) * mom(qx,qy,1,is)  &
                                      - Al(px,py) * mom(mx,my,0,is) * mom(qx,qy,1,is)  &
                                      + Al(mx,my) * mom(px,py,1,is) * mom(qx,qy,0,is)  &
                                      - Al(px,py) * mom(mx,my,1,is) * mom(qx,qy,0,is), kind=DP )
          end do
        end do
      imom = 1
        ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          qy = - py - my
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            qx = - px - mx
            jkpq_em(px,py,imom,is) = - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                              * real(   Al(px,py) * mom(qx,qy,1,is) * mom(mx,my,2,is)  &
                                      - Al(qx,qy) * mom(px,py,1,is) * mom(mx,my,2,is)  &
                                      + Al(px,py) * mom(qx,qy,2,is) * mom(mx,my,1,is)  &
                                      - Al(qx,qy) * mom(px,py,2,is) * mom(mx,my,1,is), kind=DP )
            jpqk_em(px,py,imom,is) = - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                              * real(   Al(qx,qy) * mom(mx,my,1,is) * mom(px,py,2,is)  &
                                      - Al(mx,my) * mom(qx,qy,1,is) * mom(px,py,2,is)  &
                                      + Al(qx,qy) * mom(mx,my,2,is) * mom(px,py,1,is)  &
                                      - Al(mx,my) * mom(qx,qy,2,is) * mom(px,py,1,is), kind=DP )
            jqkp_em(px,py,imom,is) = - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                              * real(   Al(mx,my) * mom(px,py,1,is) * mom(qx,qy,2,is)  &
                                      - Al(px,py) * mom(mx,my,1,is) * mom(qx,qy,2,is)  &
                                      + Al(mx,my) * mom(px,py,2,is) * mom(qx,qy,1,is)  &
                                      - Al(px,py) * mom(mx,my,2,is) * mom(qx,qy,1,is), kind=DP )
          end do
        end do
      imom = 2
        ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is))) * 0.5_DP
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          qy = - py - my
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            qx = - px - mx
            jkpq_em(px,py,imom,is) = - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                              * real(   Al(px,py) * mom(qx,qy,2,is) * mom(mx,my,4,is)  &
                                      - Al(qx,qy) * mom(px,py,2,is) * mom(mx,my,4,is)  &
                                      + Al(px,py) * mom(qx,qy,4,is) * mom(mx,my,2,is)  &
                                      - Al(qx,qy) * mom(px,py,4,is) * mom(mx,my,2,is), kind=DP )
            jpqk_em(px,py,imom,is) = - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                              * real(   Al(qx,qy) * mom(mx,my,2,is) * mom(px,py,4,is)  &
                                      - Al(mx,my) * mom(qx,qy,2,is) * mom(px,py,4,is)  &
                                      + Al(qx,qy) * mom(mx,my,4,is) * mom(px,py,2,is)  &
                                      - Al(mx,my) * mom(qx,qy,4,is) * mom(px,py,2,is), kind=DP )
            jqkp_em(px,py,imom,is) = - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                              * real(   Al(mx,my) * mom(px,py,2,is) * mom(qx,qy,4,is)  &
                                      - Al(px,py) * mom(mx,my,2,is) * mom(qx,qy,4,is)  &
                                      + Al(mx,my) * mom(px,py,4,is) * mom(qx,qy,2,is)  &
                                      - Al(px,py) * mom(mx,my,4,is) * mom(qx,qy,2,is), kind=DP )
          end do
        end do
      imom = 3
        ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          qy = - py - my
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            qx = - px - mx
            jkpq_em(px,py,imom,is) = - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                              * real(   Al(px,py) * mom(qx,qy,3,is) * mom(mx,my,5,is)  &
                                      - Al(qx,qy) * mom(px,py,3,is) * mom(mx,my,5,is)  &
                                      + Al(px,py) * mom(qx,qy,5,is) * mom(mx,my,3,is)  &
                                      - Al(qx,qy) * mom(px,py,5,is) * mom(mx,my,3,is), kind=DP )
            jpqk_em(px,py,imom,is) = - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                              * real(   Al(qx,qy) * mom(mx,my,3,is) * mom(px,py,5,is)  &
                                      - Al(mx,my) * mom(qx,qy,3,is) * mom(px,py,5,is)  &
                                      + Al(qx,qy) * mom(mx,my,5,is) * mom(px,py,3,is)  &
                                      - Al(mx,my) * mom(qx,qy,5,is) * mom(px,py,3,is), kind=DP )
            jqkp_em(px,py,imom,is) = - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                              * real(   Al(mx,my) * mom(px,py,3,is) * mom(qx,qy,5,is)  &
                                      - Al(px,py) * mom(mx,my,3,is) * mom(qx,qy,5,is)  &
                                      + Al(mx,my) * mom(px,py,5,is) * mom(qx,qy,3,is)  &
                                      - Al(px,py) * mom(mx,my,5,is) * mom(qx,qy,3,is), kind=DP )
          end do
        end do

    end do

    tk_es(:,:) = 0._DP
    tk_es_pos(:,:) = 0._DP
    tk_es_neg(:,:) = 0._DP
    tk_em(:,:) = 0._DP
    tk_em_pos(:,:) = 0._DP
    tk_em_neg(:,:) = 0._DP
    do is = 0, ns-1
      do imom = 0, nmom-1
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            tk_es(imom,is) = tk_es(imom,is) + jkpq_es(px,py,imom,is)
            if ( jkpq_es(px,py,imom,is) < 0._DP ) then
              tk_es_neg(imom,is) = tk_es_neg(imom,is) + jkpq_es(px,py,imom,is)
            else
              tk_es_pos(imom,is) = tk_es_pos(imom,is) + jkpq_es(px,py,imom,is)
            end if
            tk_em(imom,is) = tk_em(imom,is) + jkpq_em(px,py,imom,is)
            if ( jkpq_em(px,py,imom,is) < 0._DP ) then
              tk_em_neg(imom,is) = tk_em_neg(imom,is) + jkpq_em(px,py,imom,is)
            else
              tk_em_pos(imom,is) = tk_em_pos(imom,is) + jkpq_em(px,py,imom,is)
            end if
          end do
        end do
      end do
    end do
!----------------------------

    mx = diag_mx
    my = diag_my
    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) my
    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop
    write( crow, '(i0)' ) 2*nx+1
!--- ascii output for gnuplot 1 ---
    open( ofldtrninkxky, file="./data/fluiddetailtransinkxky_x"//cmx//"y"//cmy//"z"//ciz//"_t"//cloop//".dat" )
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop,", time=",time
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,  ",   kx=",kx(mx)
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "#  gmy=",my,  ",   ky=",gky(my)
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "#  giz=",giz, ",   zz=",gzz(giz)
      write( ofldtrninkxky, "(99a17)" ) "#          [1] px", "[2] py",  &
                                        "[3] jkpq_es(e)", "[4] jpqk_es(e)", "[5] jqkp_es(e)", &
                                        "[6] jkpq_es(i)", "[7] jpqk_es(i)", "[8] jqkp_es(i)", &
                                        "[9] jkpq_em(e)", "[10] jpqk_em(e)", "[11] jqkp_em(e)", &
                                        "[12] jkpq_em(i)", "[13] jpqk_em(i)", "[14] jqkp_em(i)"
      do py = -global_ny, global_ny
        do px = -nx, nx
          write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) kx(px), wky(py)
          do is = 0, ns-1
            write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) sum(jkpq_es(px,py,:,is)),  &
                                                                sum(jpqk_es(px,py,:,is)),  &
                                                                sum(jqkp_es(px,py,:,is))
          end do
          do is = 0, ns-1
            write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) sum(jkpq_em(px,py,:,is)),  &
                                                                sum(jpqk_em(px,py,:,is)),  &
                                                                sum(jqkp_em(px,py,:,is))
          end do
          write( ofldtrninkxky, * )
        end do
        write( ofldtrninkxky, * )
      end do
    close( ofldtrninkxky )
!--- ascii output for gnuplot 2 ---
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" ) time
    do is = 0, ns-1
      write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" )  &
                                sum(tk_es(:,is)), sum(tk_es_pos(:,is)), sum(tk_es_neg(:,is))
    end do
    do is = 0, ns-1
      write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" )  &
                                sum(tk_em(:,is)), sum(tk_em_pos(:,is)), sum(tk_em_neg(:,is))
    end do
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, * )
!--------------------------------

    !%%% debug
    !is = 0
    !px = 2
    !py = 3
    !qx = - px - mx
    !qy = - py - my
    !write(89,*) time, sum(jkpq_es(px,py,:,is)), sum(jkpq_es(qx,qy,:,is)), &
    !                  sum(jpqk_es(px,py,:,is)), sum(jqkp_es(qx,qy,:,is)), &
    !                  sum(jpqk_es(qx,qy,:,is)), sum(jqkp_es(px,py,:,is))
    !is = 1
    !write(88,*) time, sum(jkpq_es(px,py,:,is)), sum(jkpq_es(qx,qy,:,is)), &
    !                  sum(jpqk_es(px,py,:,is)), sum(jqkp_es(qx,qy,:,is)), &
    !                  sum(jpqk_es(qx,qy,:,is)), sum(jqkp_es(px,py,:,is))

END SUBROUTINE fluiddetailtrans_izloop_calc


SUBROUTINE fluiddetailtrans_izloop_open( diag_mx, diag_my, giz )
!-------------------------------------------------------------------------------
!
!     file open
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: diag_mx, diag_my, giz

  character(len=4) :: cmx, cmy, ciz

    write( cmx, fmt="(i4.4)" ) diag_mx
    write( cmy, fmt="(i4.4)" ) diag_my
    write( ciz, fmt="(i4.4)" ) giz

    open( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my,  &
              file="./data/fluiddetailtransinkxky_x"//cmx//"y"//cmy//"z"//ciz//".dat" )
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99a17)" ) "#        [1] time",   &
                       "[2] tk_es(ele)", "[3] tk_es_pos(ele)", "[4] tk_es_neg(ele)",  &
                       "[5] tk_es(ion)", "[6] tk_es_pos(ion)", "[7] tk_es_neg(ion)",  &
                      "[8] tk_em(ele)", "[9] tk_em_pos(ele)", "[10] tk_em_neg(ele)",  &
                      "[11] tk_em(ion)", "[12] tk_em_pos(ion)", "[13] tk_em_neg(ion)"

END SUBROUTINE fluiddetailtrans_izloop_open


SUBROUTINE fluiddetailtrans_izloop_close( diag_mx, diag_my )
!-------------------------------------------------------------------------------
!
!     file close
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: diag_mx, diag_my

    close( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my )

END SUBROUTINE fluiddetailtrans_izloop_close


SUBROUTINE fluiddetailtrans_loop( loop )
!-------------------------------------------------------------------------------
!
!     Detailed triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop, rb_Al_izloop, rb_mom_izimomisloop
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, gksq, gomg, g0, g1, gky

  integer, intent(in) :: loop 

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_ny:global_ny,-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP),  &
    dimension(-nx:nx,-global_ny:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom
  real(kind=DP) :: bb
  integer :: mx, my, iz, imom, is

!--- initialize ---
    call rb_phi_gettime( loop, time )
    do iz = -global_nz, global_nz-1
      call rb_phi_izloop( iz, loop, phi(:,0:global_ny,iz) )
      call rb_Al_izloop( iz, loop, Al(:,0:global_ny,iz) )
    end do
    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do imom = 0, nmom-1
          call rb_mom_izimomisloop( iz, imom, is, loop, mom(:,0:global_ny,iz,imom,is) )
        end do
      end do
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
    do iz = -global_nz, global_nz-1
      do my = 1, global_ny
        do mx = -nx, nx
          phi(-mx,-my,iz) = conjg( phi(mx,my,iz) )
          Al(-mx,-my,iz) = conjg( Al(mx,my,iz) )
        end do
      end do
    end do
    do is = 0, ns-1
      do imom = 0, nmom-1
        do iz = -global_nz, global_nz-1
          do my = 1, global_ny
            do mx = -nx, nx
              mom(-mx,-my,iz,imom,is) = conjg( mom(mx,my,iz,imom,is) )
            end do
          end do
        end do
      end do
    end do
    do my = 0, global_ny 
      wky(-my) = - gky(my)
      wky(my) = gky(my)
    end do
!------------------

    call fluiddetailtrans_loop_calc( loop,  0, 150, time, phi, Al, mom ) ! ETG
    call fluiddetailtrans_loop_calc( loop,  0,  66, time, phi, Al, mom ) ! Streamer
    call fluiddetailtrans_loop_calc( loop,  0,   5, time, phi, Al, mom ) ! ITG
    call fluiddetailtrans_loop_calc( loop,  0,   3, time, phi, Al, mom ) ! Turb
    call fluiddetailtrans_loop_calc( loop, 66,   0, time, phi, Al, mom ) ! High ZF
    call fluiddetailtrans_loop_calc( loop, 24,   0, time, phi, Al, mom ) ! Med ZF
    call fluiddetailtrans_loop_calc( loop,  1,   0, time, phi, Al, mom ) ! Low ZF


END SUBROUTINE fluiddetailtrans_loop


SUBROUTINE fluiddetailtrans_loop_calc( loop, diag_mx, diag_my, time, phi, Al, mom )
!-------------------------------------------------------------------------------
!
!     Detailed triad transfer in fluid approximation
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : Anum, Znum, sgn, fcs, tau, kx, gky, grootg, cfsrf

  integer, intent(in) :: loop, diag_mx, diag_my
  real(kind=DP), intent(in) :: time
  complex(kind=DP), intent(in),  &
    dimension(-nx:nx,-global_ny:global_ny,-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP), intent(in),  &
    dimension(-nx:nx,-global_ny:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom

  real(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:nmom-1,0:ns-1) ::  &
                             jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em
  real(kind=DP), dimension(0:nmom-1,0:ns-1) :: tk_es, tk_es_pos, tk_es_neg,  &
                                               tk_em, tk_em_pos, tk_em_neg
  real(kind=DP) :: ceff, fct
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop, crow
  integer :: mx, my, px, py, qx, qy, iz, imom, is

!--- calc. detail transfer ---
    jkpq_es(:,:,:,:) = 0._DP
    jpqk_es(:,:,:,:) = 0._DP
    jqkp_es(:,:,:,:) = 0._DP
    mx = diag_mx
    my = diag_my
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf

      do is = 0, ns-1
        do imom = 0, nmom-1

          if ( imom == 0 .or. imom == 1 .or. imom == 3 .or. imom == 5 ) then
            ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is)
          else if ( imom == 2 ) then
            ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * 0.5_DP
          else if ( imom == 4 ) then
            ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * 0.166666666666666_DP
          else
            write( *, * ) "nmom is wrong.", nmom, imom
            stop
          end if
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              jkpq_es(px,py,imom,is) = jkpq_es(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))        &
                                       * real( (   phi(px,py,iz) * mom(qx,qy,iz,imom,is)     &
                                                 - phi(qx,qy,iz) * mom(px,py,iz,imom,is) )   &
                                                 * mom(mx,my,iz,imom,is), kind=DP ) )
              jpqk_es(px,py,imom,is) = jpqk_es(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))        &
                                       * real( (   phi(qx,qy,iz) * mom(mx,my,iz,imom,is)     &
                                                 - phi(mx,my,iz) * mom(qx,qy,iz,imom,is) )   &
                                                 * mom(px,py,iz,imom,is), kind=DP ) )
              jqkp_es(px,py,imom,is) = jqkp_es(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))        &
                                       * real( (   phi(mx,my,iz) * mom(px,py,iz,imom,is)     &
                                                 - phi(px,py,iz) * mom(mx,my,iz,imom,is) )   &
                                                 * mom(qx,qy,iz,imom,is), kind=DP ) )
            end do
          end do

        end do
      end do

    end do

    jkpq_em(:,:,:,:) = 0._DP
    jpqk_em(:,:,:,:) = 0._DP
    jqkp_em(:,:,:,:) = 0._DP
    mx = diag_mx
    my = diag_my
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf

      do is = 0, ns-1
  
        imom = 0
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              jkpq_em(px,py,imom,is) = jkpq_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                     * real(   Al(px,py,iz) * mom(qx,qy,iz,0,is) * mom(mx,my,iz,1,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,0,is) * mom(mx,my,iz,1,is)  &
                             + Al(px,py,iz) * mom(qx,qy,iz,1,is) * mom(mx,my,iz,0,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,1,is) * mom(mx,my,iz,0,is), kind=DP ) )
              jpqk_em(px,py,imom,is) = jpqk_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                     * real(   Al(qx,qy,iz) * mom(mx,my,iz,0,is) * mom(px,py,iz,1,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,0,is) * mom(px,py,iz,1,is)  &
                             + Al(qx,qy,iz) * mom(mx,my,iz,1,is) * mom(px,py,iz,0,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,1,is) * mom(px,py,iz,0,is), kind=DP ) )
              jqkp_em(px,py,imom,is) = jqkp_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                     * real(   Al(mx,my,iz) * mom(px,py,iz,0,is) * mom(qx,qy,iz,1,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,0,is) * mom(qx,qy,iz,1,is)  &
                             + Al(mx,my,iz) * mom(px,py,iz,1,is) * mom(qx,qy,iz,0,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,1,is) * mom(qx,qy,iz,0,is), kind=DP ) )
            end do
          end do
        imom = 1
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              jkpq_em(px,py,imom,is) = jkpq_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                     * real(   Al(px,py,iz) * mom(qx,qy,iz,1,is) * mom(mx,my,iz,2,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,1,is) * mom(mx,my,iz,2,is)  &
                             + Al(px,py,iz) * mom(qx,qy,iz,2,is) * mom(mx,my,iz,1,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,2,is) * mom(mx,my,iz,1,is), kind=DP ) )
              jpqk_em(px,py,imom,is) = jpqk_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                     * real(   Al(qx,qy,iz) * mom(mx,my,iz,1,is) * mom(px,py,iz,2,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,1,is) * mom(px,py,iz,2,is)  &
                             + Al(qx,qy,iz) * mom(mx,my,iz,2,is) * mom(px,py,iz,1,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,2,is) * mom(px,py,iz,1,is), kind=DP ) )
              jqkp_em(px,py,imom,is) = jqkp_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                     * real(   Al(mx,my,iz) * mom(px,py,iz,1,is) * mom(qx,qy,iz,2,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,1,is) * mom(qx,qy,iz,2,is)  &
                             + Al(mx,my,iz) * mom(px,py,iz,2,is) * mom(qx,qy,iz,1,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,2,is) * mom(qx,qy,iz,1,is), kind=DP ) )
            end do
          end do
        imom = 2
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is))) * 0.5_DP
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              jkpq_em(px,py,imom,is) = jkpq_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                     * real(   Al(px,py,iz) * mom(qx,qy,iz,2,is) * mom(mx,my,iz,4,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,2,is) * mom(mx,my,iz,4,is)  &
                             + Al(px,py,iz) * mom(qx,qy,iz,4,is) * mom(mx,my,iz,2,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,4,is) * mom(mx,my,iz,2,is), kind=DP ) )
              jpqk_em(px,py,imom,is) = jpqk_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                     * real(   Al(qx,qy,iz) * mom(mx,my,iz,2,is) * mom(px,py,iz,4,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,2,is) * mom(px,py,iz,4,is)  &
                             + Al(qx,qy,iz) * mom(mx,my,iz,4,is) * mom(px,py,iz,2,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,4,is) * mom(px,py,iz,2,is), kind=DP ) )
              jqkp_em(px,py,imom,is) = jqkp_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                     * real(   Al(mx,my,iz) * mom(px,py,iz,2,is) * mom(qx,qy,iz,4,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,2,is) * mom(qx,qy,iz,4,is)  &
                             + Al(mx,my,iz) * mom(px,py,iz,4,is) * mom(qx,qy,iz,2,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,4,is) * mom(qx,qy,iz,2,is), kind=DP ) )
            end do
          end do
        imom = 3
          ceff = 0.5_DP * fcs(is) * tau(is) * Znum(is) * (- sqrt(tau(is) / Anum(is)))
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              jkpq_em(px,py,imom,is) = jkpq_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(px) * wky(qy) + wky(py) * kx(qx))  &
                     * real(   Al(px,py,iz) * mom(qx,qy,iz,3,is) * mom(mx,my,iz,5,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,3,is) * mom(mx,my,iz,5,is)  &
                             + Al(px,py,iz) * mom(qx,qy,iz,5,is) * mom(mx,my,iz,3,is)  &
                             - Al(qx,qy,iz) * mom(px,py,iz,5,is) * mom(mx,my,iz,3,is), kind=DP ) )
              jpqk_em(px,py,imom,is) = jpqk_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(qx) * wky(my) + wky(qy) * kx(mx))  &
                     * real(   Al(qx,qy,iz) * mom(mx,my,iz,3,is) * mom(px,py,iz,5,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,3,is) * mom(px,py,iz,5,is)  &
                             + Al(qx,qy,iz) * mom(mx,my,iz,5,is) * mom(px,py,iz,3,is)  &
                             - Al(mx,my,iz) * mom(qx,qy,iz,5,is) * mom(px,py,iz,3,is), kind=DP ) )
              jqkp_em(px,py,imom,is) = jqkp_em(px,py,imom,is) + fct * (  &! flux-surface average
                                     - ceff * (- kx(mx) * wky(py) + wky(my) * kx(px))  &
                     * real(   Al(mx,my,iz) * mom(px,py,iz,3,is) * mom(qx,qy,iz,5,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,3,is) * mom(qx,qy,iz,5,is)  &
                             + Al(mx,my,iz) * mom(px,py,iz,5,is) * mom(qx,qy,iz,3,is)  &
                             - Al(px,py,iz) * mom(mx,my,iz,5,is) * mom(qx,qy,iz,3,is), kind=DP ) )
            end do
          end do
  
      end do

    end do

    tk_es(:,:) = 0._DP
    tk_es_pos(:,:) = 0._DP
    tk_es_neg(:,:) = 0._DP
    tk_em(:,:) = 0._DP
    tk_em_pos(:,:) = 0._DP
    tk_em_neg(:,:) = 0._DP
    do is = 0, ns-1
      do imom = 0, nmom-1
        do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
          do px = max(-nx-mx,-nx), min(nx,nx-mx)
            tk_es(imom,is) = tk_es(imom,is) + jkpq_es(px,py,imom,is)
            if ( jkpq_es(px,py,imom,is) < 0._DP ) then
              tk_es_neg(imom,is) = tk_es_neg(imom,is) + jkpq_es(px,py,imom,is)
            else
              tk_es_pos(imom,is) = tk_es_pos(imom,is) + jkpq_es(px,py,imom,is)
            end if
            tk_em(imom,is) = tk_em(imom,is) + jkpq_em(px,py,imom,is)
            if ( jkpq_em(px,py,imom,is) < 0._DP ) then
              tk_em_neg(imom,is) = tk_em_neg(imom,is) + jkpq_em(px,py,imom,is)
            else
              tk_em_pos(imom,is) = tk_em_pos(imom,is) + jkpq_em(px,py,imom,is)
            end if
          end do
        end do
      end do
    end do
!----------------------------

    mx = diag_mx
    my = diag_my
    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) my
    write( cloop, fmt="(i8.8)" ) loop
    write( crow, '(i0)' ) 2*nx+1
!--- ascii output for gnuplot 1 ---
    open( ofldtrninkxky, file="./data/fluiddetailtransinkxky_x"//cmx//"y"//cmy//"_t"//cloop//".dat" )
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop,", time=",time
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,  ",   kx=",kx(mx)
      write( ofldtrninkxky, "(a,i17,a,g17.7e3)" ) "#  gmy=",my,  ",   ky=",gky(my)
      write( ofldtrninkxky, "(99a17)" ) "#          [1] px", "[2] py",  &
                                        "[3] jkpq_es(e)", "[4] jpqk_es(e)", "[5] jqkp_es(e)", &
                                        "[6] jkpq_es(i)", "[7] jpqk_es(i)", "[8] jqkp_es(i)", &
                                        "[9] jkpq_em(e)", "[10] jpqk_em(e)", "[11] jqkp_em(e)", &
                                        "[12] jkpq_em(i)", "[13] jpqk_em(i)", "[14] jqkp_em(i)"
      do py = -global_ny, global_ny
        do px = -nx, nx
          write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) kx(px), wky(py)
          do is = 0, ns-1
            write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) sum(jkpq_es(px,py,:,is)),  &
                                                                sum(jpqk_es(px,py,:,is)),  &
                                                                sum(jqkp_es(px,py,:,is))
          end do
          do is = 0, ns-1
            write( ofldtrninkxky, "(99g17.7e3)", advance="NO" ) sum(jkpq_em(px,py,:,is)),  &
                                                                sum(jpqk_em(px,py,:,is)),  &
                                                                sum(jqkp_em(px,py,:,is))
          end do
          write( ofldtrninkxky, * )
        end do
        write( ofldtrninkxky, * )
      end do
    close( ofldtrninkxky )
!--- ascii output for gnuplot 2 ---
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" ) time
    do is = 0, ns-1
      write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" )  &
                                sum(tk_es(:,is)), sum(tk_es_pos(:,is)), sum(tk_es_neg(:,is))
    end do
    do is = 0, ns-1
      write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99e17.7e3)", advance="NO" )  &
                                sum(tk_em(:,is)), sum(tk_em_pos(:,is)), sum(tk_em_neg(:,is))
    end do
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, * )
!--------------------------------


END SUBROUTINE fluiddetailtrans_loop_calc


SUBROUTINE fluiddetailtrans_loop_open( diag_mx, diag_my )
!-------------------------------------------------------------------------------
!
!     file open
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: diag_mx, diag_my

  character(len=4) :: cmx, cmy

    write( cmx, fmt="(i4.4)" ) diag_mx
    write( cmy, fmt="(i4.4)" ) diag_my

    open( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my,  &
              file="./data/fluiddetailtransinkxky_x"//cmx//"y"//cmy//".dat" )
    write( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my, "(99a17)" ) "#        [1] time",   &
                       "[2] tk_es(ele)", "[3] tk_es_pos(ele)", "[4] tk_es_neg(ele)",  &
                       "[5] tk_es(ion)", "[6] tk_es_pos(ion)", "[7] tk_es_neg(ion)",  &
                      "[8] tk_em(ele)", "[9] tk_em_pos(ele)", "[10] tk_em_neg(ele)",  &
                      "[11] tk_em(ion)", "[12] tk_em_pos(ion)", "[13] tk_em_neg(ion)"

END SUBROUTINE fluiddetailtrans_loop_open


SUBROUTINE fluiddetailtrans_loop_close( diag_mx, diag_my )
!-------------------------------------------------------------------------------
!
!     file close
!                                                   (S. Maeyama, 8 July 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: diag_mx, diag_my

    close( ofldtrninkxky+diag_mx+(2*nx+1)*diag_my )

END SUBROUTINE fluiddetailtrans_loop_close


END MODULE out_fluiddetailtrans
