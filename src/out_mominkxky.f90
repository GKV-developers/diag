MODULE out_mominkxky
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinkxky, Alinkxky, mominkxky, fenergyinkxky, fluxinkxky, fpintinkxky
         

 CONTAINS


SUBROUTINE phiinkxky( loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (kx,ky) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : kx, gky
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: phi2
  character(len=8) :: cloop
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = 0.5_DP * abs(phi(mx,my,iz))**2
        end do
      end do
    end do
    call intgrl_thet(wr3, phi2)

    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/phiinkxky_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(99a17)" ) "#              kx","ky","<|phi|^2>"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my), phi2(mx,my)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE phiinkxky


SUBROUTINE Alinkxky( loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (kx,ky) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_geom, only : kx, gky
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: Al2
  character(len=8) :: cloop
  integer :: mx, my, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = 0.5_DP * abs(Al(mx,my,iz))**2
        end do
      end do
    end do
    call intgrl_thet(wr3, Al2)

    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/Alinkxky_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(99a17)" ) "#              kx","ky","<|Al|^2>"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my), Al2(mx,my)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE Alinkxky


SUBROUTINE mominkxky( is, loop )
!-------------------------------------------------------------------------------
!
!     Output moments(0:nmom-1) in (kx,ky) at is, loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_isloop
  use diag_geom, only : kx, gky
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: mom
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:nmom-1) :: mom2
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, iz, imom

    call rb_mom_gettime( loop, time )
    call rb_mom_isloop( is, loop, mom )

    do imom = 0, nmom-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            wr3(mx,my,iz) = 0.5_DP * abs(mom(mx,my,iz,imom))**2
          end do
        end do
      end do
      call intgrl_thet(wr3, mom2(:,:,imom))
    end do

    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/mominkxky_s"//cis//"_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "#   is=",is,   ", imom=",imom
      write( omominkxky, "(99a17)" ) "#              kx","ky",                 &
                                     "<|dens|^2>","<|upara|^2>","<|ppara|^2>", &
                                     "<|pperp|^2>","<|qlpara|^2>","<|qlperp|^2>"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my), mom2(mx,my,:)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE mominkxky


SUBROUTINE fenergyinkxky( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (kx,ky) at loop
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop
  use diag_geom, only : kx, gky, fct_e_energy, fct_m_energy
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, Al
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: e_energy, m_energy
  character(len=8) :: cloop
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call rb_Al_loop( loop, Al )

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = 0.5_DP * fct_e_energy(mx,my,iz) * abs(phi(mx,my,iz))**2
        end do
      end do
    end do
    call intgrl_thet(wr3, e_energy)

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = 0.5_DP * fct_m_energy(mx,my,iz) * abs(Al(mx,my,iz))**2
        end do
      end do
    end do
    call intgrl_thet(wr3, m_energy)

    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/fenergyinkxky_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(99a17)" ) "#              kx","ky","W_E","W_M"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my),    &
                                             e_energy(mx,my), m_energy(mx,my)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE fenergyinkxky


SUBROUTINE fluxinkxky( is, loop )
!-------------------------------------------------------------------------------
!
!     Output flux in (kx,ky) at loop
!                                                   (S. Maeyama, 24 Aug 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_isloop
  use diag_geom, only : kx, gky
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1) :: mom
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: pflux_es, pflux_em, eflux_es, eflux_em
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call rb_Al_loop( loop, Al )
    call rb_mom_isloop( is, loop, mom )

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = real( - ui * gky(my) * phi(mx,my,iz)  &
                                * conjg( mom(mx,my,iz,0) ), kind=DP )
        end do
      end do
    end do
    call intgrl_thet(wr3, pflux_es)

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = real( ui * gky(my) * Al(mx,my,iz)  &
                                * conjg( mom(mx,my,iz,1) ), kind=DP )
        end do
      end do
    end do
    call intgrl_thet(wr3, pflux_em)

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = real( - ui * gky(my) * phi(mx,my,iz) &
                                * conjg( mom(mx,my,iz,2)   &
                                       + mom(mx,my,iz,3) ), kind=DP )
        end do
      end do
    end do
    call intgrl_thet(wr3, eflux_es)

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          wr3(mx,my,iz) = real( ui * gky(my) * Al(mx,my,iz)  &
                                * conjg( mom(mx,my,iz,4) &
                                       + mom(mx,my,iz,5) ), kind=DP )
        end do
      end do
    end do
    call intgrl_thet(wr3, eflux_em)

    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/fluxinkxky_s"//cis//"_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(a,i17)" ) "#   is=",is
      write( omominkxky, "(99a17)" ) "#              kx","ky","G_E","G_M","Q_E","Q_M"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my),  &
                                             pflux_es(mx,my), pflux_em(mx,my), &
                                             eflux_es(mx,my), eflux_em(mx,my)
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE fluxinkxky


SUBROUTINE fpintinkxky( loop )
!-------------------------------------------------------------------------------
!
!     Output field-particle interaction terms in (kx,ky) at loop
!                                                   (S. Maeyama, 17 Sep 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_imomisloop, loop_phi_sta, loop_phi_end
  use diag_geom, only : kx, gky, sgn, Znum
  use diag_intgrl, only : intgrl_thet

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, dAldt
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: dndt, upara
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,-2:2) :: wk
  real(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:ns-1) :: fpint_es, fpint_em
  character(len=8) :: cloop
  real(kind=DP) :: fct
  integer :: mx, my, iz, is

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    do is = 0, ns-1
      call rb_mom_imomisloop( 1, is, loop, upara(:,:,:,is) )
    end do

    if (loop == loop_phi_sta(snum)) then
      fct = 1._DP / (2._DP * dtout_ptn)
      call rb_Al_loop( loop  , wk(:,:,:, 0) )
      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            dAldt(mx,my,iz) = fct * ( -3._DP * wk(mx,my,iz, 0)  &
                                      +4._DP * wk(mx,my,iz,+1)  &
                                      -        wk(mx,my,iz,+2) )
          end do
        end do
      end do
      do is = 0, ns-1
        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              dndt(mx,my,iz,is) = fct * ( -3._DP * wk(mx,my,iz, 0)  &
                                          +4._DP * wk(mx,my,iz,+1)  &
                                          -        wk(mx,my,iz,+2) )
            end do
          end do
        end do
      end do
    else if (loop == loop_phi_sta(snum)+1) then
      fct = 1._DP / (6._DP * dtout_ptn)
      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
      call rb_Al_loop( loop  , wk(:,:,:, 0) )
      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            dAldt(mx,my,iz) = fct * ( -2._DP * wk(mx,my,iz,-1)  &
                                      -3._DP * wk(mx,my,iz, 0)  &
                                      +6._DP * wk(mx,my,iz,+1)  &
                                      -        wk(mx,my,iz,+2) )
          end do
        end do
      end do
      do is = 0, ns-1
        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              dndt(mx,my,iz,is) = fct * ( -2._DP * wk(mx,my,iz,-1)  &
                                          -3._DP * wk(mx,my,iz, 0)  &
                                          +6._DP * wk(mx,my,iz,+1)  &
                                          -        wk(mx,my,iz,+2) )
            end do
          end do
        end do
      end do
    else if (loop == loop_phi_end(enum)-1) then
      fct = 1._DP / (6._DP * dtout_ptn)
      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
      call rb_Al_loop( loop  , wk(:,:,:, 0) )
      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
                                      -6._DP * wk(mx,my,iz,-1)  &
                                      +3._DP * wk(mx,my,iz, 0)  &
                                      +2._DP * wk(mx,my,iz,+1) )
          end do
        end do
      end do
      do is = 0, ns-1
        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
                                          -6._DP * wk(mx,my,iz,-1)  &
                                          +3._DP * wk(mx,my,iz, 0)  &
                                          +2._DP * wk(mx,my,iz,+1) )
            end do
          end do
        end do
      end do
    else if (loop == loop_phi_end(enum)) then
      fct = 1._DP / (2._DP * dtout_ptn)
      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
      call rb_Al_loop( loop  , wk(:,:,:, 0) )
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
                                      -4._DP * wk(mx,my,iz,-1)  &
                                      +3._DP * wk(mx,my,iz, 0) )
          end do
        end do
      end do
      do is = 0, ns-1
        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
                                          -4._DP * wk(mx,my,iz,-1)  &
                                          +3._DP * wk(mx,my,iz, 0) )
            end do
          end do
        end do
      end do
    else ! loop_phi_sta(snum)+2 <= loop <= loop_phi_end(enum)
      fct = 1._DP / (12._DP * dtout_ptn)
      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
                                      -8._DP * wk(mx,my,iz,-1)  &
                                      +8._DP * wk(mx,my,iz,+1)  &
                                      -        wk(mx,my,iz,+2) )
          end do
        end do
      end do
      do is = 0, ns-1
        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
                                          -8._DP * wk(mx,my,iz,-1)  &
                                          +8._DP * wk(mx,my,iz,+1)  &
                                          -        wk(mx,my,iz,+2) )
            end do
          end do
        end do
      end do
    end if

    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            wr3(mx,my,iz) = real( - conjg( phi(mx,my,iz) ) * sgn(is) * Znum(is) &
                                   * dndt(mx,my,iz,is), kind=DP )
          end do
        end do
      end do
      call intgrl_thet(wr3, fpint_es(:,:,is))

      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            wr3(mx,my,iz) = real( - conjg( upara(mx,my,iz,is) ) * sgn(is) * Znum(is) &
                                   * dAldt(mx,my,iz), kind=DP )
          end do
        end do
      end do
      call intgrl_thet(wr3, fpint_em(:,:,is))
    end do

    write( cloop, fmt="(i8.8)" ) loop
    open( omominkxky, file="./data/fpintinkxky_t"//cloop//".dat" )
      write( omominkxky, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominkxky, "(99a17)" ) "#              kx","ky","R_sE", "R_sM"
      do my = 0, global_ny
        do mx = -nx, nx
          write( omominkxky, "(99g17.7e3)" ) kx(mx), gky(my),    &
                                           ( fpint_es(mx,my,is), &
                                             fpint_em(mx,my,is), is=0, ns-1 )
        end do
        write( omominkxky, * )
      end do
    close( omominkxky )

END SUBROUTINE fpintinkxky


END MODULE out_mominkxky
