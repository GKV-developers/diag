MODULE out_momintky
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!  != example to write moments in tky =
!    call momintky_open( "./data/wesintky.dat" )
!    do loop = loop_phi_sta(snum), loop_phi_end(enum), 10
!      call wesintky( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/wemintky.dat" )
!    do loop = loop_phi_sta(snum), loop_phi_end(enum), 10
!      call wemintky( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/engintky.dat" )
!    do loop = loop_phi_sta(snum), loop_phi_end(enum), 10
!      call engintky( loop )
!    end do
!    call momintky_close
!    call momintky_open( "./data/menintky.dat" )
!    do loop = loop_phi_sta(snum), loop_phi_end(enum), 10
!      call menintky( loop )
!    end do
!    call momintky_close
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public momintky_open, momintky_close, &
         wesintky, wemintky, engintky, menintky, &!, fpintintky
         wesintky_parity, wemintky_parity, engintky_parity, menintky_parity


 CONTAINS


SUBROUTINE momintky_open( filename )
  character(len=*), intent(in) :: filename
    open( omomintky, file=filename )
END SUBROUTINE momintky_open


SUBROUTINE momintky_close
    close( omomintky )
END SUBROUTINE momintky_close


SUBROUTINE wesintky( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!                                                   (S. Maeyama, 17 June 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : grootg, cfsrf, fct_e_energy

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: e_energy
  real(kind=DP), dimension(0:global_ny) :: wke
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )

    e_energy(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          e_energy(mx,my) = e_energy(mx,my)  &
                          + fct * fct_e_energy(mx,my,iz) * abs(phi(mx,my,iz))**2
        end do
      end do
    end do

    wke(:) = 0._DP
    my = 0
      do mx = 0, nx
        wke(my) = wke(my) + e_energy(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        wke(my) = wke(my) + e_energy(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(wke), wke(:)

END SUBROUTINE wesintky


SUBROUTINE wemintky( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!                                                   (S. Maeyama, 17 June 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_geom, only : grootg, cfsrf, fct_m_energy

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: m_energy
  real(kind=DP), dimension(0:global_ny) :: wkm
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )

    m_energy(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          m_energy(mx,my) = m_energy(mx,my)  &
                          + fct * fct_m_energy(mx,my,iz) * abs(Al(mx,my,iz))**2
        end do
      end do
    end do

    wkm(:) = 0._DP
    my = 0
      do mx = 0, nx
        wkm(my) = wkm(my) + m_energy(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        wkm(my) = wkm(my) + m_energy(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(wkm), wkm(:)

END SUBROUTINE wemintky


SUBROUTINE engintky( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!                                                   (S. Maeyama, 17 June 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : grootg, cfsrf

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: wr2
  real(kind=DP), dimension(0:global_ny) :: mode_e
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )

    wr2(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          wr2(mx,my) = wr2(mx,my) + fct * abs(phi(mx,my,iz))**2
        end do
      end do
    end do

    mode_e(:) = 0._DP
    my = 0
      do mx = 0, nx
        mode_e(my) = mode_e(my) + wr2(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        mode_e(my) = mode_e(my) + wr2(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(mode_e), mode_e(:)

END SUBROUTINE engintky


SUBROUTINE menintky( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!                                                   (S. Maeyama, 17 June 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_geom, only : grootg, cfsrf

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: wr2
  real(kind=DP), dimension(0:global_ny) :: mode_b
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )

    wr2(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          wr2(mx,my) = wr2(mx,my) + fct * abs(Al(mx,my,iz))**2
        end do
      end do
    end do

    mode_b(:) = 0._DP
    my = 0
      do mx = 0, nx
        mode_b(my) = mode_b(my) + wr2(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        mode_b(my) = mode_b(my) + wr2(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(mode_b), mode_b(:)

END SUBROUTINE menintky


!SUBROUTINE fpintintky( loop )
!!-------------------------------------------------------------------------------
!!
!!     Output field-particle interaction terms in (t,ky)
!!                                                   (S. Maeyama, 16 June 2016)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_imomisloop, loop_phi_sta, loop_phi_end
!  use diag_geom, only : grootg, cfsrf, sgn, Znum
!
!  integer, intent(in) :: loop
!
!  real(kind=DP) :: time
!  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, dAldt
!  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: dndt, upara
!  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,-2:2) :: wk
!  real(kind=DP), dimension(-nx:nx,0:global_ny,0:ns-1) :: fpint_es, fpint_em
!  real(kind=DP), dimension(0:global_ny,0:ns-1) :: wke, wkm
!  real(kind=DP) :: fct
!  integer :: mx, my, iz, is
!
!    call rb_phi_gettime( loop, time )
!    call rb_phi_loop( loop, phi )
!    do is = 0, ns-1
!      call rb_mom_imomisloop( 1, is, loop, upara(:,:,:,is) )
!    end do
!
!    if (loop == loop_phi_sta(snum)) then
!      fct = 1._DP / (2._DP * dtout_ptn)
!      call rb_Al_loop( loop  , wk(:,:,:, 0) )
!      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
!      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
!      do iz = -global_nz, global_nz-1
!        do my = 0, global_ny
!          do mx = -nx, nx
!            dAldt(mx,my,iz) = fct * ( -3._DP * wk(mx,my,iz, 0)  &
!                                      +4._DP * wk(mx,my,iz,+1)  &
!                                      -        wk(mx,my,iz,+2) )
!          end do
!        end do
!      end do
!      do is = 0, ns-1
!        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
!        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
!        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
!        do iz = -global_nz, global_nz-1
!          do my = 0, global_ny
!            do mx = -nx, nx
!              dndt(mx,my,iz,is) = fct * ( -3._DP * wk(mx,my,iz, 0)  &
!                                          +4._DP * wk(mx,my,iz,+1)  &
!                                          -        wk(mx,my,iz,+2) )
!            end do
!          end do
!        end do
!      end do
!    else if (loop == loop_phi_sta(snum)+1) then
!      fct = 1._DP / (6._DP * dtout_ptn)
!      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
!      call rb_Al_loop( loop  , wk(:,:,:, 0) )
!      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
!      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
!      do iz = -global_nz, global_nz-1
!        do my = 0, global_ny
!          do mx = -nx, nx
!            dAldt(mx,my,iz) = fct * ( -2._DP * wk(mx,my,iz,-1)  &
!                                      -3._DP * wk(mx,my,iz, 0)  &
!                                      +6._DP * wk(mx,my,iz,+1)  &
!                                      -        wk(mx,my,iz,+2) )
!          end do
!        end do
!      end do
!      do is = 0, ns-1
!        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
!        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
!        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
!        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
!        do iz = -global_nz, global_nz-1
!          do my = 0, global_ny
!            do mx = -nx, nx
!              dndt(mx,my,iz,is) = fct * ( -2._DP * wk(mx,my,iz,-1)  &
!                                          -3._DP * wk(mx,my,iz, 0)  &
!                                          +6._DP * wk(mx,my,iz,+1)  &
!                                          -        wk(mx,my,iz,+2) )
!            end do
!          end do
!        end do
!      end do
!    else if (loop == loop_phi_end(enum)-1) then
!      fct = 1._DP / (6._DP * dtout_ptn)
!      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
!      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
!      call rb_Al_loop( loop  , wk(:,:,:, 0) )
!      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
!      do iz = -global_nz, global_nz-1
!        do my = 0, global_ny
!          do mx = -nx, nx
!            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
!                                      -6._DP * wk(mx,my,iz,-1)  &
!                                      +3._DP * wk(mx,my,iz, 0)  &
!                                      +2._DP * wk(mx,my,iz,+1) )
!          end do
!        end do
!      end do
!      do is = 0, ns-1
!        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
!        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
!        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
!        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
!        do iz = -global_nz, global_nz-1
!          do my = 0, global_ny
!            do mx = -nx, nx
!              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
!                                          -6._DP * wk(mx,my,iz,-1)  &
!                                          +3._DP * wk(mx,my,iz, 0)  &
!                                          +2._DP * wk(mx,my,iz,+1) )
!            end do
!          end do
!        end do
!      end do
!    else if (loop == loop_phi_end(enum)) then
!      fct = 1._DP / (2._DP * dtout_ptn)
!      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
!      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
!      call rb_Al_loop( loop  , wk(:,:,:, 0) )
!      do iz = -global_nz, global_nz-1
!        do my = 0, global_ny
!          do mx = -nx, nx
!            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
!                                      -4._DP * wk(mx,my,iz,-1)  &
!                                      +3._DP * wk(mx,my,iz, 0) )
!          end do
!        end do
!      end do
!      do is = 0, ns-1
!        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
!        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
!        call rb_mom_imomisloop( 0, is, loop  , wk(:,:,:, 0) )
!        do iz = -global_nz, global_nz-1
!          do my = 0, global_ny
!            do mx = -nx, nx
!              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
!                                          -4._DP * wk(mx,my,iz,-1)  &
!                                          +3._DP * wk(mx,my,iz, 0) )
!            end do
!          end do
!        end do
!      end do
!    else ! loop_phi_sta(snum)+2 <= loop <= loop_phi_end(enum)
!      fct = 1._DP / (12._DP * dtout_ptn)
!      call rb_Al_loop( loop-2, wk(:,:,:,-2) )
!      call rb_Al_loop( loop-1, wk(:,:,:,-1) )
!      call rb_Al_loop( loop+1, wk(:,:,:,+1) )
!      call rb_Al_loop( loop+2, wk(:,:,:,+2) )
!      do iz = -global_nz, global_nz-1
!        do my = 0, global_ny
!          do mx = -nx, nx
!            dAldt(mx,my,iz) = fct * ( +        wk(mx,my,iz,-2)  &
!                                      -8._DP * wk(mx,my,iz,-1)  &
!                                      +8._DP * wk(mx,my,iz,+1)  &
!                                      -        wk(mx,my,iz,+2) )
!          end do
!        end do
!      end do
!      do is = 0, ns-1
!        call rb_mom_imomisloop( 0, is, loop-2, wk(:,:,:,-2) )
!        call rb_mom_imomisloop( 0, is, loop-1, wk(:,:,:,-1) )
!        call rb_mom_imomisloop( 0, is, loop+1, wk(:,:,:,+1) )
!        call rb_mom_imomisloop( 0, is, loop+2, wk(:,:,:,+2) )
!        do iz = -global_nz, global_nz-1
!          do my = 0, global_ny
!            do mx = -nx, nx
!              dndt(mx,my,iz,is) = fct * ( +        wk(mx,my,iz,-2)  &
!                                          -8._DP * wk(mx,my,iz,-1)  &
!                                          +8._DP * wk(mx,my,iz,+1)  &
!                                          -        wk(mx,my,iz,+2) )
!            end do
!          end do
!        end do
!      end do
!    end if
!
!    fpint_es(:,:,:) = 0._DP
!    fpint_em(:,:,:) = 0._DP
!    do is = 0, ns-1
!      do iz = -global_nz, global_nz-1
!        fct = grootg(iz) / cfsrf      ! flux-surface average
!        do my = 0, global_ny
!          do mx = -nx, nx
!            fpint_es(mx,my,is) = fpint_es(mx,my,is) + fct  &
!                               * real( - conjg( phi(mx,my,iz) ) * sgn(is) * Znum(is) &
!                                       * dndt(mx,my,iz,is), kind=DP )
!            fpint_em(mx,my,is) = fpint_em(mx,my,is) + fct  &
!                               * real( - conjg( upara(mx,my,iz,is) ) * sgn(is) * Znum(is) &
!                                       * dAldt(mx,my,iz), kind=DP )
!          end do
!        end do
!      end do
!    end do
!
!    wke(:,:) = 0._DP
!    wkm(:,:) = 0._DP
!    do is = 0, ns-1
!      my = 0
!        do mx = 0, nx
!          wke(my,is) = wke(my,is) + 2._DP * fpint_es(mx,my,is)
!          wkm(my,is) = wkm(my,is) + 2._DP * fpint_em(mx,my,is)
!        end do
!      do my = 1, global_ny
!        do mx = -nx, nx
!          wke(my,is) = wke(my,is) + 2._DP * fpint_es(mx,my,is)
!          wkm(my,is) = wkm(my,is) + 2._DP * fpint_em(mx,my,is)
!        end do
!      end do
!    end do
!
!    is = 0
!    write( 96, "(999g17.7e3)" ) time, sum(wke(:,is)), wke(:,is)
!    write( 97, "(999g17.7e3)" ) time, sum(wkm(:,is)), wkm(:,is)
!    is = 1
!    write( 98, "(999g17.7e3)" ) time, sum(wke(:,is)), wke(:,is)
!    write( 99, "(999g17.7e3)" ) time, sum(wkm(:,is)), wkm(:,is)
!
!END SUBROUTINE fpintintky


SUBROUTINE wesintky_parity( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!
!     NOTE: Up-down symmetry grootg(-z) = grootg(z) is assumed, so that
!        \int dz grootg(z) even(z) odd(z) = 0.
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : grootg, cfsrf, fct_e_energy

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, phi_even, phi_odd
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: e_even, e_odd
  real(kind=DP), dimension(0:global_ny) :: wke_even, wke_odd
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call parity_mom3(phi, phi_even, phi_odd)

    e_even(:,:) = 0._DP
    e_odd(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          e_even(mx,my) = e_even(mx,my)  &
                          + fct * fct_e_energy(mx,my,iz) * abs(phi_even(mx,my,iz))**2
          e_odd(mx,my) = e_odd(mx,my)  &
                          + fct * fct_e_energy(mx,my,iz) * abs(phi_odd(mx,my,iz))**2
        end do
      end do
    end do

    wke_even(:) = 0._DP
    wke_odd(:) = 0._DP
    my = 0
      do mx = 0, nx
        wke_even(my) = wke_even(my) + e_even(mx,my)
        wke_odd(my) = wke_odd(my) + e_odd(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        wke_even(my) = wke_even(my) + e_even(mx,my)
        wke_odd(my) = wke_odd(my) + e_odd(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(wke_even), wke_even(:), sum(wke_odd), wke_odd(:)

END SUBROUTINE wesintky_parity


SUBROUTINE wemintky_parity( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!
!     NOTE: Up-down symmetry grootg(-z) = grootg(z) is assumed, so that
!        \int dz grootg(z) even(z) odd(z) = 0.
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : grootg, cfsrf, fct_m_energy

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al, Al_even, Al_odd
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: m_even, m_odd
  real(kind=DP), dimension(0:global_ny) :: wkm_even, wkm_odd
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )
    call parity_mom3(Al, Al_even, Al_odd)

    m_even(:,:) = 0._DP
    m_odd(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          m_even(mx,my) = m_even(mx,my)  &
                          + fct * fct_m_energy(mx,my,iz) * abs(Al_even(mx,my,iz))**2
          m_odd(mx,my) = m_odd(mx,my)  &
                          + fct * fct_m_energy(mx,my,iz) * abs(Al_odd(mx,my,iz))**2
        end do
      end do
    end do

    wkm_even(:) = 0._DP
    wkm_odd(:) = 0._DP
    my = 0
      do mx = 0, nx
        wkm_even(my) = wkm_even(my) + m_even(mx,my)
        wkm_odd(my) = wkm_odd(my) + m_odd(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        wkm_even(my) = wkm_even(my) + m_even(mx,my)
        wkm_odd(my) = wkm_odd(my) + m_odd(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(wkm_even), wkm_even(:), sum(wkm_odd), wkm_odd(:)

END SUBROUTINE wemintky_parity


SUBROUTINE engintky_parity( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!
!     NOTE: Up-down symmetry grootg(-z) = grootg(z) is assumed, so that
!        \int dz grootg(z) even(z) odd(z) = 0.
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : grootg, cfsrf

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi, phi_even, phi_odd
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: wr2_even, wr2_odd
  real(kind=DP), dimension(0:global_ny) :: mode_even, mode_odd
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )
    call parity_mom3(phi, phi_even, phi_odd)

    wr2_even(:,:) = 0._DP
    wr2_odd(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          wr2_even(mx,my) = wr2_even(mx,my) + fct * abs(phi_even(mx,my,iz))**2
          wr2_odd(mx,my) = wr2_odd(mx,my) + fct * abs(phi_odd(mx,my,iz))**2
        end do
      end do
    end do

    mode_even(:) = 0._DP
    mode_odd(:) = 0._DP
    my = 0
      do mx = 0, nx
        mode_even(my) = mode_even(my) + wr2_even(mx,my)
        mode_odd(my) = mode_odd(my) + wr2_odd(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        mode_even(my) = mode_even(my) + wr2_even(mx,my)
        mode_odd(my) = mode_odd(my) + wr2_odd(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(mode_even), mode_even(:), sum(mode_odd), mode_odd(:)

END SUBROUTINE engintky_parity


SUBROUTINE menintky_parity( loop )
!-------------------------------------------------------------------------------
!
!     Output field energy in (t,ky)
!
!     NOTE: Up-down symmetry grootg(-z) = grootg(z) is assumed, so that
!        \int dz grootg(z) even(z) odd(z) = 0.
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop
  use diag_parity, only : parity_mom3
  use diag_geom, only : grootg, cfsrf

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al, Al_even, Al_odd
  real(kind=DP), dimension(-nx:nx,0:global_ny) :: wr2_even, wr2_odd
  real(kind=DP), dimension(0:global_ny) :: mode_even, mode_odd
  real(kind=DP) :: fct
  integer :: mx, my, iz

    call rb_Al_gettime( loop, time )
    call rb_Al_loop( loop, Al )
    call parity_mom3(Al, Al_even, Al_odd)

    wr2_even(:,:) = 0._DP
    wr2_odd(:,:) = 0._DP
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf      ! flux-surface average
      do my = 0, global_ny
        do mx = -nx, nx
          wr2_even(mx,my) = wr2_even(mx,my) + fct * abs(Al_even(mx,my,iz))**2
          wr2_odd(mx,my) = wr2_odd(mx,my) + fct * abs(Al_odd(mx,my,iz))**2
        end do
      end do
    end do

    mode_even(:) = 0._DP
    mode_odd(:) = 0._DP
    my = 0
      do mx = 0, nx
        mode_even(my) = mode_even(my) + wr2_even(mx,my)
        mode_odd(my) = mode_odd(my) + wr2_odd(mx,my)
      end do
    do my = 1, global_ny
      do mx = -nx, nx
        mode_even(my) = mode_even(my) + wr2_even(mx,my)
        mode_odd(my) = mode_odd(my) + wr2_odd(mx,my)
      end do
    end do

    write( omomintky, "(999g17.7e3)" ) time, sum(mode_even), mode_even(:), sum(mode_odd), mode_odd(:)

END SUBROUTINE menintky_parity


END MODULE out_momintky
