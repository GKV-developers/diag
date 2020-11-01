MODULE out_zfdensity
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public zfd_open, zfd_close, zfdensity


 CONTAINS


SUBROUTINE zfd_open
    open( ozfdensity, file="./data/tzfdensity.dat" )
END SUBROUTINE zfd_open
SUBROUTINE zfd_close
    close( ozfdensity )
END SUBROUTINE zfd_close

SUBROUTINE zfdensity( loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, rb_Al_loop, rb_mom_imomisloop
  use diag_fft, only : fft_backward_x1
  use diag_geom, only : grootg, cfsrf, xx

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wc3
  complex(kind=DP), dimension(-nx:nx) :: phik, Alk
  complex(kind=DP), dimension(-nx:nx,0:nmom-1,0:ns-1) :: momk
  complex(kind=DP), dimension(0:2*nxw-1) :: phi, Al
  complex(kind=DP), dimension(0:2*nxw-1,0:nmom-1,0:ns-1) :: mom
  real(kind=DP) :: fct
  integer :: mx, my, iz, imom, is

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, wc3 )
    phik(:) = (0._DP, 0._DP)
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf
      my = 0
        do mx = -nx, nx
          phik(mx) = phik(mx) + fct * wc3(mx,my,iz)
        end do
    end do
    call rb_Al_loop( loop, wc3 )
    Alk(:) = (0._DP, 0._DP)
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf
      my = 0
        do mx = 0, nx
          Alk(mx) = Alk(mx) + fct * wc3(mx,my,iz)
        end do
    end do
    momk(:,:,:) = (0._DP, 0._DP)
    do is = 0, ns-1
      do imom = 0, nmom-1
        call rb_mom_imomisloop( imom, is, loop, wc3 )
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / cfsrf
          my = 0
            do mx = 0, nx
              momk(mx,imom,is) = momk(mx,imom,is) + fct * wc3(mx,my,iz)
            end do
        end do
      end do
    end do
    call fft_backward_x1( phik, phi )
    call fft_backward_x1( Alk, Al )
    do is = 0, ns-1
      do imom = 0, nmom-1
        call fft_backward_x1( momk(:,imom,is), mom(:,imom,is) )
      end do
    end do

    do mx = 0, 2*nxw-1
      write( ozfdensity, "(999g17.7e3)" ) time, xx(mx),  &
                                 real(phi(mx), kind=DP), &
                                 real(Al(mx), kind=DP),  &
                                 real(mom(mx,:,:), kind=DP)
    end do
    write( ozfdensity, * )


END SUBROUTINE zfdensity


END MODULE out_zfdensity
