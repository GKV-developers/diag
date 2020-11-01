MODULE out_realsp
!-------------------------------------------------------------------------------
!
!     Output phi(x,y) from pk(kx,ky)
!                                                   (S. Maeyama, 26 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public realsp_create, realsp_open, realsp_close,  &
         realsp_gettime, realsp_getalltime,         &
         realsp_phiirec, realsp_phiixiy, nrec_realsp

  integer :: nrec_realsp

 CONTAINS


SUBROUTINE realsp_create( giz, loop_sta, loop_end )
!-------------------------------------------------------------------------------
!
!     Output phi(x,y) from pk(kx,ky)
!                                                   (S. Maeyama, 26 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop
  use diag_fft, only : fft_backward_xy

  integer, intent(in) :: giz, loop_sta, loop_end

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: phikxky
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)   :: phixy
  character(len=4) :: ciz
  integer :: loop

    write( ciz, fmt="(i4.4)" ) giz
    open( ophireal, file="./data/phiintxy_z"//ciz//".dat", status="replace",  &
                    action="write", form="unformatted", access="stream" )
      do loop = loop_sta, loop_end
        call rb_phi_gettime( loop, time )
        call rb_phi_izloop( giz, loop, phikxky )
        call fft_backward_xy( phikxky, phixy )
        write( ophireal ) time, phixy
      end do
    close( ophireal )

END SUBROUTINE realsp_create


SUBROUTINE realsp_open( giz )
!-------------------------------------------------------------------------------
!
!     Open ophireal
!                                                   (S. Maeyama, 28 Feb. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: giz

  character(len=4) :: ciz
  integer(kind=8) :: filesize!, statb(13)

    write( ciz, fmt="(i4.4)" ) giz
    open( ophireal, file="./data/phiintxy_z"//ciz//".dat", status="old",  &
                    action="read", form="unformatted", access="stream" )

   !%%% for GNU fortran   %%%
   ! call stat( "./data/phiintxy_z"//ciz//".dat", statb )
   ! filesize = statb(8)
   !%%% for Intel fortran %%%
    inquire( file="./data/phiintxy_z"//ciz//".dat", size=filesize )
   !%%%%%%%%%%%%%%%%%%%%%%%%%
    nrec_realsp = int(filesize/(DP+(2*nxw)*(2*nyw)*DP), kind=4)
                                        !--- for debug ---
                                        write(*,*)  &
                                          "# filesize = ", filesize, &
                                          "# nrec_realsp = ", nrec_realsp
                                        !-----------------

END SUBROUTINE realsp_open


SUBROUTINE realsp_close
!-------------------------------------------------------------------------------
!
!     Close ophireal
!                                                   (S. Maeyama, 28 Feb. 2014)
!
!-------------------------------------------------------------------------------

    close( ophireal )

END SUBROUTINE realsp_close


SUBROUTINE realsp_gettime( irec, time )
!-------------------------------------------------------------------------------
!
!     Get time at irec
!                                                   (S. Maeyama, 27 Feb. 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: irec
  real(kind=DP), intent(out) :: time

  integer(kind=8) :: skipbyte

    skipbyte = (DP+(2*nxw)*(2*nyw)*DP) *(irec-1)
    read( ophireal, pos=skipbyte+1 ) time

END SUBROUTINE realsp_gettime


SUBROUTINE realsp_getalltime( time )
!-------------------------------------------------------------------------------
!
!     Get time(loop)
!                                                   (S. Maeyama, 27 Feb. 2014)
!
!-------------------------------------------------------------------------------
  real(kind=DP), dimension(1:nrec_realsp), intent(out) :: time

  integer(kind=8) :: skipbyte
  integer :: irec

    do irec = 1, nrec_realsp
      skipbyte = (DP+(2*nxw)*(2*nyw)*DP) *(irec-1)
      read( ophireal, pos=skipbyte+1 ) time(irec)
    end do
                                        !--- for debug ---
                                        !write(*,"(9f8.2)") time(1:9)
                                        !-----------------

END SUBROUTINE realsp_getalltime


SUBROUTINE realsp_phiirec( irec, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(ix,iy) at irec
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: irec
  real(kind=DP), intent(out), dimension(0:2*nxw-1,0:2*nyw-1) :: phi

  integer(kind=8) :: skipbyte

      skipbyte = (DP+(2*nxw)*(2*nyw)*DP) *(irec-1)  &
               + DP
      read( ophireal, pos=skipbyte+1 ) phi
                                        !--- for debug ---
                                        !do iy = 0, 2*nyw-1
                                        !do ix = 0, 2*nxw-1
                                        !write(98,*) ix, iy, phi(ix,iy)
                                        !end do
                                        !write(98,*)
                                        !end do
                                        !-----------------

END SUBROUTINE realsp_phiirec


SUBROUTINE realsp_phiixiy( ix, iy, phi )
!-------------------------------------------------------------------------------
!
!     Get phi(loop) at xx(ix), yy(iy)
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: ix, iy
  real(kind=DP), intent(out), dimension(1:nrec_realsp) :: phi

  integer(kind=8) :: skipbyte
  integer :: irec

    do irec = 1, nrec_realsp
      skipbyte = (DP+(2*nxw)*(2*nyw)*DP) *(irec-1)  &
               + DP + ( ix + (2*nxw)*iy )*DP
      read( ophireal, pos=skipbyte+1 ) phi(irec)
    end do

END SUBROUTINE realsp_phiixiy


END MODULE out_realsp
