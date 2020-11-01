MODULE out_mominfreq
!-------------------------------------------------------------------------------
!
!     Frequency analysis
!                                                   (S. Maeyama, 3 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinfreq


 CONTAINS


SUBROUTINE phiinfreq( mx, gmy, giz, trange )
!-------------------------------------------------------------------------------
!
!     Frequency analysis for (mx,gmy,giz)
!                                                   (S. Maeyama, 3 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_getalltime, rb_phi_mxmyiz, loop_phi_sta, loop_phi_end
  use diag_geom, only : kx, gky, gzz
  use diag_fft, only : fft_time2freq_init, fft_time2freq

  integer, intent(in) :: mx, gmy, giz
  real(kind=DP), intent(in) :: trange

  real(kind=DP), dimension(:), allocatable :: time, omega
  complex(kind=DP), dimension(:), allocatable :: phimxmyiz, sample, spectrum
  real(kind=DP) :: period, omega_nyq
  integer :: nsample
  integer :: isample, loop
  character(len=4) :: cmx, cmy, ciz

    allocate( time(loop_phi_sta(snum):loop_phi_end(enum)) )
    allocate( phimxmyiz(loop_phi_sta(snum):loop_phi_end(enum)) )

    call rb_phi_getalltime( time )
    call rb_phi_mxmyiz( mx, gmy, giz, phimxmyiz )

   !--- Normalize ---
   ! do loop = loop_phi_sta(snum), loop_phi_end(enum)
   !   phimxmyiz(loop) = phimxmyiz(loop) / abs(phimxmyiz(loop))
   ! end do
   !-----------------
   !--- Resampling with equi-distant grid may be favorable. ---
   !
   !-----------------------------------------------------------

    nsample = int(trange / dtout_ptn)
    period = dtout_ptn * real( nsample, kind=DP )
    omega_nyq = pi / dtout_ptn
    allocate( sample(0:nsample-1) )
    allocate( spectrum(0:nsample-1) )
    allocate( omega(0:nsample-1) )

    do isample = 0, nsample-1
      omega(isample) = 2._DP * pi * real(isample, kind=DP) / period
                    != real(isample) * ( 2._DP * omega_nyq / real(nsample) )
    end do

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( ciz, fmt="(i4.4)" ) giz
    open(omominfreq, file="./data/phiinfreq_x"//cmx//"y"//cmy//"z"//ciz//".dat")
      write( omominfreq, "(a,i17,a,g17.7e3)" ) "#  mx=",mx,  ", kx=",kx(mx)
      write( omominfreq, "(a,i17,a,g17.7e3)" ) "# gmy=",gmy, ", ky=",gky(gmy)
      write( omominfreq, "(a,i17,a,g17.7e3)" ) "# giz=",giz, ", zz=",gzz(giz)
      write( omominfreq, "(a,i17,a,g17.7e3)" )  &
                                 "# nsample=",nsample, ", period=",period
      write( omominfreq, "(99a17)" ) "#            time","omega","spectrum"

      call fft_time2freq_init( nsample )
      do loop = loop_phi_sta(snum), loop_phi_end(enum)-nsample, 1+nsample/10
        sample(0:nsample-1) = phimxmyiz(loop:loop+nsample-1)

       !call fft_time2freq( nsample, sample, spectrum )
       !call fft_time2freq( nsample, gaussian_window(nsample,sample), spectrum )
       !call fft_time2freq( nsample, hanning_window(nsample,sample), spectrum )
        call fft_time2freq( nsample, hamming_window(nsample,sample), spectrum )
       !call fft_time2freq( nsample, blackman_window(nsample,sample), spectrum )

       !--- Assume expansions in exp(-i*omega*t)  ---
       !--- where -omega_nyq < omega < omega_nyq. ---
        do isample = nsample/2, nsample-1
          write( omominfreq, "(99g17.7e3)" ) time(loop+nsample/2),  &
                                   2._DP*omega_nyq-omega(isample),  &
                                        abs(spectrum(isample))**2
        end do
        do isample = 0, nsample/2-1
          write( omominfreq, "(99g17.7e3)" ) time(loop+nsample/2),  &
                                                  -omega(isample),  &
                                        abs(spectrum(isample))**2
        end do
        write( omominfreq, * )
      end do

    close( omominfreq )
                                        !--- for debug ---
                                        !do loop=loop_phi_sta(snum),loop_phi_end(enum)
                                        !  write(*,*) time(loop),             &
                                        !             real(phimxmyiz(loop)),  &
                                        !             aimag(phimxmyiz(loop))
                                        !end do
                                        !-----------------
    deallocate( time )
    deallocate( phimxmyiz )
    deallocate( sample )
    deallocate( spectrum )
    deallocate( omega )

END SUBROUTINE phiinfreq


!SUBROUTINE testwindow
!!-------------------------------------------------------------------------------
!!
!!     Test window functions
!!                                                   (S. Maeyama, 4 March 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_fft, only : fft_time2freq
!
!  integer, parameter :: nsample = 100
!  complex(kind=DP), dimension(0:nsample-1) :: sample, window, spectrum
!
!  integer :: isample
!
!    sample(:) = ( 1._DP, 0._DP )
!    open(6, file="rectw.dat")
!    window(:) = sample(:)
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="gausw.dat")
!    window(:) = gaussian_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="hannw.dat")
!    window(:) = hanning_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="hammw.dat")
!    window(:) = hamming_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="blckw.dat")
!    window(:) = blackman_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!
!    do isample = 0, nsample-1 
!      sample(isample) = sin( 3.5_DP *(2._DP*pi*real(isample)/real(nsample)) ) &
!                      + cos( 1._DP *(2._DP*pi*real(isample)/real(nsample)) ) &
!                      + sin( 10.2_DP *(2._DP*pi*real(isample)/real(nsample)) )
!    end do
!    open(6, file="rectd.dat")
!    window(:) = sample(:)
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="gausd.dat")
!    window(:) = gaussian_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="hannd.dat")
!    window(:) = hanning_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="hammd.dat")
!    window(:) = hamming_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!    open(6, file="blckd.dat")
!    window(:) = blackman_window( nsample, sample )
!    do isample = 0, nsample-1
!      write(6,*) isample, real(window(isample))
!    end do
!    close(6)
!
!    open(6, file="rects.dat")
!    call fft_time2freq( nsample, sample, spectrum )
!    do isample = 0, nsample-1
!      write(6,*) isample, abs(spectrum(isample))
!    end do
!    close(6)
!    open(6, file="gauss.dat")
!    call fft_time2freq( nsample, gaussian_window(nsample, sample), spectrum )
!    do isample = 0, nsample-1
!      write(6,*) isample, abs(spectrum(isample))
!    end do
!    close(6)
!    open(6, file="hanns.dat")
!    call fft_time2freq( nsample, hanning_window(nsample, sample), spectrum )
!    do isample = 0, nsample-1
!      write(6,*) isample, abs(spectrum(isample))
!    end do
!    close(6)
!    open(6, file="hamms.dat")
!    call fft_time2freq( nsample, hamming_window(nsample, sample), spectrum )
!    do isample = 0, nsample-1
!      write(6,*) isample, abs(spectrum(isample))
!    end do
!    close(6)
!    open(6, file="blcks.dat")
!    call fft_time2freq( nsample, blackman_window(nsample, sample), spectrum )
!    do isample = 0, nsample-1
!      write(6,*) isample, abs(spectrum(isample))
!    end do
!    close(6)
!
!END SUBROUTINE testwindow


FUNCTION gaussian_window(nsample, sample)
!-------------------------------------------------------------------------------
!
!     Gaussian window
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  complex(kind=DP), dimension(0:nsample-1) :: gaussian_window

  real(kind=DP) :: x, sigma
  integer :: isample

    sigma = 0.2_DP * real( nsample, kind=DP )
    do isample = 0, nsample-1
      x = real( - nsample + 2 * isample, kind=DP ) / 2._DP
      gaussian_window(isample) = exp( - 0.5_DP * (x/sigma)**2 )  &
                               * sample(isample)
    end do

    RETURN

END FUNCTION gaussian_window


FUNCTION hanning_window(nsample, sample)
!-------------------------------------------------------------------------------
!
!     Hanning window
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  complex(kind=DP), dimension(0:nsample-1) :: hanning_window

  real(kind=DP) :: x
  integer :: isample

    do isample = 0, nsample-1
      x = 2._DP * pi * real( isample, kind=DP ) / real( nsample, kind=DP )
      hanning_window(isample) = ( 0.5_DP - 0.5_DP * cos( x ) )  &
                              * sample(isample)
    end do

    RETURN

END FUNCTION hanning_window


FUNCTION hamming_window(nsample, sample)
!-------------------------------------------------------------------------------
!
!     Hamming window
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  complex(kind=DP), dimension(0:nsample-1) :: hamming_window

  real(kind=DP) :: x
  integer :: isample

    do isample = 0, nsample-1
      x = 2._DP * pi * real( isample, kind=DP ) / real( nsample, kind=DP )
      hamming_window(isample) = ( 0.54_DP - 0.46_DP * cos( x ) )  &
                              * sample(isample)
    end do

    RETURN

END FUNCTION hamming_window


FUNCTION blackman_window(nsample, sample)
!-------------------------------------------------------------------------------
!
!     Blackman window
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  complex(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  complex(kind=DP), dimension(0:nsample-1) :: blackman_window

  real(kind=DP) :: x
  integer :: isample

    do isample = 0, nsample-1
      x = 2._DP * pi * real( isample, kind=DP ) / real( nsample, kind=DP )
      blackman_window(isample) = ( 0.42_DP - 0.5_DP * cos( x )             &
                                           + 0.08_DP * cos( 2._DP * x ) )  &
                               * sample(isample)
    end do

    RETURN

END FUNCTION blackman_window


END MODULE out_mominfreq
