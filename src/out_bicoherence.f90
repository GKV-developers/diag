MODULE out_bicoherence
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 13 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public bicoh_phicheck, bicoh_freqcheck, bicoh_bicoherence


 CONTAINS


SUBROUTINE bicoh_phicheck( loop )
!-------------------------------------------------------------------------------
!
!     Check phi(ix,iy) at irec
!                                                   (S. Maeyama, 3 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : xx, yy
  use out_realsp, only : realsp_open, realsp_close, realsp_gettime, realsp_phiirec

  integer, intent(in) :: loop 

  real(kind=DP) :: time
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: phi
  integer :: ix, iy, irec
  character(len=8) :: cloop

    call realsp_open( 0 )

    irec = loop + 1
    call realsp_gettime( irec, time )
    call realsp_phiirec( irec, phi )

    write( cloop, fmt="(i8.8)" ) loop
    open( omominxy, file="./data/phiinxy_t"//cloop//".dat" )
      write( omominxy, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominxy, "(99a17)" ) "#              xx","yy","phi"
      do iy = 0, 2*nyw-1
        do ix = 0, 2*nxw-1
          write( omominxy, "(99g17.7e3)" ) xx(ix), yy(iy), phi(ix,iy)
        end do
        write( omominxy, * )
      end do
    close( omominxy )

    call realsp_close

END SUBROUTINE bicoh_phicheck


SUBROUTINE bicoh_freqcheck( ix, iy, trange )
!-------------------------------------------------------------------------------
!
!     Frequency analysis for (ix,iy)
!                                                   (S. Maeyama, 3 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : xx, yy
  use diag_fft, only : fft_time2freq_r2c_init, fft_time2freq_r2c
  use out_realsp, only : realsp_open, realsp_close, nrec_realsp,  &
                         realsp_getalltime, realsp_phiixiy

  integer, intent(in) :: ix, iy
  real(kind=DP), intent(in) :: trange

  real(kind=DP), dimension(:), allocatable :: time, phi, sample, omega
  complex(kind=DP), dimension(:), allocatable :: spectrum
  real(kind=DP) :: period, omega_nyq
  integer :: nsample
  integer :: isample, irec
  character(len=4) :: cix, ciy

    call realsp_open( 0 )

    allocate( time(1:nrec_realsp) )
    allocate( phi(1:nrec_realsp) )

    call realsp_getalltime( time )
    call realsp_phiixiy( ix, iy, phi )

    nsample = int(trange / dtout_ptn); nsample = nsample - mod(nsample, 2) ! even number
    period = dtout_ptn * real( nsample, kind=DP )
    omega_nyq = pi / dtout_ptn
    allocate( sample(0:nsample-1) )
    allocate( omega(-nsample/2:nsample/2) )
    allocate( spectrum(-nsample/2:nsample/2) )

    do isample = -nsample/2, nsample/2
      omega(isample) = -2._DP * pi * real(isample, kind=DP) / period
                    != -real(isample) * ( 2._DP * omega_nyq / real(nsample) )
                    !--- Assume expansions in exp(-i*omega*t)  ---
                    !--- where -omega_nyq < omega < omega_nyq. ---
    end do

    write( cix, fmt="(i4.4)" ) ix
    write( ciy, fmt="(i4.4)" ) iy
    open(omominfreq, file="./data/phiinfreq_ix"//cix//"iy"//ciy//".dat")
      write( omominfreq, "(a,i17,a,g17.7e3)" ) "#  ix=",ix, ", xx=",xx(ix)
      write( omominfreq, "(a,i17,a,g17.7e3)" ) "#  iy=",iy, ", yy=",yy(iy)
      write( omominfreq, "(a,i17,a,g17.7e3)" )  &
                                 "# nsample=",nsample, ", period=",period
      write( omominfreq, "(99a17)" ) "#            time","omega","spectrum"

      call fft_time2freq_r2c_init( nsample )
      do irec = 1, nrec_realsp-nsample, 1+nsample/10
        sample(0:nsample-1) = phi(irec:irec+nsample-1)

       !call fft_time2freq_r2c( nsample, sample, spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, gaussian_window(nsample,sample), spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, hanning_window(nsample,sample), spectrum(0:nsample/2) )
        call fft_time2freq_r2c( nsample, hamming_window(nsample,sample), spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, blackman_window(nsample,sample), spectrum(0:nsample/2) )

        do isample = 1, nsample/2
          spectrum(-isample) = conjg( spectrum(isample) )
        end do

        do isample = -nsample/2, nsample/2
          write( omominfreq, "(99g17.7e3)" ) time(irec+nsample/2),  &
                                                   omega(isample),  &
                                        abs(spectrum(isample))**2
        end do
        write( omominfreq, * )
      end do

    close( omominfreq )
                                        !--- for debug ---
                                        !do loop=loop_phi_sta(snum),loop_phi_end(enum)
                                        !  write(*,*) time(loop),             &
                                        !             real(phi(loop)),  &
                                        !             aimag(phi(loop))
                                        !end do
                                        !-----------------
    deallocate( time )
    deallocate( phi )
    deallocate( sample )
    deallocate( spectrum )
    deallocate( omega )

    call realsp_close

END SUBROUTINE bicoh_freqcheck


FUNCTION gaussian_window(nsample, sample)
!-------------------------------------------------------------------------------
!
!     Gaussian window
!                                                   (S. Maeyama, 4 March 2014)
!
!-------------------------------------------------------------------------------
  integer, intent(in) :: nsample
  real(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  real(kind=DP), dimension(0:nsample-1) :: gaussian_window

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
  real(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  real(kind=DP), dimension(0:nsample-1) :: hanning_window

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
  real(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  real(kind=DP), dimension(0:nsample-1) :: hamming_window

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
  real(kind=DP), dimension(0:nsample-1), intent(in) :: sample
  real(kind=DP), dimension(0:nsample-1) :: blackman_window

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


SUBROUTINE bicoh_bicoherence( tsta, tend )
!-------------------------------------------------------------------------------
!
!     Bicoherence analysis ( statistical average is replaced by volume average )
!                                                   (S. Maeyama, 3 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_fft, only : fft_time2freq_r2c_init, fft_time2freq_r2c
  use out_realsp, only : realsp_open, realsp_close, nrec_realsp,  &
                         realsp_getalltime, realsp_phiixiy

  real(kind=DP), intent(in) :: tsta, tend

  real(kind=DP), dimension(:), allocatable :: time, phi, sample, omega
  complex(kind=DP), dimension(:), allocatable :: spectrum
  real(kind=DP), dimension(:), allocatable :: ampl
  real(kind=DP), dimension(:,:), allocatable :: biampl, bicoherence
  complex(kind=DP), dimension(:,:), allocatable :: bispectrum
  real(kind=DP) :: trange, period, omega_nyq
  integer :: nsample
  integer :: isample, irec, ix, iy, w, p, q, inum
  character(len=8) :: cnum
  character(len=4) :: csta, cend

    call realsp_open( 0 )

!--- initialize
    allocate( time(1:nrec_realsp) )
    allocate( phi(1:nrec_realsp) )

    call realsp_getalltime( time )

    trange = tend - tsta
    nsample = int(trange / dtout_ptn); nsample = nsample - mod(nsample, 2) ! even number
    period = dtout_ptn * real( nsample, kind=DP )
    omega_nyq = pi / dtout_ptn
    allocate( sample(0:nsample-1) )
    allocate( omega(-nsample/2:nsample/2) )
    allocate( spectrum(-nsample/2:nsample/2) )
    allocate( ampl(-nsample/2:nsample/2) )
    allocate( biampl(-nsample/2:nsample/2,-nsample/2:nsample/2) )
    allocate( bicoherence(-nsample/2:nsample/2,-nsample/2:nsample/2) )
    allocate( bispectrum(-nsample/2:nsample/2,-nsample/2:nsample/2) )

    bispectrum(:,:) = (0._DP, 0._DP)
    ampl(:) = 0._DP
    biampl(:,:) = 0._DP

    do isample = -nsample/2, nsample/2
      omega(isample) = -2._DP * pi * real(isample, kind=DP) / period
                    != -real(isample) * ( 2._DP * omega_nyq / real(nsample) )
                    !--- Assume expansions in exp(-i*omega*t)  ---
                    !--- where -omega_nyq < omega < omega_nyq. ---
    end do

    call fft_time2freq_r2c_init( nsample )
    irec = 1+int(tsta/dtout_ptn)  ! do irec = 1, nrec_realsp-nsample, 1+nsample/10

!--- volume average (instead of statistical average)
    inum = 0
    do iy = 0, 2*nyw-1, 10
      do ix = 0, 2*nxw-1, 10
        inum = inum+1
        call realsp_phiixiy( ix, iy, phi )
        sample(0:nsample-1) = phi(irec:irec+nsample-1)

       !call fft_time2freq_r2c( nsample, sample, spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, gaussian_window(nsample,sample), spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, hanning_window(nsample,sample), spectrum(0:nsample/2) )
        call fft_time2freq_r2c( nsample, hamming_window(nsample,sample), spectrum(0:nsample/2) )
       !call fft_time2freq_r2c( nsample, blackman_window(nsample,sample), spectrum(0:nsample/2) )
    
        do isample = 1, nsample/2
          spectrum(-isample) = conjg( spectrum(isample) )
        end do

        do w = -nsample/2, nsample/2
          ampl(w) = ampl(w) + abs(spectrum(w))**2
        end do
        do q = -nsample/2, nsample/2
          do p = -nsample/2, nsample/2
            biampl(p,q) = biampl(p,q) + abs(spectrum(p)*spectrum(q))**2
          end do
        end do
        do q = -nsample/2, nsample/2
          do p = -nsample/2, nsample/2
            w = -p-q
            if ( -nsample/2 <= w .and. w <= nsample/2 ) then
              bispectrum(p,q) = bispectrum(p,q) + spectrum(p) * spectrum(q) * spectrum(w)
            end if
          end do
        end do

                                        !--- for debug ---
                                        ! if ( mod(inum,1000) == 0 ) then
                                        !   write(cnum,'(i8.8)') inum
                                        !   write(csta,'(i4.4)') int(tsta)
                                        !   write(cend,'(i4.4)') int(tend)
                                        !   open(omominfreq,                            &
                                        !        file="./checkbicoh_tsta"//csta//&
                                        !             "tend"//cend//"num"//cnum//".dat", &
                                        !        status="replace", action="write",      &
                                        !        form="unformatted", access="stream",   &
                                        !        convert="LITTLE_ENDIAN" )
                                        !   do q = -nsample/2, nsample/2
                                        !     do p = -nsample/2, nsample/2
                                        !       w = -p-q
                                        !       if ( -nsample/2 <= w .and.  &
                                        !                 w <= nsample/2 ) then
                                        !         bicoherence(p,q) =         &
                                        !           abs(bispectrum(p,q))**2  &
                                        !           / (biampl(p,q) * ampl(w))
                                        !       end if
                                        !       write(omominfreq) &
                                        !         real(omega(p),kind=4),  &
                                        !         real(omega(q),kind=4),  &
                                        !         real(bicoherence(p,q),kind=4)
                                        !     end do
                                        !   end do
                                        !   close( omominfreq )
                                        ! end if
    
      end do
    end do

    do q = -nsample/2, nsample/2
      do p = -nsample/2, nsample/2
        w = -p-q
        if ( -nsample/2 <= w .and. w <= nsample/2 ) then
          bicoherence(p,q) = abs(bispectrum(p,q))**2 / (biampl(p,q) * ampl(w))
        end if
      end do
    end do
    
    write(cnum,'(i8.8)') inum
    write(csta,'(i4.4)') int(tsta)
    write(cend,'(i4.4)') int(tend)
    open(omominfreq, file="./data/phiinbicoherence_tsta"//csta//"tend"//cend//"num"//cnum//".dat")
      write( omominfreq, "(a,i17,a,g17.7e3,a,g17.7e3)" )  &
                     "# nsample=",nsample, ", tsta=",tsta, ", tend=",tend
      write( omominfreq, "(99a17)" ) "#        omega(p)","omega(q)","bicoherence(p,q)"
      do q = -nsample/2, nsample/2
        do p = -nsample/2, nsample/2
          write(omominfreq,"(99g17.7e3)") omega(p), omega(q), bicoherence(p,q)
        end do
      end do
    close( omominfreq )
    open(omominfreq, file="./data/phiinfreq_tsta"//csta//"tend"//cend//"num"//cnum//".dat")
      write( omominfreq, "(a,i17,a,g17.7e3,a,g17.7e3)" )  &
                     "# nsample=",nsample, ", tsta=",tsta, ", tend=",tend
      write( omominfreq, "(99a17)" ) "#        omega","ave. spectrum"
      do w = -nsample/2, nsample/2
        write(omominfreq,"(99g17.7e3)") omega(w), ampl(w)
      end do
    close( omominfreq )
                                        !--- for debug ---
                                        !do loop=loop_phi_sta(snum),loop_phi_end(enum)
                                        !  write(*,*) time(loop),             &
                                        !             real(phi(loop)),  &
                                        !             aimag(phi(loop))
                                        !end do
                                        !-----------------
    deallocate( time )
    deallocate( phi )
    deallocate( sample )
    deallocate( omega )
    deallocate( spectrum )
    deallocate( ampl )
    deallocate( biampl )
    deallocate( bispectrum )
    deallocate( bicoherence )

    call realsp_close

END SUBROUTINE bicoh_bicoherence


END MODULE out_bicoherence
