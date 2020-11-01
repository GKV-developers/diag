PROGRAM diag
!-------------------------------------------------------------------------------
!
!     Data diagnostics from binary output
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  use diag_geom, only : geom_init
  use diag_rb, only : loop_phi_sta, loop_phi_end
  use out_mominxy, only : phiinxy
  implicit none
  integer :: giz, loop


!- Initialization -
    call initialize
    call geom_init( snum )


!- output -
    giz = 0
    loop = 10
    call phiinxy( giz, loop )


!- Finalization -
    call finalize

END PROGRAM diag



















SUBROUTINE initialize
!-------------------------------------------------------------------------------
!
!     Initialize diag
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_clock, only : clock_init
  use diag_rb, only : rb_phi_fileopen, rb_phi_check, rb_phi_setnloop,  &
                      rb_Al_fileopen,  rb_Al_check,  rb_Al_setnloop,   &
                      rb_mom_fileopen, rb_mom_check, rb_mom_setnloop,  &
                      rb_trn_fileopen, rb_trn_check, rb_trn_setnloop,  &
                      rb_tri_fileopen, rb_tri_check, rb_tri_setnloop,  &
                      rb_cnt_fileopen, rb_cnt_check, rb_cnt_setnloop
  use diag_fft, only : fft_pre
  implicit none

    call clock_init
    call rb_phi_fileopen
    call rb_phi_check
    call rb_phi_setnloop
    call rb_Al_fileopen
    call rb_Al_check
    call rb_Al_setnloop
    call rb_mom_fileopen
    call rb_mom_check
    call rb_mom_setnloop
    call rb_trn_fileopen
    call rb_trn_check
    call rb_trn_setnloop
    call rb_tri_fileopen
    call rb_tri_check
    call rb_tri_setnloop
    !call rb_cnt_fileopen
    !call rb_cnt_check
    !call rb_cnt_setnloop
    call fft_pre

END SUBROUTINE initialize


SUBROUTINE finalize
!-------------------------------------------------------------------------------
!
!     Finalize diag
!                                                   (S. Maeyama, 9 Dec. 2016)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_fileclose,  &
                      rb_Al_fileclose,   &
                      rb_mom_fileclose,  &
                      rb_trn_fileclose,  &
                      rb_tri_fileclose,  &
                      rb_cnt_fileclose
  implicit none

    call rb_phi_fileclose
    call rb_Al_fileclose
    call rb_mom_fileclose
    call rb_trn_fileclose
    call rb_tri_fileclose
    !call rb_cnt_fileclose

END SUBROUTINE finalize
