MODULE diag_clock
!-------------------------------------------------------------------------------
!
!    clock module
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header, only : DP
  implicit none

  private

  public   clock_init, clock_sta, clock_end, elt

  real(kind=DP), dimension(0:100), save :: elt
  integer, dimension(0:100), save :: s_count, e_count


 CONTAINS


SUBROUTINE clock_init
!-------------------------------------------------------------------------------
!
!    Initialize clock
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------

    elt(:) = 0._DP
    s_count(:) = 0
    e_count(:) = 0

END SUBROUTINE clock_init


SUBROUTINE clock_sta( id )
!-------------------------------------------------------------------------------
!
!    Get elapsed time
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------

  integer, intent(in) :: id
  !integer :: count_per_sec, count_max

    call system_clock(s_count(id))
    !call system_clock(s_count(id), count_per_sec, count_max)

END SUBROUTINE clock_sta


SUBROUTINE clock_end( id )
!-------------------------------------------------------------------------------
!
!    Get elapsed time
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------

  integer, intent(in) :: id
  integer :: count_per_sec, count_max

    call system_clock(e_count(id), count_per_sec, count_max)

    if ( e_count(id) < s_count(id) ) then
      elt(id) = elt(id)  &
              + real(count_max + e_count(id) - s_count(id), kind=DP)  &
              / real(count_per_sec, kind=DP)
    else
      elt(id) = elt(id)  &
              + real(e_count(id) - s_count(id), kind=DP)  &
              / real(count_per_sec, kind=DP)
    end if

END SUBROUTINE clock_end


END MODULE diag_clock

