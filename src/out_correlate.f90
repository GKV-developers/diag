MODULE out_correlate
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 13 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public correlate


 CONTAINS


SUBROUTINE correlate( giz, loop_sta, loop_end )
!-------------------------------------------------------------------------------
!
!     Output correlation of moments in (kx,ky)
!                                                   (S. Maeyama, 13 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_izloop
  use diag_geom, only : kx, gky, gzz

  integer, intent(in) :: giz, loop_sta, loop_end

  real(kind=DP) :: time_sta, time_end
  complex(kind=DP), &
    dimension(-nx:nx,0:global_ny) :: gphi
  complex(kind=DP), &
    dimension(-nx:nx,-global_ny:global_ny) :: wkphi
!  complex(kind=DP), &
!    dimension(-nx:nx,-global_ny:global_ny,-nx:nx,-global_ny:global_ny) :: corr
  complex(kind=DP), &
    dimension(-nx:nx,-global_ny:global_ny) :: intcorr
  real(kind=DP), &
    dimension(-nx:nx,-global_ny:global_ny) :: norm1, norm2
  real(kind=DP), &
    dimension(-global_ny:global_ny) :: ky
  character(len=4) :: ciz
  character(len=8) :: cloop
  integer :: nsample
  integer :: mx, my, px, py, qx, qy, loop

    nsample = loop_end - loop_sta + 1
    call rb_phi_gettime( loop_sta, time_sta )
    call rb_phi_gettime( loop_end, time_end )

    intcorr(:,:) = ( 0._DP, 0._DP )
    norm1(:,:) = 0._DP
    norm2(:,:) = 0._DP

    do my = 0, global_ny
      ky(-my) = -gky(my)
      ky(my) = gky(my)
    end do

    do loop = loop_sta, loop_end
      call rb_phi_izloop( giz, loop, gphi )
!$OMP parallel do
      do my = 0, global_ny
        do mx = -nx, nx
          wkphi(mx,my) = gphi(mx,my)
        end do
      end do
!$OMP parallel do
      do my = 1, global_ny
        do mx = -nx, nx
          wkphi(-mx,-my) = conjg( gphi(mx,my) )
        end do
      end do

!$OMP parallel do private(qx,qy) reduction(+:intcorr)
      do my = -global_ny, global_ny
        do mx = -nx, nx
          !NOTE: The other corr(px,py) should be zero.
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - py - my
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - px - mx
              !corr(px,py,mx,my) = wkphi(mx,my)*wkphi(px,py)*wkphi(qx,qy)
              intcorr(mx,my) = intcorr(mx,my)  &
                + wkphi(mx,my)*conjg(wkphi(px,py))*wkphi(qx,qy)
                !+ (kx(px)*ky(qy) - ky(py)*kx(qx))*wkphi(mx,my)*wkphi(px,py)*wkphi(qx,qy)
              norm1(mx,my) = norm1(mx,my) + abs(wkphi(mx,my))**2
              norm2(mx,my) = norm2(mx,my) + abs(wkphi(px,py)*wkphi(qx,qy))**2
            end do
          end do
          !--- same as follow ---
          !do py = -global_ny, global_ny
          !  qy = - py - my
          !  do px = -nx, nx  ! px = max(-nx,-nx-mx), min(nx,nx-mx) like that?
          !    qx = - px - mx
          !    if ( abs(qy) <= global_ny .and. abs(qx) <= nx ) then
          !      intcorr(mx,my) = intcorr(mx,my) + wkphi(mx,my)*wkphi(px,py)*wkphi(qx,qy)
          !    else
          !      intcorr(mx,my) = intcorr(mx,my) + ( 0._DP, 0._DP )
          !    end if
          !  end do
          !end do
          !----------------------
        end do
      end do

    end do

!$OMP parallel do
    do my = -global_ny, global_ny
      do mx = -nx, nx
        !intcorr(mx,my) = intcorr(mx,my) / real( nsample, kind=DP )
        intcorr(mx,my) = abs(intcorr(mx,my))**2 / (norm1(mx,my) * norm2(mx,my))
      end do
    end do

    write( ciz, fmt="(i4.4)" ) giz
    write( cloop, fmt="(i8.8)" ) loop_sta
    open( ocorrelate, file="./data/correlate_z"//ciz//"_t"//cloop//".dat" )
      write( ocorrelate, "(a,i17,a,g17.7e3)" )  &
                         "# loop_sta=",loop_sta, ", time_sta=",time_sta
      write( ocorrelate, "(a,i17,a,g17.7e3)" )  &
                         "# loop_end=",loop_end, ", time_end=",time_end
      write( ocorrelate, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",   zz=",gzz(giz)
      write( ocorrelate, "(99a17)" ) "#              kx","ky","corr"
      do my = -global_ny, -1
        do mx = -nx, nx
          write( ocorrelate, "(99g17.7e3)" )kx(mx),-gky(-my),abs(intcorr(mx,my))
        end do
        write( ocorrelate, * )
      end do
      do my = 0, global_ny
        do mx = -nx, nx
          write( ocorrelate, "(99g17.7e3)" ) kx(mx), gky(my),abs(intcorr(mx,my))
        end do
        write( ocorrelate, * )
      end do
    close( ocorrelate )

END SUBROUTINE correlate


!SUBROUTINE thet_ave_z ( wn, wa )
!!-------------------------------------------------------------------------------
!!
!!     Average of a complex variable wn in the theta space
!!                                                   (S. Maeyama, 13 March 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_geom, only : grootg, cfsrf
!  complex(kind=DP), intent(in),  &
!    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wn
!  complex(kind=DP), intent(out), &
!    dimension(-nx:nx,0:global_ny)                        :: wa
!
!  real(kind=DP) :: fct
!  integer ::  mx, my, iz
!
!    wa   = ( 0._DP, 0._DP )
!    do iz = -global_nz, global_nz-1
!      fct = grootg(iz) / cfsrf
!      do my = 0, global_ny
!        do mx = -nx, nx
!          wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
!        end do
!      end do
!    end do
!
!END SUBROUTINE thet_ave_z
!
!
!SUBROUTINE thet_ave_r ( wn, wa )
!!-------------------------------------------------------------------------------
!!
!!     Average of a real variable wn in the theta space
!!                                                   (S. Maeyama, 13 March 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_geom, only : grootg, cfsrf
!  real(kind=DP), intent(in),  &
!    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wn
!  real(kind=DP), intent(out), &
!    dimension(-nx:nx,0:global_ny)                        :: wa
!  
!  real(kind=DP) :: fct
!  integer ::  mx, my, iz
!  
!    wa   = 0._DP
!    do iz = -global_nz, global_nz-1
!      fct = grootg(iz) / cfsrf
!      do my = 0, global_ny
!        do mx = -nx, nx
!          wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
!        end do
!      end do
!    end do
!
!END SUBROUTINE thet_ave_r


END MODULE out_correlate
