MODULE out_zfshearing
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y)
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public zfs_open, zfs_close, zfshearing, zfshearing_timeaverage


 CONTAINS


SUBROUTINE zfs_open
    character(len=1) :: cis
    integer :: is
    open( ozfshearing, file="./data/tzfshearing.dat" )
    do is = 0, ns-1
      write( cis, '(i1.1)' ) is
      open( ozfshearing+1000*(is+1), file="./data/tzfshearing_is"//cis//".dat" )
    end do
END SUBROUTINE zfs_open
SUBROUTINE zfs_close
    integer :: is
    close( ozfshearing )
    do is = 0, ns-1
      close( ozfshearing+1000*(is+1) )
    end do
END SUBROUTINE zfs_close

SUBROUTINE zfshearing( loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : gksq, grootg, cfsrf, tau, Anum, Znum, gomg

  integer, intent(in) :: loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  complex(kind=DP), dimension(0:nx) :: omegak
  real(kind=DP) :: omega, fct, bb
  integer :: mx, my, iz, is

    call rb_phi_gettime( loop, time )
    call rb_phi_loop( loop, phi )

    omegak(:) = (0._DP, 0._DP)
    do iz = -global_nz, global_nz-1
      fct = grootg(iz) / cfsrf
      my = 0
        do mx = 0, nx
          omegak(mx) = omegak(mx) + fct * gksq(mx,my,iz) * phi(mx,my,iz)
        end do                                           ! flux-surface average
    end do
    omega = sqrt(sum(abs(omegak(:))**2))

    write( ozfshearing, "(999g17.7e3)" ) time, omega, abs(omegak(0:nx))

    do is = 0, ns-1
      omegak(:) = (0._DP, 0._DP)
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / cfsrf
        my = 0
          do mx = 0, nx
            bb = gksq(mx,my,iz)*tau(is)*Anum(is)/(Znum(is)**2*gomg(iz)**2)
            omegak(mx) = omegak(mx) + fct * gksq(mx,my,iz) * phi(mx,my,iz) * exp(-0.5_DP*bb)
          end do
      end do
      omega = sqrt(sum(abs(omegak(:)**2)))
      write( ozfshearing+1000*(is+1), "(999g17.7e3)" ) time, omega, abs(omegak(0:nx))
    end do


END SUBROUTINE zfshearing


SUBROUTINE zfshearing_timeaverage( lsta, lend, lskip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy) at giz, loop
!                                                   (S. Maeyama, 11 June 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop
  use diag_geom, only : kx, gksq, grootg, cfsrf, tau, Anum, Znum, gomg

  integer, intent(in) :: lsta, lend, lskip

  real(kind=DP) :: tsta, tend
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  complex(kind=DP), dimension(0:nx) :: omegak
  complex(kind=DP), dimension(0:nx,0:ns-1) :: omegask
  real(kind=DP), dimension(0:ns-1) :: omegas
  real(kind=DP) :: omega, fct, bb
  integer :: mx, my, iz, is, nloop, loop

    nloop = int((lend - lsta) / lskip) + 1
    call rb_phi_gettime( lsta, tsta )
    call rb_phi_gettime( lend, tend )


    omegak(:) = (0._DP, 0._DP)
    omegask(:,:) = (0._DP, 0._DP)
    do loop = lsta, lend, lskip
      call rb_phi_loop( loop, phi )
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / (cfsrf * real(nloop, kind=DP))
        my = 0
          do mx = 0, nx
            omegak(mx) = omegak(mx) + fct * gksq(mx,my,iz) * phi(mx,my,iz)
          end do                                           ! flux-surface average
      end do
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / (cfsrf * real(nloop, kind=DP))
          my = 0
            do mx = 0, nx
              bb = gksq(mx,my,iz)*tau(is)*Anum(is)/(Znum(is)**2*gomg(iz)**2)
              omegask(mx,is) = omegask(mx,is) + fct * gksq(mx,my,iz) * phi(mx,my,iz) * exp(-0.5_DP*bb)
            end do
        end do
      end do
    end do
    omega = sqrt(sum(abs(omegak(:))**2))
    do is = 0, ns-1
      omegas(is) = sqrt(sum(abs(omegask(:,is)**2)))
    end do

    open( ozfsave, file="./data/ave_zfshearing.dat" )
      write( ozfsave, '(a,i17,a,i17)' ) "# lsta=", lsta, ", lend=", lend
      write( ozfsave, '(a,i17,a,i17)' ) "#lskip=",lskip, ",nloop=",nloop
      write( ozfsave, '(a,f17.7,a,f17.7)' ) "# tsta=", tsta, ", tend=", tend
      write( ozfsave, * )
      write( ozfsave, '(a17,99g17.7e3)' ) "#           Total", omega, omegas(:)
      write( ozfsave, '(4a17)' ) "#              kx", "omega_k", "omega_ek", "omega_ik"
      do mx = 0, nx
        write( ozfsave, "(4g17.7e3)" ) kx(mx), abs(omegak(mx)), (abs(omegask(mx,is)), is = 0, ns-1)
      end do
    close( ozfsave )


END SUBROUTINE zfshearing_timeaverage


END MODULE out_zfshearing
