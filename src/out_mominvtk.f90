MODULE out_mominvtk
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y,z)
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinvtk

  integer, parameter :: n_alp = 6
  ! Adapt a flux tube as 1/n_alp torus for visualization.
  ! Then, Larmor radius rho/r_major = pi*eps_r/(q_0*ly*n_alp)

  integer, parameter :: nzw = 3 * global_nz
  ! Linear interpolation increases field-aligned resolution, nzw >= global_nz.

 CONTAINS


SUBROUTINE phiinvtk( flag, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in VTK binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly
  use diag_rb, only : rb_phi_loop
  implicit none
  integer, intent(in) :: flag ! 1:coord&var&header(fluxtube)
                              ! 3:coord&var&header(fulltorus)
                              ! 5:coord&var&header(fieldalignedcoord)
  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP), dimension(:,:,:,:), allocatable :: coords
  real(kind=DP), dimension(:,:,:), allocatable :: phi_xyz
  complex(kind=DP), dimension(:,:,:), allocatable :: phi
  character(len=8) :: cloop
  character(len=3) :: c_alp
  integer(kind=8) :: offset, head ! UInt64
  integer :: loop, i_alp

    allocate(coords(3,0:2*nxw,0:2*nyw,-nzw:nzw))
    allocate(phi_xyz(0:2*nxw,0:2*nyw,-nzw:nzw))
    allocate(phi(-nx:nx,0:global_ny,-global_nz:global_nz-1))

    if ( flag == 1 ) then ! output variable:fluxtube

      != set coordinates =
      i_alp = 0
      call cartesian_coordinates( i_alp, coords )

      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop

        != set values =
        call rb_phi_loop( loop, phi )
        call phi_kxkyz2xyz_fluxtube( phi, phi_xyz )

        != output VTK-structured-grid file =
        open( omominxyz, file="./data/phiinvtk_tube_t"//cloop//".vts", status="replace", &
                         action="write", form="formatted", access="sequential" )
          write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
          write( omominxyz, '(a)' ) '<VTKFile type="StructuredGrid" version="1.0" '//&
                                     'byte_order="LittleEndian" header_type="UInt64">'
          write( omominxyz, '(a,i0,a,i0,a,i0,a)' ) '<StructuredGrid WholeExtent="0 ',2*nxw,' 0 ',2*nyw,' 0 ',2*nzw,'">'
          write( omominxyz, '(a,i0,a,i0,a,i0,a)' ) '<Piece Extent="0 ',2*nxw,' 0 ',2*nyw,' 0 ',2*nzw,'">'
          write( omominxyz, '(a)' ) '<Points>'
          write( omominxyz, '(a)' ) '<DataArray Name="coordinates" NumberOfComponents="3" '//&
                                     'type="Float32" format="appended" offset="0"/>'
          write( omominxyz, '(a)' ) '</Points>'
          write( omominxyz, '(a)' ) '<PointData scalars="phi">'
          offset = 8 + 3*(2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! head(Uint64)+coordinates(Float32)
          write( omominxyz, '(a,i0,a)' ) '<DataArray Name="phi" NumberOfComponents="1" type="Float32" '//&
                                          'format="appended" offset="',offset,'"/>'
          write( omominxyz, '(a)' ) '</PointData>'
          write( omominxyz, '(a)' ) '</Piece>'
          write( omominxyz, '(a)' ) '</StructuredGrid>'
          write( omominxyz, '(a)' ) '<AppendedData encoding="raw">'
        close( omominxyz )

        open( omominxyz, file="./data/phiinvtk_tube_t"//cloop//".vts", status="old", &
                         action="write", form="unformatted", access="stream", position="append", convert="LITTLE_ENDIAN" )
          write( omominxyz ) '_' ! Literal underscore is required.
          head = 3*(2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! coordinates(Float32)
          write( omominxyz ) head, real(coords, kind=4)
          head = (2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! phi(Float32)
          write( omominxyz ) head, real(phi_xyz, kind=4)
        close( omominxyz )

        open( omominxyz, file="./data/phiinvtk_tube_t"//cloop//".vts", status="old", &
                         action="write", form="formatted", access="sequential", position="append" )
          write( omominxyz, '(a)' ) ''
          write( omominxyz, '(a)' ) '</AppendedData>'
          write( omominxyz, '(a)' ) '</VTKFile>'
        close( omominxyz )

      end do

    else if ( flag == 3 ) then ! output variable:full torus

      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop

        != set values =
        call rb_phi_loop( loop, phi )
        call phi_kxkyz2xyz_fluxtube( phi, phi_xyz )

        do i_alp = 0, n_alp-1
          write( c_alp, '(i3.3)' ) i_alp

          != set coordinates =
          call cartesian_coordinates( i_alp, coords )

          != output VTK-structured-grid file =
          open( omominxyz, file="./data/phiinvtk_full_t"//cloop//"_alp"//c_alp//".vts", status="replace", &
                           action="write", form="formatted", access="sequential" )
            write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
            write( omominxyz, '(a)' ) '<VTKFile type="StructuredGrid" version="1.0" '//&
                                       'byte_order="LittleEndian" header_type="UInt64">'
            write( omominxyz, '(a,i0,a,i0,a,i0,a,i0,a)' ) '<StructuredGrid WholeExtent="0 ',2*nxw,' ',&
                                                          2*nyw*i_alp,' ',2*nyw*(i_alp+1),' 0 ',2*nzw,'">'
            write( omominxyz, '(a,i0,a,i0,a,i0,a,i0,a)' ) '<Piece Extent="0 ',2*nxw,' ',&
                                                          2*nyw*i_alp,' ',2*nyw*(i_alp+1),' 0 ',2*nzw,'">'
            write( omominxyz, '(a)' ) '<Points>'
            write( omominxyz, '(a)' ) '<DataArray Name="coordinates" NumberOfComponents="3" '//&
                                       'type="Float32" format="appended" offset="0"/>'
            write( omominxyz, '(a)' ) '</Points>'
            write( omominxyz, '(a)' ) '<PointData scalars="phi">'
            offset = 8 + 3*(2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! head(Uint64)+coordinates(Float32)
            write( omominxyz, '(a,i0,a)' ) '<DataArray Name="phi" NumberOfComponents="1" type="Float32" '//&
                                            'format="appended" offset="',offset,'"/>'
            write( omominxyz, '(a)' ) '</PointData>'
            write( omominxyz, '(a)' ) '</Piece>'
            write( omominxyz, '(a)' ) '</StructuredGrid>'
            write( omominxyz, '(a)' ) '<AppendedData encoding="raw">'
          close( omominxyz )

          open( omominxyz, file="./data/phiinvtk_full_t"//cloop//"_alp"//c_alp//".vts", status="old", &
                           action="write", form="unformatted", access="stream", position="append", convert="LITTLE_ENDIAN" )
            write( omominxyz ) '_' ! Literal underscore is required.
            head = 3*(2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! coordinates(Float32)
            write( omominxyz ) head, real(coords, kind=4)
            head = (2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! phi(Float32)
            write( omominxyz ) head, real(phi_xyz, kind=4)
          close( omominxyz )

          open( omominxyz, file="./data/phiinvtk_full_t"//cloop//"_alp"//c_alp//".vts", status="old", &
                           action="write", form="formatted", access="sequential", position="append" )
            write( omominxyz, '(a)' ) ''
            write( omominxyz, '(a)' ) '</AppendedData>'
            write( omominxyz, '(a)' ) '</VTKFile>'
          close( omominxyz )

        end do

        open( omominxyz, file="./data/phiinvtk_full_t"//cloop//".pvts", status="replace", &
                         action="write", form="formatted", access="sequential" )
          write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
          write( omominxyz, '(a)' ) '<VTKFile type="PStructuredGrid" version="1.0" '//&
                                     'byte_order="LittleEndian" header_type="UInt64">'
          write( omominxyz, '(a,i0,a,i0,a,i0,a,i0,a)' ) '<PStructuredGrid WholeExtent="0 ',2*nxw,' ',&
                                                        0,' ',2*nyw*n_alp,' 0 ',2*nzw,'" GhostLevel="#">'
          write( omominxyz, '(a)' ) '<PPoints>'
          write( omominxyz, '(a)' ) '<PDataArray Name="coordinates" NumberOfComponents="3" type="Float32"/>'
          write( omominxyz, '(a)' ) '</PPoints>'
          write( omominxyz, '(a)' ) '<PPointData scalars="phi">'
          write( omominxyz, '(a,i0,a)' ) '<PDataArray Name="phi" NumberOfComponents="1" type="Float32"/>'
          write( omominxyz, '(a)' ) '</PPointData>'
          do i_alp = 0, n_alp-1
            write( c_alp, '(i3.3)' ) i_alp
            write( omominxyz, '(a,i0,a,i0,a,i0,a,i0,a)' )  &
                                    '<Piece Extent="0 ',2*nxw,' ',2*nyw*i_alp,' ',2*nyw*(i_alp+1),' 0 ',2*nzw,'" '//&
                                     'Source="./phiinvtk_full_t'//cloop//'_alp'//c_alp//'.vts"/>'
          end do
          write( omominxyz, '(a)' ) '</PStructuredGrid>'
          write( omominxyz, '(a)' ) '</VTKFile>'
        close( omominxyz )

      end do

    else if ( flag == 5 ) then ! output variable: Field aligned

      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop

        != set values =
        call rb_phi_loop( loop, phi )
        call phi_kxkyz2xyz_fluxtube( phi, phi_xyz )

        != output VTK-structured-grid file =
        open( omominxyz, file="./data/phiinvtk_align_t"//cloop//".vti", status="replace", &
                         action="write", form="formatted", access="sequential" )
          write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
          write( omominxyz, '(a)' ) '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
          write( omominxyz, '(a,i0,a,i0,a,i0,a,g17.7,a,g17.7,a,g17.7,a)' ) &
                                    '<ImageData WholeExtent="0 ',2*nxw,' 0 ',2*nyw,' 0 ',2*nzw,'" '//&
                                    'Origin="0 0 0" Spacing="',lx/real(nxw),' ',ly/real(nyw),' ',4._DP*ly/real(nzw),'">'
          write( omominxyz, '(a,i0,a,i0,a,i0,a)' ) '<Piece Extent="0 ',2*nxw,' 0 ',2*nyw,' 0 ',2*nzw,'">'
          write( omominxyz, '(a)' ) '<PointData scalars="phi">'
          write( omominxyz, '(a,i0,a)' ) '<DataArray Name="phi" NumberOfComponents="1" type="Float32" '//&
                                        'format="appended" offset="0"/>'
          write( omominxyz, '(a)' ) '</PointData>'
          write( omominxyz, '(a)' ) '</Piece>'
          write( omominxyz, '(a)' ) '</ImageData>'
          write( omominxyz, '(a)' ) '<AppendedData encoding="raw">'
        close( omominxyz )

        open( omominxyz, file="./data/phiinvtk_align_t"//cloop//".vti", status="old", &
                         action="write", form="unformatted", access="stream", position="append", convert="LITTLE_ENDIAN" )
          write( omominxyz ) '_' ! Literal underscore is required.
          head = (2*nxw+1)*(2*nyw+1)*(2*nzw+1)*4 ! phi(Float32)
          write( omominxyz ) head, real(phi_xyz, kind=4)
        close( omominxyz )

        open( omominxyz, file="./data/phiinvtk_align_t"//cloop//".vti", status="old", &
                         action="write", form="formatted", access="sequential", position="append" )
          write( omominxyz, '(a)' ) ''
          write( omominxyz, '(a)' ) '</AppendedData>'
          write( omominxyz, '(a)' ) '</VTKFile>'
        close( omominxyz )

      end do

    else

      write(*,*) "flag in phiinvtk is invalid !  flag =", flag
      stop

    end if

END SUBROUTINE phiinvtk


SUBROUTINE phi_kxkyz2xyz_fluxtube( phi, phi_xyz )
!-------------------------------------------------------------------------------
!
!    Calculate real space quantity
!
!-------------------------------------------------------------------------------
  use diag_geom, only : ck, dj, lz, dz, gzz
  use diag_fft, only : fft_backward_xy
  implicit none
  complex(kind=DP), intent(in), &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  real(kind=DP), intent(out), &
    dimension(0:2*nxw,0:2*nyw,-nzw:nzw) :: phi_xyz

  complex(kind=DP), dimension(:,:,:), allocatable :: phi_kxkyz
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: phi_interp
  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wr2
  real(kind=DP) :: wdz, wzz, alpha
  integer :: mx, my, iz, izw, mwp

    allocate(phi_kxkyz(-nx:nx,0:global_ny,-global_nz:global_nz))
  
    !- pseudo-periodic boundary in z -
    phi_kxkyz(:,:,-global_nz:global_nz-1) = phi(:,:,:)
    iz = global_nz
      do my = 0, global_ny
        do mx = -nx, nx
          mwp = mx - dj(my)          ! --- mw = mx - dj for the positive-z 
            if( abs(mwp) > nx ) then
              phi_kxkyz(mx,my,iz) = ( 0._DP, 0._DP )
            else
              phi_kxkyz(mx,my,iz) = conjg( ck(my) ) * phi_kxkyz(mwp,my,-global_nz)
            end if
        end do
      end do

    wdz = lz / real( nzw, kind=DP )
!$OMP parallel do default(none) &
!$OMP shared(wdz,dz,gzz,phi_kxkyz,phi_xyz) &
!$OMP private(iz,izw,my,mx,wzz,alpha,phi_interp,wr2)
    do izw = -nzw, nzw

      !- interpolate along z -
      wzz = wdz * real(izw, kind=DP)
      if (izw == -nzw) then
        phi_interp(:,:) = phi_kxkyz(:,:,-global_nz)
      else if (izw == nzw) then
        phi_interp(:,:) = phi_kxkyz(:,:,global_nz)
      else
        iz = int(wzz / dz) ! find position zz(iz)<= zzw < zz(iz+1)
        alpha = (wzz - gzz(iz)) / dz
        phi_interp(:,:) = (1._DP - alpha) * phi_kxkyz(:,:,iz) + alpha * phi_kxkyz(:,:,iz+1)
      end if

      !- backward fft in x,y -
      call fft_backward_xy( phi_interp, wr2 )
      phi_xyz(0:2*nxw-1,0:2*nyw-1,izw) =  wr2(:,:)
      phi_xyz(2*nxw,0:2*nyw-1,izw) = phi_xyz(0,0:2*nyw-1,izw)
      phi_xyz(0:2*nxw,2*nyw,izw) = phi_xyz(0:2*nxw,0,izw)

    end do

    deallocate(phi_kxkyz)


END SUBROUTINE phi_kxkyz2xyz_fluxtube


SUBROUTINE cartesian_coordinates( i_alp, coords )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : q_0, s_hat, eps_r
  implicit none
  integer, intent(in) :: i_alp
  real(kind=DP), intent(out), &
    dimension(3,0:2*nxw,0:2*nyw,-nzw:nzw) :: coords

    call cartesian_coordinates_salpha(coords,i_alp,q_0,s_hat,eps_r)
    !call cartesian_coordinates_miller(coords,i_alp,q_0,s_hat,eps_r,&
    !       dRmildr=-0.1d0,dZmildr=0.0,kappa=1.5d0,s_kappa=0.7d0,delta=0.4d0,s_delta=1.3d0,zetasq=0.d0,s_zetasq=0.d0)

END SUBROUTINE cartesian_coordinates


SUBROUTINE cartesian_coordinates_salpha( coords, i_alp, q_0, s_hat, eps_r )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz
  implicit none
  real(kind=DP), intent(out), &
    dimension(3,0:2*nxw,0:2*nyw,-nzw:nzw) :: coords
  integer, intent(in) :: i_alp
  real(kind=DP), intent(in) :: q_0, s_hat, eps_r
  real(kind=DP) :: rho, q_r
  real(kind=DP) :: wxx, wyy, wzz, wsr, wmr, wdx, wdy, wdz, wtheta, wzeta, wx_car, wy_car, wz_car
  integer :: mx, my, iz

    rho = pi * eps_r / (q_0 * ly * n_alp) ! = Larmor radius rho/r_major
    write(*,*) "# Larmor radius rho/r_major = ", rho
    if (lx*rho > eps_r) then
      write(*,*) "# WARNING in out_mominvtk. lx*rho < eps_r is recommended. Set larger n_alp."
      write(*,*) "# lx=",lx,", rho=",rho,", eps_r=",eps_r,", n_alp=",n_alp 
    end if

    wdx = lx / real( nxw, kind=DP )
    wdy = ly / real( nyw, kind=DP )
    wdz = lz / real( nzw, kind=DP )
!$OMP parallel do default(none) &
!$OMP shared(i_alp,wdx,wdy,wdz,lx,ly,q_0,s_hat,rho,eps_r,coords) &
!$OMP private(iz,my,mx,wzz,wtheta,wyy,wxx,q_r,wsr,wmr,wzeta,wx_car,wy_car,wz_car)
    do iz = -nzw, nzw
      wzz = wdz * real( iz, kind=DP )
      wtheta = wzz
      do my = 0, 2*nyw
        wyy = -ly + wdy * real( my, kind=DP )
        do mx = 0, 2*nxw
          wxx = -lx + wdx * real( mx, kind=DP )
          q_r = q_0 * ( 1._DP + s_hat * wxx * rho / eps_r )

          wsr = eps_r + rho * wxx
          wmr = 1._DP + wsr * cos( wtheta )
          wzeta = q_r * wtheta - q_0 * wyy * rho / eps_r  &
                - 2._DP * pi * real(i_alp, kind=DP) / real(n_alp, kind=DP)
          wx_car = wmr * cos( wzeta )
          wy_car = wmr * sin( wzeta )
          wz_car = wsr * sin( wtheta )
          coords(1,mx,my,iz) = wx_car
          coords(2,mx,my,iz) = wy_car
          coords(3,mx,my,iz) = wz_car
        end do
      end do
    end do

END SUBROUTINE cartesian_coordinates_salpha


SUBROUTINE cartesian_coordinates_miller( coords, i_alp, q_0, s_hat, eps_r, &
             dRmildr, dZmildr, kappa, s_kappa, delta, s_delta, zetasq, s_zetasq )
!-------------------------------------------------------------------------------
!
!    Calculate Cartesian coordinates for s-alpha model
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz
  implicit none
  real(kind=DP), intent(out), &
    dimension(3,0:2*nxw,0:2*nyw,-nzw:nzw) :: coords
  integer, intent(in) :: i_alp
  real(kind=DP), intent(in) :: q_0, s_hat, eps_r, dRmildr, dZmildr, kappa, s_kappa, delta, s_delta, zetasq, s_zetasq
  real(kind=DP) :: rho, q_r, kappa_r, delta_r, zetasq_r, Rmil_r, Zmil_r
  real(kind=DP) :: wxx, wyy, wzz, wsr, wmr, wdx, wdy, wdz, wtheta, wzeta, wx_car, wy_car, wz_car
  integer :: mx, my, iz

    rho = pi * eps_r / (q_0 * ly * n_alp) ! = Larmor radius rho/r_major
    write(*,*) "# Larmor radius rho/r_major = ", rho
    if (lx*rho > eps_r) then
      write(*,*) "# WARNING in out_mominvtk. lx*rho < eps_r is recommended. Set larger n_alp."
      write(*,*) "# lx=",lx,", rho=",rho,", eps_r=",eps_r,", n_alp=",n_alp 
    end if

    wdx = lx / real( nxw, kind=DP )
    wdy = ly / real( nyw, kind=DP )
    wdz = lz / real( nzw, kind=DP )
    do iz = -nzw, nzw
      wzz = wdz * real( iz, kind=DP )
      wtheta = wzz
      do my = 0, 2*nyw
        wyy = -ly + wdy * real( my, kind=DP )
        do mx = 0, 2*nxw
          wxx = -lx + wdx * real( mx, kind=DP )
          q_r = q_0 * ( 1._DP + s_hat * wxx * rho / eps_r )
          kappa_r = kappa * ( 1._DP + s_kappa * wxx * rho / eps_r )
          delta_r = delta + sqrt( 1._DP - delta**2) * s_delta * wxx * rho / eps_r
          zetasq_r = zetasq + s_zetasq * wxx * rho / eps_r
          Rmil_r = 1.d0 + dRmildr * wxx * rho
          Zmil_r = 0.d0 + dZmildr * wxx * rho

          wsr = eps_r + rho * wxx
          wmr = Rmil_r + wsr * cos( wtheta + asin(delta_r) * sin(wtheta) )
          wzeta = q_r * wtheta - q_0 * wyy * rho / eps_r  &
                - 2._DP * pi * real(i_alp, kind=DP) / real(n_alp, kind=DP)
          wx_car = wmr * cos( wzeta )
          wy_car = wmr * sin( wzeta )
          wz_car = Zmil_r + wsr * kappa_r * sin( wtheta + zetasq_r * sin(2._DP*wtheta) )
          coords(1,mx,my,iz) = wx_car
          coords(2,mx,my,iz) = wy_car
          coords(3,mx,my,iz) = wz_car
        end do
      end do
    end do

END SUBROUTINE cartesian_coordinates_miller


END MODULE out_mominvtk
