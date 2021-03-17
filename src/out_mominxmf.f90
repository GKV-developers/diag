MODULE out_mominxmf
!-------------------------------------------------------------------------------
!
!     Output moments in (x,y,z)
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public  mominxmf_coord, &
          mominxmf_var_phi, mominxmf_header_phi, &
          mominxmf_var_Al,  mominxmf_header_Al,  &
          mominxmf_var_mom, mominxmf_header_mom

  integer, parameter :: n_alp = 6
  ! Adapt a flux tube as 1/n_alp torus for visualization.
  ! Then, Larmor radius rho/r_major = pi*eps_r/(q_0*ly*n_alp)

  integer, parameter :: nzw = 3 * global_nz
  ! Linear interpolation increases field-aligned resolution, nzw >= global_nz.

 CONTAINS


SUBROUTINE mominxmf_coord
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  real(kind=DP) :: time
  real(kind=DP), dimension(:,:,:,:), allocatable :: coords
  character(len=3) :: c_alp
  integer :: i_alp

    allocate(coords(3,0:2*nxw,0:2*nyw,-nzw:nzw))

      do i_alp = 0, n_alp-1
        write( c_alp, '(i3.3)' ) i_alp

        != set coordinates =
        call cartesian_coordinates( i_alp, coords )

        open( omominxyz, file='./data/mominxmf_alp'//c_alp//'_xcoord.bin', status='replace',  &
                         action='write', form='unformatted', access='stream',       &
                         convert='LITTLE_ENDIAN' )
          write( omominxyz ) real(coords(1,:,:,:), kind=4)
        close( omominxyz )
        open( omominxyz, file='./data/mominxmf_alp'//c_alp//'_ycoord.bin', status='replace',  &
                         action='write', form='unformatted', access='stream',       &
                         convert='LITTLE_ENDIAN' )
          write( omominxyz ) real(coords(2,:,:,:), kind=4)
        close( omominxyz )
        open( omominxyz, file='./data/mominxmf_alp'//c_alp//'_zcoord.bin', status='replace',  &
                         action='write', form='unformatted', access='stream',       &
                         convert='LITTLE_ENDIAN' )
          write( omominxyz ) real(coords(3,:,:,:), kind=4)
        close( omominxyz )

      end do

    deallocate(coords)

END SUBROUTINE mominxmf_coord


SUBROUTINE mominxmf_var_phi(loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_loop

  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  real(kind=DP), dimension(:,:,:), allocatable :: phi_xyz
  complex(kind=DP), dimension(:,:,:), allocatable :: phi
  character(len=8) :: cloop
  integer :: loop

    allocate(phi_xyz(0:2*nxw,0:2*nyw,-nzw:nzw))
    allocate(phi(-nx:nx,0:global_ny,-global_nz:global_nz-1))

    do loop = loop_sta, loop_end, loop_skip
      write( cloop, '(i8.8)' ) loop
      call rb_phi_loop( loop, phi )
      call phi_kxkyz2xyz_fluxtube( phi, phi_xyz )
      open( omominxyz, file='./data/phiinxmf_var'//cloop//'.bin', status='replace',  &
                       action='write', form='unformatted', access='stream',               &
                       convert='LITTLE_ENDIAN' )
        write( omominxyz ) real(phi_xyz, kind=4)
      close( omominxyz )
    end do

    deallocate(phi_xyz)
    deallocate(phi)

END SUBROUTINE mominxmf_var_phi


SUBROUTINE mominxmf_header_phi( loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz
  use diag_rb, only : rb_phi_gettime

  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  character(len=8) :: cloop
  character(len=3) :: c_alp
  integer :: loop, i_alp

    open( omominxyz, file='./data/phiinxmf_header_align.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="phi" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DCORECTMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="ORIGIN_DXDYDZ">'
        write( omominxyz, '(a)' ) '<DataItem Name="Origin" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) 0.d0, 0.d0, 0.d0
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '<DataItem Name="DxDyDz" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) max(2*ly,2*lx)/(2*nzw), (2*ly)/(2*nyw), (2*lx)/(2*nxw)
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="phi">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  phiinxmf_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/phiinxmf_header_tube.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="phi" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_xcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_ycoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_zcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="phi">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  phiinxmf_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/phiinxmf_header_full.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="phi" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        do i_alp = 0, n_alp-1
          write( c_alp, '(i3.3)' ) i_alp
          write( omominxyz, '(a)' ) '<Grid Name="phi_var'//cloop//'_alp'//c_alp//'">'
          write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
          write( omominxyz, '(a)' ) '</Topology>'
          write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_xcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_ycoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_zcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Geometry>'
          write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="phi">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  phiinxmf_var'//cloop//'.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Attribute>'
          write( omominxyz, '(a)' ) '</Grid>'
        end do
        write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

END SUBROUTINE mominxmf_header_phi


SUBROUTINE mominxmf_var_Al(loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output Al in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_loop

  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  real(kind=DP), dimension(:,:,:), allocatable :: Al_xyz
  complex(kind=DP), dimension(:,:,:), allocatable :: Al
  character(len=8) :: cloop
  integer :: loop

    allocate(Al_xyz(0:2*nxw,0:2*nyw,-nzw:nzw))
    allocate(Al(-nx:nx,0:global_ny,-global_nz:global_nz-1))

    do loop = loop_sta, loop_end, loop_skip
      write( cloop, '(i8.8)' ) loop
      call rb_Al_loop( loop, Al )
      call phi_kxkyz2xyz_fluxtube( Al, Al_xyz )
      open( omominxyz, file='./data/Alinxmf_var'//cloop//'.bin', status='replace',  &
                       action='write', form='unformatted', access='stream',               &
                       convert='LITTLE_ENDIAN' )
        write( omominxyz ) real(Al_xyz, kind=4)
      close( omominxyz )
    end do

    deallocate(Al_xyz)
    deallocate(Al)

END SUBROUTINE mominxmf_var_Al


SUBROUTINE mominxmf_header_Al( loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output Al in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz
  use diag_rb, only : rb_Al_gettime

  integer, intent(in) :: loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  character(len=8) :: cloop
  character(len=3) :: c_alp
  integer :: loop, i_alp

    open( omominxyz, file='./data/Alinxmf_header_align.xmf' )
      call rb_Al_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="Al" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DCORECTMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="ORIGIN_DXDYDZ">'
        write( omominxyz, '(a)' ) '<DataItem Name="Origin" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) 0.d0, 0.d0, 0.d0
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '<DataItem Name="DxDyDz" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) max(2*ly,2*lx)/(2*nzw), (2*ly)/(2*nyw), (2*lx)/(2*nxw)
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="Al">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  Alinxmf_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/Alinxmf_header_tube.xmf' )
      call rb_Al_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="Al" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_xcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_ycoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_zcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="Al">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  Alinxmf_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/Alinxmf_header_full.xmf' )
      call rb_Al_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="Al" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        do i_alp = 0, n_alp-1
          write( c_alp, '(i3.3)' ) i_alp
          write( omominxyz, '(a)' ) '<Grid Name="Al_var'//cloop//'_alp'//c_alp//'">'
          write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
          write( omominxyz, '(a)' ) '</Topology>'
          write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_xcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_ycoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_zcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Geometry>'
          write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="Al">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  Alinxmf_var'//cloop//'.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Attribute>'
          write( omominxyz, '(a)' ) '</Grid>'
        end do
        write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

END SUBROUTINE mominxmf_header_Al


SUBROUTINE mominxmf_var_mom( imom, is, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output mom in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_imomisloop

  integer, intent(in) :: imom, is, loop_sta, loop_end, loop_skip

  real(kind=DP), dimension(:,:,:), allocatable :: mom_xyz
  complex(kind=DP), dimension(:,:,:), allocatable :: mom
  character(len=8) :: cloop
  character(len=1) :: cimom
  character(len=1) :: cis
  integer :: loop

    allocate(mom_xyz(0:2*nxw,0:2*nyw,-nzw:nzw))
    allocate(mom(-nx:nx,0:global_ny,-global_nz:global_nz-1))
    write( cimom, '(i1.1)' ) imom
    write( cis, '(i1.1)' ) is

    do loop = loop_sta, loop_end, loop_skip
      write( cloop, '(i8.8)' ) loop
      call rb_mom_imomisloop( imom, is, loop, mom )
      call phi_kxkyz2xyz_fluxtube( mom, mom_xyz )
      open( omominxyz, file='./data/mominxmf_mom'//cimom//'s'//cis//'_var'//cloop//'.bin', status='replace',  &
                       action='write', form='unformatted', access='stream',               &
                       convert='LITTLE_ENDIAN' )
        write( omominxyz ) real(mom_xyz, kind=4)
      close( omominxyz )
    end do

    deallocate(mom_xyz)
    deallocate(mom)

END SUBROUTINE mominxmf_var_mom


SUBROUTINE mominxmf_header_mom( imom, is, loop_sta, loop_end, loop_skip )
!-------------------------------------------------------------------------------
!
!     Output phi in (xx,yy,zz,time) in XMF binary format
!                                                   (S. Maeyama, 30 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_geom, only : lx, ly, lz
  use diag_rb, only : rb_phi_gettime

  integer, intent(in) :: imom, is, loop_sta, loop_end, loop_skip

  real(kind=DP) :: time
  character(len=8) :: cloop
  character(len=3) :: c_alp
  character(len=1) :: cimom
  character(len=1) :: cis
  integer :: loop, i_alp

    write( cimom, '(i1.1)' ) imom
    write( cis, '(i1.1)' ) is

    open( omominxyz, file='./data/mominxmf_mom'//cimom//'s'//cis//'_header_align.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DCORECTMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="ORIGIN_DXDYDZ">'
        write( omominxyz, '(a)' ) '<DataItem Name="Origin" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) 0.d0, 0.d0, 0.d0
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '<DataItem Name="DxDyDz" Format="XML" NumberType="Float" Dimensions="3">'
        write( omominxyz, '(3g17.7e3)' ) max(2*ly,2*lx)/(2*nzw), (2*ly)/(2*nyw), (2*lx)/(2*nxw)
        write( omominxyz, '(a)' ) '</DataItem>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="mom'//cimom//'s'//cis//'">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_mom'//cimom//'s'//cis//'_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/mominxmf_mom'//cimom//'s'//cis//'_header_tube.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        !write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'">'
        write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
        write( omominxyz, '(a)' ) '</Topology>'
        write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_xcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_ycoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_alp000_zcoord.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Geometry>'
        write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="mom'//cimom//'s'//cis//'">'
        write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                         2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
        write( omominxyz, '(a)' ) '  mominxmf_mom'//cimom//'s'//cis//'_var'//cloop//'.bin'
        write( omominxyz, '(a)' ) '</DataStructure>'
        write( omominxyz, '(a)' ) '</Attribute>'
        write( omominxyz, '(a)' ) '</Grid>'
        !write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

    open( omominxyz, file='./data/mominxmf_mom'//cimom//'s'//cis//'_header_full.xmf' )
      call rb_phi_gettime( loop_sta, time )
      write( omominxyz, '(a)' ) '<?xml version="1.0"?>'
      write( omominxyz, '(a)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Xdmf>'
      write( omominxyz, '(a)' ) '<Domain>'
      write( omominxyz, '(a)' ) ''
      write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'" GridType="Collection" CollectionType="Temporal">'
      write( omominxyz, '(a)' ) '<Time TimeType="HyperSlab">'
      write( omominxyz, '(a)' ) '<DataItem Name="Time" Format="XML" NumberType="Float" Dimensions="3">'
      write( omominxyz, '(2g17.7e3,i17)' ) time, dtout_ptn*real(loop_skip, kind=DP), int((loop_end-loop_sta)/loop_skip)+1
      write( omominxyz, '(a)' ) '</DataItem>'
      write( omominxyz, '(a)' ) '</Time>'
      write( omominxyz, '(a)' ) ''
      do loop = loop_sta, loop_end, loop_skip
        write( cloop, '(i8.8)' ) loop
        write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'" GridType="Collection" CollectionType="Spatial">'
        do i_alp = 0, n_alp-1
          write( c_alp, '(i3.3)' ) i_alp
          write( omominxyz, '(a)' ) '<Grid Name="mom'//cimom//'s'//cis//'_var'//cloop//'_alp'//c_alp//'">'
          write( omominxyz, '(a,3i17,a)' ) '<Topology TopologyType="3DSMesh" Dimensions="', 2*nzw+1, 2*nyw+1, 2*nxw+1, '">'
          write( omominxyz, '(a)' ) '</Topology>'
          write( omominxyz, '(a)' ) '<Geometry GeometryType="X_Y_Z">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_xcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_ycoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_alp'//c_alp//'_zcoord.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Geometry>'
          write( omominxyz, '(a)' ) '<Attribute Active="1" Type="Scalar" Center="Node" Name="mom'//cimom//'s'//cis//'">'
          write( omominxyz, '(a,3i17,a)' ) '<DataStructure DataType="Float" Precision="4" Dimensions="',  &
                                           2*nzw+1, 2*nyw+1, 2*nxw+1, '" Format="Binary" Endian="Little">'
          write( omominxyz, '(a)' ) '  mominxmf_mom'//cimom//'s'//cis//'_var'//cloop//'.bin'
          write( omominxyz, '(a)' ) '</DataStructure>'
          write( omominxyz, '(a)' ) '</Attribute>'
          write( omominxyz, '(a)' ) '</Grid>'
        end do
        write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Spatial" -->'
        write( omominxyz, '(a)' ) ''
      end do
      write( omominxyz, '(a)' ) '</Grid><!-- End GridType="Collection" CollectionType="Temporal" -->'
      write( omominxyz, '(a)' ) '</Domain>'
      write( omominxyz, '(a)' ) '</Xdmf>'
    close( omominxyz )

END SUBROUTINE mominxmf_header_mom




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
    write(*,*) '# Larmor radius rho/r_major = ', rho
    if (lx*rho > eps_r) then
      write(*,*) '# WARNING in out_mominvtk. lx*rho < eps_r is recommended. Set larger n_alp.'
      write(*,*) '# lx=',lx,', rho=',rho,', eps_r=',eps_r,', n_alp=',n_alp 
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
    write(*,*) '# Larmor radius rho/r_major = ', rho
    if (lx*rho > eps_r) then
      write(*,*) '# WARNING in out_mominvtk. lx*rho < eps_r is recommended. Set larger n_alp.'
      write(*,*) '# lx=',lx,', rho=',rho,', eps_r=',eps_r,', n_alp=',n_alp 
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


END MODULE out_mominxmf
