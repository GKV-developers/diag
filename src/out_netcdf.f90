MODULE out_netcdf
!-------------------------------------------------------------------------------
!
!     Output NetCDF data
!                                                   (S. Maeyama, 17 June 2020)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!    fileio by using netCDF4 (enables parallel IO through HDF5 or PnetCDF)
!
!                      Create file: nf90_create  ("ncid" specify file.)
!        Define dimensions in file: nf90_def_dim
!          Define variable in file: nf90_def_var ("varid" specify variable.)
!        End of definition of file: nf90_enddef
!        Parallel data access type: nf90_var_par_access
!                   Write variable: nf90_put_var
!                       Close file: nf90_close
!
!-------------------------------------------------------------------------------
  use diag_header
  use netcdf
  implicit none

  private

  public phiinnetcdf,  Alinnetcdf, mominnetcdf, &
         fxvinnetcdf, cntinnetcdf, trninnetcdf, triinnetcdf


 CONTAINS


SUBROUTINE phiinnetcdf
!-------------------------------------------------------------------------------
!
!     Output phi(kx,ky,zz,time) in netcdf
!                                                   (S. Maeyama, 17 June 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_loop, loop_phi_sta, loop_phi_end
  use diag_geom, only : kx, gky, gzz

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: phi
  character(len=3) :: cnum
  integer :: loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:4), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_tt, varid_rephi, varid_imphi
  integer(kind=4) :: start_time(1), count_time(1), start_phi(1:4), count_phi(1:4)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/phi."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/phi."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx", len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky", len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="zz", len=int(2*global_nz,kind=4), dimid=dimids(3))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",  len=NF90_UNLIMITED,          dimid=dimids(4))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",    xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",    xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="zz",    xtype=NF90_DOUBLE, dimids=dimids(3),   varid=varid_zz)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",     xtype=NF90_DOUBLE, dimids=dimids(4),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="rephi", xtype=NF90_DOUBLE, dimids=dimids(1:4), varid=varid_rephi)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="imphi", xtype=NF90_DOUBLE, dimids=dimids(1:4), varid=varid_imphi)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx, values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky, values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz, values=gzz(-global_nz:global_nz-1))
      !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx, values=kx(-nx:nx), start=(/1/), count=(/2*nx+1/))
      !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky, values=gky(0:global_ny), start=(/1/), count=(/global_ny+1/))
      !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz, values=gzz(-global_nz:global_nz-1), start=(/1/), count=(/2*global_nz/))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(1) = 1 !!! count_time(1) = (/ 1 /) 
      count_phi(:) = (/ 2*nx+1,global_ny+1,2*global_nz,1 /)
      do loop = loop_phi_sta(inum), loop_phi_end(inum)
        call rb_phi_gettime( loop, time )
        call rb_phi_loop( loop, phi )
        !start_time(:) = (/ 1+loop /) 
        !start_phi(:) = (/ 1,1,1,1+loop /)
        start_time(1) =  1+loop-loop_phi_sta(inum) !!! start_time(:) = (/ int(1+loop-loop_phi_sta(inum),kind=4) /) 
        start_phi(:) = (/ 1,1,1,1+loop-loop_phi_sta(inum) /) 
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_rephi, values=dble(phi(:,:,:)), start=start_phi, count=count_phi)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imphi, values=aimag(phi(:,:,:)), start=start_phi, count=count_phi)

        !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=(/1+loop/), count=(/1/))
        !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_rephi, values=dble(phi(:,:,:)), &
        !                                       start=(/1,1,1,1+loop/), count=(/2*nx+1,global_ny+1,2*global_nz,1/))
        !ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imphi, values=aimag(phi(:,:,:)), &
        !                                       start=(/1,1,1,1+loop/), count=(/2*nx+1,global_ny+1,2*global_nz,1/))
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")

      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")

    end do

END SUBROUTINE phiinnetcdf


SUBROUTINE Alinnetcdf
!-------------------------------------------------------------------------------
!
!     Output Al(kx,ky,zz,time) in netcdf
!                                                   (T. Ibuki, 31 August 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_loop, loop_Al_sta, loop_Al_end
  use diag_geom, only : kx, gky, gzz

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: Al
  character(len=3) :: cnum
  integer :: loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:4), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_tt, varid_reAl, varid_imAl
  integer(kind=4) :: start_time(1), count_time(1), start_Al(1:4), count_Al(1:4)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum
      != Create file
      ierr_nf90=nf90_create(path="./data/Al."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/Al."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx", len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky", len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="zz", len=int(2*global_nz,kind=4), dimid=dimids(3))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",  len=NF90_UNLIMITED,          dimid=dimids(4))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",    xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",    xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="zz",    xtype=NF90_DOUBLE, dimids=dimids(3),   varid=varid_zz)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",     xtype=NF90_DOUBLE, dimids=dimids(4),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="reAl", xtype=NF90_DOUBLE, dimids=dimids(1:4), varid=varid_reAl)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="imAl", xtype=NF90_DOUBLE, dimids=dimids(1:4), varid=varid_imAl)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx, values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky, values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz, values=gzz(-global_nz:global_nz-1))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 !!! count_time(:) = (/ 1 /) 
      count_Al(:) = (/ 2*nx+1,global_ny+1,2*global_nz,1 /)
      do loop = loop_Al_sta(inum), loop_Al_end(inum)
        call rb_Al_gettime( loop, time )
        call rb_Al_loop( loop, Al )
        start_time(:) = 1+loop-loop_Al_sta(inum) !!! start_time(:) = (/ 1+loop-loop_Al_sta(inum) /) 
        start_Al(:) = (/ 1,1,1,1+loop-loop_Al_sta(inum) /) 
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_reAl, values=dble(Al(:,:,:)), start=start_Al, count=count_Al)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imAl, values=aimag(Al(:,:,:)), start=start_Al, count=count_Al)
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")

      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")

    end do

END SUBROUTINE Alinnetcdf


SUBROUTINE mominnetcdf
!-------------------------------------------------------------------------------
!
!     Output mom(kx,ky,zz,imom,is,time) in netcdf
!                                                    (T. Ibuki, 31 August 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_loop, loop_mom_sta, loop_mom_end
  use diag_geom, only : kx, gky, gzz

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom
  character(len=3) :: cnum
  integer :: imom, is, loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:6), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_imom, varid_is, varid_tt, varid_remom, varid_immom
  integer(kind=4) :: start_time(1), count_time(1), start_mom(1:6), count_mom(1:6)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/mom."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/mom."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx",   len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky",   len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="zz",   len=int(2*global_nz,kind=4), dimid=dimids(3))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="imom", len=int(nmom,kind=4),        dimid=dimids(4)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="is",   len=int(ns,kind=4),          dimid=dimids(5)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",    len=NF90_UNLIMITED,          dimid=dimids(6))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",    xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",    xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="zz",    xtype=NF90_DOUBLE, dimids=dimids(3),   varid=varid_zz)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="imom",  xtype=NF90_INT,    dimids=dimids(4),   varid=varid_imom)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="is",    xtype=NF90_INT,    dimids=dimids(5),   varid=varid_is)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",     xtype=NF90_DOUBLE, dimids=dimids(6),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="remom", xtype=NF90_DOUBLE, dimids=dimids(1:6), varid=varid_remom)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="immom", xtype=NF90_DOUBLE, dimids=dimids(1:6), varid=varid_immom)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx,   values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky,   values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz,   values=gzz(-global_nz:global_nz-1))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imom, values=(/ (imom, imom=0,nmom-1) /))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_is,   values=(/ (is, is=0,ns-1) /))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 
      count_mom(:) = (/ 2*nx+1,global_ny+1,2*global_nz,nmom,ns,1 /)
      do loop = loop_mom_sta(inum), loop_mom_end(inum)
        call rb_mom_gettime( loop, time )
        call rb_mom_loop( loop, mom )
        start_time(:) = 1+loop-loop_mom_sta(inum)
        start_mom(:) = (/ 1,1,1,1,1,1+loop-loop_mom_sta(inum) /) 
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_remom, values=dble(mom(:,:,:,:,:)), start=start_mom, count=count_mom)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_immom, values=aimag(mom(:,:,:,:,:)), start=start_mom, count=count_mom)
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")
      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")
    end do
END SUBROUTINE mominnetcdf


SUBROUTINE fxvinnetcdf
!-------------------------------------------------------------------------------
!
!     Output fxv(kx,ky,zz,vl,mu,is,time) in netcdf
!     Note that zz is thinned out, only iz=0 for each rankz.
!                                                    (S. Maeyama, 15 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_fxv_gettime, rb_fxv_rankzivimisloop, loop_fxv_sta, loop_fxv_end
  use diag_geom, only : kx, gky, gzz, gvl, gmu

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,0:nprocz-1) :: fxv
  character(len=3) :: cnum
  integer :: rankz, iv, im, is, loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:7), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is, varid_tt, varid_refxv, varid_imfxv
  integer(kind=4) :: start_time(1), count_time(1), start_fxv(1:7), count_fxv(1:7)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/fxv."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/fxv."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx", len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky", len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="zz", len=int(nprocz,kind=4),      dimid=dimids(3))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="vl", len=int(2*global_nv,kind=4), dimid=dimids(4)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="mu", len=int(global_nm+1,kind=4), dimid=dimids(5)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="is", len=int(ns,kind=4),          dimid=dimids(6)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",  len=NF90_UNLIMITED,          dimid=dimids(7))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",    xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",    xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="zz",    xtype=NF90_DOUBLE, dimids=dimids(3),   varid=varid_zz)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="vl",    xtype=NF90_DOUBLE, dimids=dimids(4),   varid=varid_vl)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="mu",    xtype=NF90_DOUBLE, dimids=dimids(5),   varid=varid_mu)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="is",    xtype=NF90_INT,    dimids=dimids(6),   varid=varid_is)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",     xtype=NF90_DOUBLE, dimids=dimids(7),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="refxv", xtype=NF90_DOUBLE, dimids=dimids(1:7), varid=varid_refxv)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="imfxv", xtype=NF90_DOUBLE, dimids=dimids(1:7), varid=varid_imfxv)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx, values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky, values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz, values=gzz(-global_nz:global_nz-1:2*nz))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_vl, values=gvl(1:2*global_nv))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_mu, values=gmu(0:global_nm))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_is, values=(/ (is, is=0,ns-1) /))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 
      count_fxv(:) = int((/ 2*nx+1,global_ny+1,nprocz,1,1,1,1 /), kind=4)
      do loop = loop_fxv_sta(inum), loop_fxv_end(inum)
        call rb_fxv_gettime( loop, time )
        start_time(:) = 1+loop-loop_fxv_sta(inum)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        do is = 0, ns-1
          do im = 0, global_nm
            do iv = 1, 2*global_nv
              do rankz = 0, nprocz-1
                call rb_fxv_rankzivimisloop( rankz, iv, im, is, loop, fxv(:,:,rankz) )
              end do
              start_fxv(:) = (/ 1,1,1,iv,1+im,1+is,1+loop-loop_fxv_sta(inum) /) 
              ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_refxv, values=dble(fxv(:,:,:)), start=start_fxv, count=count_fxv)
              ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imfxv, values=aimag(fxv(:,:,:)), start=start_fxv, count=count_fxv)
            end do
          end do
        end do
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")
      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")
    end do
END SUBROUTINE fxvinnetcdf


SUBROUTINE cntinnetcdf
!-------------------------------------------------------------------------------
!
!     Output cnt(kx,ky,zz,vl,mu,is,time) in netcdf
!                                                    (S. Maeyama, 15 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_ivimisloop, loop_cnt_sta, loop_cnt_end
  use diag_geom, only : kx, gky, gzz, gvl, gmu

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: cnt
  character(len=3) :: cnum
  integer :: iv, im, is, loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:7), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_zz, varid_vl, varid_mu, varid_is, varid_tt, varid_recnt, varid_imcnt
  integer(kind=4) :: start_time(1), count_time(1), start_cnt(1:7), count_cnt(1:7)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/cnt."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/cnt."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx", len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky", len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="zz", len=int(2*global_nz,kind=4), dimid=dimids(3))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="vl", len=int(2*global_nv,kind=4), dimid=dimids(4)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="mu", len=int(global_nm+1,kind=4), dimid=dimids(5)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="is", len=int(ns,kind=4),          dimid=dimids(6)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",  len=NF90_UNLIMITED,          dimid=dimids(7))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",    xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",    xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="zz",    xtype=NF90_DOUBLE, dimids=dimids(3),   varid=varid_zz)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="vl",    xtype=NF90_DOUBLE, dimids=dimids(4),   varid=varid_vl)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="mu",    xtype=NF90_DOUBLE, dimids=dimids(5),   varid=varid_mu)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="is",    xtype=NF90_INT,    dimids=dimids(6),   varid=varid_is)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",     xtype=NF90_DOUBLE, dimids=dimids(7),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="recnt", xtype=NF90_DOUBLE, dimids=dimids(1:7), varid=varid_recnt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="imcnt", xtype=NF90_DOUBLE, dimids=dimids(1:7), varid=varid_imcnt)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx, values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky, values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_zz, values=gzz(-global_nz:global_nz-1))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_vl, values=gvl(1:2*global_nv))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_mu, values=gmu(0:global_nm))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_is, values=(/ (is, is=0,ns-1) /))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 
      count_cnt(:) = int((/ 2*nx+1,global_ny+1,2*global_nz,1,1,1,1 /), kind=4)
      do loop = loop_cnt_sta(inum), loop_cnt_end(inum)
        call rb_cnt_gettime( loop, time )
        start_time(:) = 1+loop-loop_cnt_sta(inum)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        do is = 0, ns-1
          do im = 0, global_nm
            do iv = 1, 2*global_nv
              call rb_cnt_ivimisloop( iv, im, is, loop, cnt )
              start_cnt(:) = (/ 1,1,1,iv,1+im,1+is,1+loop-loop_cnt_sta(inum) /) 
              ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_recnt, values=dble(cnt(:,:,:)), start=start_cnt, count=count_cnt)
              ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_imcnt, values=aimag(cnt(:,:,:)), start=start_cnt, count=count_cnt)
            end do
          end do
        end do
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")
      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")
    end do
END SUBROUTINE cntinnetcdf


SUBROUTINE trninnetcdf
!-------------------------------------------------------------------------------
!
!     Output trn(kx,ky,itrn,is,time) in netcdf
!     NOTE that trn is real.
!                                                    (S. Maeyama, 15 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_trn_gettime, rb_trn_loop, loop_trn_sta, loop_trn_end
  use diag_geom, only : kx, gky

  real(kind=DP) :: time
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:ntrn-1,0:ns-1) :: trn
  character(len=3) :: cnum
  integer :: itrn, is, loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:5), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_itrn, varid_is, varid_tt, varid_trn
  integer(kind=4) :: start_time(1), count_time(1), start_trn(1:5), count_trn(1:5)

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/trn."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/trn."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx",   len=int(2*nx+1,kind=4),      dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky",   len=int(global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="itrn", len=int(ntrn,kind=4),        dimid=dimids(3)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="is",   len=int(ns,kind=4),          dimid=dimids(4)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",    len=NF90_UNLIMITED,          dimid=dimids(5))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",   xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",   xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="itrn", xtype=NF90_INT,    dimids=dimids(3),   varid=varid_itrn)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="is",   xtype=NF90_INT,    dimids=dimids(4),   varid=varid_is)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",    xtype=NF90_DOUBLE, dimids=dimids(5),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="trn",  xtype=NF90_DOUBLE, dimids=dimids(1:5), varid=varid_trn)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx,   values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky,   values=gky(0:global_ny))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_itrn, values=(/ (itrn, itrn=0,ntrn-1) /))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_is,   values=(/ (is, is=0,ns-1) /))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 
      count_trn(:) = (/ 2*nx+1,global_ny+1,ntrn,ns,1 /)
      do loop = loop_trn_sta(inum), loop_trn_end(inum)
        call rb_trn_gettime( loop, time )
        call rb_trn_loop( loop, trn )
        start_time(:) = 1+loop-loop_trn_sta(inum)
        start_trn(:) = (/ 1,1,1,1,1+loop-loop_trn_sta(inum) /) 
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_trn, values=trn(:,:,:,:), start=start_trn, count=count_trn)
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")
      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")
    end do
END SUBROUTINE trninnetcdf


SUBROUTINE triinnetcdf( mxt, myt )
!-------------------------------------------------------------------------------
!
!     Output tri(kx,ky,itri,is,time) in netcdf
!     NOTE that tri is real. ky is extended to -global_ny:global_ny.
!                                                    (S. Maeyama, 15 Oct 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_tri_gettime, rb_tri_mxtmytisloop, loop_tri_sta, loop_tri_end
  use diag_geom, only : kx, gky

  integer, intent(in) :: mxt, myt

  real(kind=DP) :: time
  real(kind=DP), dimension(-nx:nx,-global_ny:global_ny,0:ntri-1,0:ns-1) :: tri
  character(len=3) :: cnum
  character(len=4) :: cmx, cmy
  integer :: itri, is, loop, inum

  integer(kind=4) :: ncid_out1, dimids(1:5), ierr_nf90
  integer(kind=4) :: varid_kx, varid_ky, varid_itri, varid_is, varid_tt, varid_tri
  integer(kind=4) :: start_time(1), count_time(1), start_tri(1:5), count_tri(1:5)

    write( cmx, fmt="(i4.4)" ) mxt
    write( cmy, fmt="(i4.4)" ) myt

    do inum = snum, enum
      write( cnum, fmt="(i3.3)" ) inum

      != Create file
      ierr_nf90=nf90_create(path="./data/tri.mx"//cmx//"my"//cmy//"."//cnum//".nc", cmode=NF90_CLOBBER, ncid=ncid_out1)!, &
      !ierr_nf90=nf90_create(path="./data/tri_mx"//cmx//"my"//cmy//"."//cnum//".nc", cmode=IOR(NF90_NETCDF4,NF90_MPIIO), ncid=ncid_out1)!, &
                            !comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
      call check_nf90err(ierr_nf90, "nf90_create")
  
      != Define dimensions in file
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="kx",   len=int(2*nx+1,kind=4),        dimid=dimids(1))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="ky",   len=int(2*global_ny+1,kind=4), dimid=dimids(2))
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="itri", len=int(ntri,kind=4),          dimid=dimids(3)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="is",   len=int(ns,kind=4),            dimid=dimids(4)) 
      ierr_nf90=nf90_def_dim(ncid=ncid_out1, name="t",    len=NF90_UNLIMITED,            dimid=dimids(5))
      call check_nf90err(ierr_nf90, "nf90_def_dim")
    
      != Define variables in file
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="kx",   xtype=NF90_DOUBLE, dimids=dimids(1),   varid=varid_kx)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="ky",   xtype=NF90_DOUBLE, dimids=dimids(2),   varid=varid_ky)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="itri", xtype=NF90_INT,    dimids=dimids(3),   varid=varid_itri)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="is",   xtype=NF90_INT,    dimids=dimids(4),   varid=varid_is)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="t",    xtype=NF90_DOUBLE, dimids=dimids(5),   varid=varid_tt)
      ierr_nf90=nf90_def_var(ncid=ncid_out1, name="tri",  xtype=NF90_DOUBLE, dimids=dimids(1:5), varid=varid_tri)
      call check_nf90err(ierr_nf90, "nf90_def_var")
  
      != End of definition of file
      ierr_nf90=nf90_enddef(ncid=ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_enddef")
  
      != Write variables: static coordinates x and y
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_kx,   values=kx(-nx:nx))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_ky,   values=(/ -gky(global_ny:1:-1),gky(0:global_ny) /))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_itri, values=(/ (itri, itri=0,ntri-1) /))
      ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_is,   values=(/ (is, is=0,ns-1) /))
      call check_nf90err(ierr_nf90, "nf90_putvar")
  
      !%%% Time step loop %%%
      count_time(:) = 1 
      count_tri(:) = (/ 2*nx+1,2*global_ny+1,ntri,ns,1 /)
      do loop = loop_tri_sta(inum), loop_tri_end(inum)
        call rb_tri_gettime( loop, time )
        do is = 0, ns-1
          call rb_tri_mxtmytisloop( mxt, myt, is, loop, tri(:,:,:,is) )
        end do
        start_time(:) = 1+loop-loop_tri_sta(inum)
        start_tri(:) = (/ 1,1,1,1,1+loop-loop_tri_sta(inum) /) 
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tt, values=(/time/), start=start_time, count=count_time)
        ierr_nf90=nf90_put_var(ncid=ncid_out1, varid=varid_tri, values=tri(:,:,:,:), start=start_tri, count=count_tri)
      end do
      call check_nf90err(ierr_nf90, "nf90_putvar")
      != Close file
      ierr_nf90=nf90_close(ncid_out1)
      call check_nf90err(ierr_nf90, "nf90_close")
    end do
END SUBROUTINE triinnetcdf



SUBROUTINE check_nf90err(werr, comment)
!--------------------------------------
!  Check error message of nf90
  integer(kind=4), intent(in) :: werr
  character(len=*), intent(in) :: comment
    
    if(werr /= nf90_noerr) then 
      write(*,*) comment//" "//trim(nf90_strerror(werr))
      stop
    end if

END SUBROUTINE check_nf90err

END MODULE out_netcdf

