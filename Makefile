### Fujitsu Fortran Compiler ###
FC = mpifrtpx
FFLAGS = -Kfast,parallel # Optimization
FFLAGS += -X9 # Fortran95
FFLAGS += -CcdII8 # Change default integer type as integer(kind=8).
FFLAGS += -Koptmsg=2 -Nlst=t # Optimization report
FFLAGS += -fw # Suppress message
FFLAGS += -Kopenmp #-Nfjomplib # OpenMP
FFLAGS += -mcmodel=large # Static memory larger than 2GB
#FFLAGS += -Haefosux -NRtrap #-O0 # Debug
OPTRPT = 'lst'
#FFLAGS += -Nfjprof # Fujitsu profiler fapp
#FFLAGS += -Ksimd_nouse_multiple_structures # Specific option for compiler tcs1.2.26 to avoid slowing down GKV
#FFLAGS += -Knosch_pre_ra # Specific option for compiler tcs1.2.26 to avoid slowing down GKV

### Usage of FFTW (module load fftw-tune)
#INC += -I$(FFTW_DIR)/include
#LIB += -L$(FFTW_DIR)/lib -lfftw3 -lm
LIB += -lfftw3 -lm
### Usage of NetCDF (module load netcdf-fortran netcdf-c phdf5)
### NetCDF does not work on the FLOW supercomputer for now, Jan 17 2021
#INC += -I$(NETCDF_FORTRAN_DIR)/include -I$(NETCDF_DIR)/include -I$(PHDF5_DIR)/include
#LIB += -L$(NETCDF_FORTRAN_DIR)/lib -L$(NETCDF_DIR)/lib -L$(PHDF5_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -Wl,-rpath-link /opt/FJSVxtclanga/tcsds-1.2.27/lib64/
LIB += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5


### Read Fortran binary GKV output
FILEIO = diag_rb_fortran
### Read NetCDF GKV output
#FILEIO = diag_rb_netcdf


SRC = ./src
PROG = diag.exe

diag:	${SRC}/diag_header.f90\
	${SRC}/diag_functions.f90\
	${SRC}/diag_clock.f90\
	${SRC}/gkvp_igs.f90\
	${SRC}/gkvp_vmecbzx.f90\
	${SRC}/diag_geom.f90\
	${SRC}/diag_intgrl.f90\
	${SRC}/diag_parity.f90\
	${SRC}/${FILEIO}.f90\
	${SRC}/diag_fft.f90\
	${SRC}/out_textfile.f90\
	${SRC}/out_netcdf.f90\
	${SRC}/out_mominxy.f90\
	${SRC}/out_mominxz.f90\
	${SRC}/out_mominkxky.f90\
	${SRC}/out_trninkxky.f90\
	${SRC}/out_triinkxky.f90\
	${SRC}/out_mominz.f90\
	${SRC}/out_momintky.f90\
	${SRC}/out_mominavs.f90\
	${SRC}/out_mominxmf.f90\
	${SRC}/out_mominvtk.f90\
	${SRC}/out_mominrz.f90\
	${SRC}/out_ffinkxky.f90\
	${SRC}/out_ffinvm.f90\
	${SRC}/out_ffinzv.f90\
	${SRC}/out_ffinzvm_vtk.f90\
	${SRC}/out_mominfreq.f90\
	${SRC}/out_linfreq.f90\
	${SRC}/out_correlate.f90\
	${SRC}/out_realsp.f90\
	${SRC}/out_bicoherence.f90\
	${SRC}/out_zfshearing.f90\
	${SRC}/out_zfdensity.f90\
	${SRC}/out_fluidtotaltrans.f90\
	${SRC}/out_fluidsubsptrans.f90\
	${SRC}/out_fluiddetailtrans.f90\
	${SRC}/diag_main.f90

	${FC} ${FFLAGS} -c ${SRC}/diag_header.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_functions.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_clock.f90
	${FC} ${FFLAGS} -c ${SRC}/gkvp_igs.f90
	${FC} ${FFLAGS} -c ${SRC}/gkvp_vmecbzx.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_geom.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_intgrl.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_parity.f90
	${FC} ${FFLAGS} -c ${SRC}/${FILEIO}.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/diag_fft.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/out_textfile.f90
	${FC} ${FFLAGS} -c ${SRC}/out_netcdf.f90 ${INC}
	${FC} ${FFLAGS} -c ${SRC}/out_mominxy.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_mominxz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_trninkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_triinkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_momintky.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_mominavs.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominxmf.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominvtk.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominrz.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinkxky.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinvm.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinzv.f90
	${FC} ${FFLAGS} -c ${SRC}/out_ffinzvm_vtk.f90
	${FC} ${FFLAGS} -c ${SRC}/out_mominfreq.f90
	${FC} ${FFLAGS} -c ${SRC}/out_linfreq.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_correlate.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_realsp.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_bicoherence.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_zfshearing.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_zfdensity.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_fluidtotaltrans.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_fluidsubsptrans.f90
#	${FC} ${FFLAGS} -c ${SRC}/out_fluiddetailtrans.f90
	${FC} ${FFLAGS} -c ${SRC}/diag_main.f90

	${FC} ${FFLAGS} *.o -o ${PROG} ${LIB}

	rm -f *.o *.mod *.${OPTRPT}

clean:
	rm -f *.o *.mod *.${OPTRPT} *.exe *.out *.err *.stats go.diag_*.out.*
