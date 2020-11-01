MODULE out_mominz
!-------------------------------------------------------------------------------
!
!     Output moments in (z)
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public phiinz, Alinz, mominz, phiinz_freq,                             &
         phiinz_connect, Alinz_connect, mominz_connect, eneinz_connect,  &
         phiinz_parity, phiinz_parity_freq, &
         Alinz_parity, Alinz_parity_freq

  character(len=9), parameter :: flag_normalize = "none" ! as-is
  !character(len=9), parameter :: flag_normalize = "phi0" ! phi(z=0) normalize
  !character(len=9), parameter :: flag_normalize = "Al0"  ! Al(z=0) normalize

 CONTAINS


SUBROUTINE phiinz( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: phi
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iz

    call rb_phi_gettime( loop, time )
    call rb_phi_mxmyloop( mx, gmy, loop, phi )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
    if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/phiinz_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[phi]","Im[phi]"
      do iz = -global_nz, global_nz-1
        write( omominz, "(99g17.7e3)" ) gzz(iz), real(phi(iz), kind=DP), aimag(phi(iz))
      end do
    close( omominz )

END SUBROUTINE phiinz


SUBROUTINE Alinz( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: Al
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iz

    call rb_Al_gettime( loop, time )
    call rb_Al_mxmyloop( mx, gmy, loop, Al )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
    if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/Alinz_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[Al]","Im[Al]"
      do iz = -global_nz, global_nz-1
        write( omominz, "(99g17.7e3)" ) gzz(iz), real(Al(iz), kind=DP), aimag(Al(iz))
      end do
    close( omominz )

END SUBROUTINE Alinz


SUBROUTINE mominz( mx, gmy, is, loop )
!-------------------------------------------------------------------------------
!
!     Output mom in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_mxmyimomisloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta

  integer, intent(in) :: mx, gmy, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1,0:nmom-1) :: mom
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iz, imom

    call rb_mom_gettime( loop, time )
    do imom = 0, nmom-1
      call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom) )
    end do
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
    if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/mominz_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[mom]","Im[mom]"
      do iz = -global_nz, global_nz-1
        write( omominz, "(99g17.7e3)" ) gzz(iz), (real(mom(iz,imom), kind=DP), &
                                                 aimag(mom(iz,imom)), imom=0,nmom-1)
      end do
    close( omominz )

END SUBROUTINE mominz


SUBROUTINE phiinz_freq( mx, gmy )
!-------------------------------------------------------------------------------
!
!     Output phi in (z) at mx, gmy
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop,  &
                      loop_phi_sta, loop_phi_end
  use diag_geom, only : kx, gky, gzz

  integer, intent(in) :: mx, gmy

  real(kind=DP) :: time, wtime
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: phi, wph
  complex(kind=DP) :: omg
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iz, loop

    loop = loop_phi_sta(snum)

      call rb_phi_gettime( loop, time )
      call rb_phi_mxmyloop( mx, gmy, loop, phi )
      wtime = time
      wph(:) = phi(:)

    do loop = loop_phi_sta(snum)+1, loop_phi_end(enum)
      call rb_phi_gettime( loop, time )

      write( cmx, fmt="(i4.4)" ) mx
      write( cmy, fmt="(i4.4)" ) gmy
      write( cloop, fmt="(i8.8)" ) loop
      open( omominz, file="./data/phiinz_freq_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
        write( omominz, "(99a17)" ) "#              zz","omega","gamma"
  
        call rb_phi_mxmyloop( mx, gmy, loop, phi )
        do iz = -global_nz+1, global_nz-1
          omg = log( phi(iz) / wph(iz) ) / ( ui * ( wtime - time ) )
          write( omominz, "(99g17.7e3)" ) gzz(iz), real(omg, kind=DP), aimag(omg)
        end do
        wtime = time
        wph(:) = phi(:)

      close( omominz )
  
    end do
  
END SUBROUTINE phiinz_freq


SUBROUTINE phiinz_connect( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta, dj, ck

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: phi
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz

    call rb_phi_gettime( loop, time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/phiinz_connect_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[phi]","Im[phi]"

      if ( dj(gmy) == 0 ) then
  
          call rb_phi_mxmyloop( mx, gmy, loop, phi )
          if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
          if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz), real(phi(iz), kind=DP), aimag(phi(iz))
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            call rb_phi_mxmyloop( mxw, gmy, loop, phi )
            if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
            if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                           real(ck(gmy)**iconnect * phi(iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * phi(iz))
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            call rb_phi_mxmyloop( mxw, gmy, loop, phi )
            if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
            if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                           real(conjg( ck(gmy)**iconnect ) * phi(iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * phi(iz))
            end do
          end do
  
      end if

    close( omominz )

END SUBROUTINE phiinz_connect


SUBROUTINE phiinz_parity( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output phi in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_myloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_parity, only : parity_mom2
  use diag_geom, only : kx, gky, gzz, beta, dj, ck

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) :: phi, phi_even, phi_odd
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz

    call rb_phi_gettime( loop, time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )
    call rb_phi_myloop( gmy, loop, phi)
    if ( trim(flag_normalize) == "phi0" ) phi(:,:) = phi(:,:) / phi0
    if ( trim(flag_normalize) == "Al0" )  phi(:,:) = phi(:,:) / (Al0/sqrt(beta))
    call parity_mom2(gmy, phi, phi_even, phi_odd)

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/phiinz_parity_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[phi]","Im[phi]"

      if ( dj(gmy) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz), real(phi_even(mx,iz), kind=DP), aimag(phi_even(mx,iz)), &
                                                     real(phi_odd(mx,iz),  kind=DP), aimag(phi_odd(mx,iz))
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                           real(ck(gmy)**iconnect * phi_even(mxw,iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * phi_even(mxw,iz)),  &
                                           real(ck(gmy)**iconnect * phi_odd(mxw,iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * phi_odd(mxw,iz))
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                           real(conjg( ck(gmy)**iconnect ) * phi_even(mxw,iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * phi_even(mxw,iz)),  &
                                           real(conjg( ck(gmy)**iconnect ) * phi_odd(mxw,iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * phi_odd(mxw,iz))
            end do
          end do
  
      end if

    close( omominz )

END SUBROUTINE phiinz_parity


SUBROUTINE phiinz_parity_freq( mx, gmy )
!-------------------------------------------------------------------------------
!
!     Output phi in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_phi_gettime, rb_phi_myloop, loop_phi_sta, loop_phi_end
  use diag_parity, only : parity_mom2
  use diag_geom, only : kx, gky, gzz, dj

  integer, intent(in) :: mx, gmy

  real(kind=DP) :: time, wtime
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) ::  &
                                     w1, w2, phi_even, phi_odd, wph_even, wph_odd
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz, loop

    loop = loop_phi_sta(snum)

      call rb_phi_gettime(loop, time)
      call rb_phi_myloop(gmy, loop, w1)
      call parity_mom2(gmy, w1, phi_even, phi_odd)
      wtime = time
      wph_even(:,:) = phi_even(:,:)
      wph_odd(:,:) = phi_odd(:,:)

    do loop = loop_phi_sta(snum)+1, loop_phi_end(enum)
      call rb_phi_gettime(loop, time)
      call rb_phi_myloop(gmy, loop, w1)
      call parity_mom2(gmy, w1, phi_even, phi_odd)
      w1(:,:) = log( phi_even(:,:) / wph_even(:,:) ) / ( ui * ( wtime - time ) ) ! omg_even
      w2(:,:) = log( phi_odd(:,:)  / wph_odd(:,:)  ) / ( ui * ( wtime - time ) ) ! omg_odd
      wtime = time
      wph_even(:,:) = phi_even(:,:)
      wph_odd(:,:) = phi_odd(:,:)

      write( cmx, fmt="(i4.4)" ) mx
      write( cmy, fmt="(i4.4)" ) gmy
      write( cloop, fmt="(i8.8)" ) loop
      open( omominz, file="./data/phiinz_parity_freq_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
        write( omominz, "(99a17)" ) "#              zz","wr_even","gamma_even","wr_odd","gamma_odd"
  
        if ( dj(gmy) == 0 ) then
    
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) gzz(iz), real(w1(mx,iz), kind=DP), aimag(w1(mx,iz)), &
                                                       real(w2(mx,iz), kind=DP), aimag(w2(mx,iz))
            end do
    
        else
    
          connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
          if ( connect_min .ne. 0 ) then
            do iconnect = connect_min, 1, -1
              mxw = mx+iconnect*dj(gmy)
              do iz = -global_nz, global_nz-1
                write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                             real(w1(mxw,iz), kind=DP), aimag(w1(mxw,iz)),  &
                                             real(w2(mxw,iz), kind=DP), aimag(w2(mxw,iz))
              end do
            end do
          end if
    
          connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(gmy)
              do iz = -global_nz, global_nz-1
                write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                             real(w1(mxw,iz), kind=DP), aimag(w1(mxw,iz)),  &
                                             real(w2(mxw,iz), kind=DP), aimag(w2(mxw,iz))
              end do
            end do
    
        end if
  
      close( omominz )

    end do

END SUBROUTINE phiinz_parity_freq


!SUBROUTINE phiinz_parity_old_forkx0( mx, gmy, loop )
!!-------------------------------------------------------------------------------
!!
!!     Output phi in (z) at mx, gmy, loop
!!                                                   (S. Maeyama, 11 Nov. 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
!  use diag_geom, only : kx, gky, gzz, beta, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy, loop
!
!  real(kind=DP) :: time
!  complex(kind=DP), dimension(:), allocatable :: phi
!  complex(kind=DP) :: phi0, Al0, phi_even, phi_odd
!  character(len=4) :: cmx, cmy
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, izw
!
!    call rb_phi_gettime( loop, time )
!    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
!    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )
!
!    write( cmx, fmt="(i4.4)" ) mx
!    write( cmy, fmt="(i4.4)" ) gmy
!    write( cloop, fmt="(i8.8)" ) loop
!    open( omominz, file="./data/phiinz_parity_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
!      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!      write( omominz, "(99a17)" ) "#              zz","Re[phi]","Im[phi]"
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( phi(-global_nz:global_nz-1) )
!
!        call rb_phi_mxmyloop( mx, gmy, loop, phi )
!        if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
!        if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))
!        do iz = -global_nz+1, global_nz-1
!          phi_even = 0.5_DP * (phi(iz) + phi(-iz))
!           phi_odd = 0.5_DP * (phi(iz) - phi(-iz))
!          write( omominz, "(99g17.7e3)" ) gzz(iz), real(phi_even, kind=DP), aimag(phi_even), &
!                                                   real(phi_odd, kind=DP), aimag(phi_odd)
!        end do
!
!        deallocate( phi )
!  
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( phi(-zmin:zmax-1) )
!
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) )
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!              phi(izw) = ck(gmy)**iconnect * phi(izw)
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) ) 
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!              phi(izw) = conjg( ck(gmy)**iconnect ) * phi(izw)
!            end do
!          end do
!
!        if ( trim(flag_normalize) == "phi0" ) phi(:) = phi(:) / phi0
!        if ( trim(flag_normalize) == "Al0" )  phi(:) = phi(:) / (Al0/sqrt(beta))
!        do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!          phi_even = 0.5_DP * (phi(iz) + phi(-iz))
!           phi_odd = 0.5_DP * (phi(iz) - phi(-iz))
!          write( omominz, "(99g17.7e3)" ) dz * real(iz), real(phi_even, kind=DP), aimag(phi_even),  &
!                                                           real(phi_odd, kind=DP), aimag(phi_odd)
!        end do
!
!        deallocate( phi )
!  
!      end if
!
!    close( omominz )
!
!END SUBROUTINE phiinz_parity_old_forkx0
!
!
!SUBROUTINE phiinz_parity_old_forkx0_freq( mx, gmy )
!!-------------------------------------------------------------------------------
!!
!!     Output phi in (z) at mx, gmy
!!                                                   (S. Maeyama, 11 Nov. 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_phi_gettime, rb_phi_mxmyloop,  &
!                      loop_phi_sta, loop_phi_end
!  use diag_geom, only : kx, gky, gzz, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy
!
!  real(kind=DP) :: time, wtime
!  complex(kind=DP), dimension(:), allocatable :: phi, wph
!  complex(kind=DP) :: phi_even, phi_odd, wph_even, wph_odd, omg_even, omg_odd
!  character(len=4) :: cmx, cmy
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, izw, loop
!
!    loop = loop_phi_sta(snum)
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( phi(-global_nz:global_nz-1) )
!        allocate( wph(-global_nz:global_nz-1) )
!
!        call rb_phi_gettime( loop, time )
!        call rb_phi_mxmyloop( mx, gmy, loop, phi )
!        wtime = time
!        wph(:) = phi(:)
!
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( phi(-zmin:zmax-1) )
!        allocate( wph(-zmin:zmax-1) )
!
!        call rb_phi_gettime( loop, time )
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) )
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!              phi(izw) = ck(gmy)**iconnect * phi(izw)
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) ) 
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!              phi(izw) = conjg( ck(gmy)**iconnect ) * phi(izw)
!            end do
!          end do
!        wtime = time
!        wph(:) = phi(:)
!
!      end if
!
!    do loop = loop_phi_sta(snum)+1, loop_phi_end(enum)
!      call rb_phi_gettime( loop, time )
!
!      write( cmx, fmt="(i4.4)" ) mx
!      write( cmy, fmt="(i4.4)" ) gmy
!      write( cloop, fmt="(i8.8)" ) loop
!      open( omominz, file="./data/phiinz_parity_freq_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
!        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!        write( omominz, "(99a17)" ) "#              zz","Re[phi]","Im[phi]"
!  
!        if ( dj(gmy) == 0 ) then
!    
!          call rb_phi_mxmyloop( mx, gmy, loop, phi )
!          do iz = -global_nz+1, global_nz-1
!            phi_even = 0.5_DP * (phi(iz) + phi(-iz))
!             phi_odd = 0.5_DP * (phi(iz) - phi(-iz))
!            wph_even = 0.5_DP * (wph(iz) + wph(-iz))
!             wph_odd = 0.5_DP * (wph(iz) - wph(-iz))
!            omg_even = log( phi_even / wph_even ) / ( ui * ( wtime - time ) )
!             omg_odd = log( phi_odd  / wph_odd  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) gzz(iz), real(omg_even, kind=DP), aimag(omg_even), &
!                                                     real(omg_odd, kind=DP), aimag(omg_odd)
!          end do
!          wtime = time
!          wph(:) = phi(:)
!  
!        else
!  
!          if ( connect_min .ne. 0 ) then
!            do iconnect = connect_min, 1, -1
!              mxw = mx+iconnect*dj(gmy)
!              izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!              call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) )
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!                phi(izw) = ck(gmy)**iconnect * phi(izw)
!              end do
!            end do
!          end if
!            do iconnect = 0, connect_max
!              mxw = mx-iconnect*dj(gmy)
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!              call rb_phi_mxmyloop( mxw, gmy, loop, phi(izw:izw+2*global_nz) ) 
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!                phi(izw) = conjg( ck(gmy)**iconnect ) * phi(izw)
!              end do
!            end do
!  
!          do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!            phi_even = 0.5_DP * (phi(iz) + phi(-iz))
!             phi_odd = 0.5_DP * (phi(iz) - phi(-iz))
!            wph_even = 0.5_DP * (wph(iz) + wph(-iz))
!             wph_odd = 0.5_DP * (wph(iz) - wph(-iz))
!            omg_even = log( phi_even / wph_even ) / ( ui * ( wtime - time ) )
!             omg_odd = log( phi_odd  / wph_odd  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) - dz * real(iz), real(omg_even, kind=DP), aimag(omg_even),  &
!                                                             real(omg_odd, kind=DP), aimag(omg_odd)
!          end do
!          wtime = time
!          wph(:) = phi(:)
!  
!        end if
!
!      close( omominz )
!
!    end do
!  
!    deallocate( phi )
!    deallocate( wph )
!
!
!END SUBROUTINE phiinz_parity_old_forkx0_freq


SUBROUTINE Alinz_connect( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 11 Nov. 2014)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta, dj, ck

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: Al
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz

    call rb_Al_gettime( loop, time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/Alinz_connect_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[Al]","Im[Al]"

      if ( dj(gmy) == 0 ) then
  
          call rb_Al_mxmyloop( mx, gmy, loop, Al )
          if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
          if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz), real(Al(iz), kind=DP), aimag(Al(iz))
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            call rb_Al_mxmyloop( mxw, gmy, loop, Al )
            if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
            if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                           real(ck(gmy)**iconnect * Al(iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * Al(iz))
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            call rb_Al_mxmyloop( mxw, gmy, loop, Al )
            if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
            if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                           real(conjg( ck(gmy)**iconnect ) * Al(iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * Al(iz))
            end do
          end do
  
      end if

    close( omominz )

END SUBROUTINE Alinz_connect


SUBROUTINE Alinz_parity( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output Al in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_myloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_parity, only : parity_mom2
  use diag_geom, only : kx, gky, gzz, beta, dj, ck

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) :: Al, Al_even, Al_odd
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz

    call rb_Al_gettime( loop, time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )
    call rb_Al_myloop( gmy, loop, Al)
    if ( trim(flag_normalize) == "phi0" ) Al(:,:) = (Al(:,:)/sqrt(beta)) / phi0
    if ( trim(flag_normalize) == "Al0" )  Al(:,:) = Al(:,:) / Al0
    call parity_mom2(gmy, Al, Al_even, Al_odd)

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/Alinz_parity_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[Al]","Im[Al]"

      if ( dj(gmy) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz), real(Al_even(mx,iz), kind=DP), aimag(Al_even(mx,iz)), &
                                                     real(Al_odd(mx,iz),  kind=DP), aimag(Al_odd(mx,iz))
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                           real(ck(gmy)**iconnect * Al_even(mxw,iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * Al_even(mxw,iz)),  &
                                           real(ck(gmy)**iconnect * Al_odd(mxw,iz), kind=DP), &
                                           aimag(ck(gmy)**iconnect * Al_odd(mxw,iz))
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                           real(conjg( ck(gmy)**iconnect ) * Al_even(mxw,iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * Al_even(mxw,iz)),  &
                                           real(conjg( ck(gmy)**iconnect ) * Al_odd(mxw,iz), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * Al_odd(mxw,iz))
            end do
          end do
  
      end if

    close( omominz )

END SUBROUTINE Alinz_parity


SUBROUTINE Alinz_parity_freq( mx, gmy )
!-------------------------------------------------------------------------------
!
!     Output Al in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 18 June 2018)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_Al_gettime, rb_Al_myloop, loop_Al_sta, loop_Al_end
  use diag_parity, only : parity_mom2
  use diag_geom, only : kx, gky, gzz, dj

  integer, intent(in) :: mx, gmy

  real(kind=DP) :: time, wtime
  complex(kind=DP), dimension(-nx:nx,-global_nz:global_nz-1) ::  &
                                     w1, w2, Al_even, Al_odd, wal_even, wal_odd
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz, loop

    loop = loop_Al_sta(snum)

      call rb_Al_gettime(loop, time)
      call rb_Al_myloop(gmy, loop, w1)
      call parity_mom2(gmy, w1, Al_even, Al_odd)
      wtime = time
      wal_even(:,:) = Al_even(:,:)
      wal_odd(:,:) = Al_odd(:,:)

    do loop = loop_Al_sta(snum)+1, loop_Al_end(enum)
      call rb_Al_gettime(loop, time)
      call rb_Al_myloop(gmy, loop, w1)
      call parity_mom2(gmy, w1, Al_even, Al_odd)
      w1(:,:) = log( Al_even(:,:) / wal_even(:,:) ) / ( ui * ( wtime - time ) ) ! omg_even
      w2(:,:) = log( Al_odd(:,:)  / wal_odd(:,:)  ) / ( ui * ( wtime - time ) ) ! omg_odd
      wtime = time
      wal_even(:,:) = Al_even(:,:)
      wal_odd(:,:) = Al_odd(:,:)

      write( cmx, fmt="(i4.4)" ) mx
      write( cmy, fmt="(i4.4)" ) gmy
      write( cloop, fmt="(i8.8)" ) loop
      open( omominz, file="./data/Alinz_parity_freq_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
        write( omominz, "(99a17)" ) "#              zz","wr_even","gamma_even","wr_odd","gamma_odd"
  
        if ( dj(gmy) == 0 ) then
    
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) gzz(iz), real(w1(mx,iz), kind=DP), aimag(w1(mx,iz)), &
                                                       real(w2(mx,iz), kind=DP), aimag(w2(mx,iz))
            end do
    
        else
    
          connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
          if ( connect_min .ne. 0 ) then
            do iconnect = connect_min, 1, -1
              mxw = mx+iconnect*dj(gmy)
              do iz = -global_nz, global_nz-1
                write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                             real(w1(mxw,iz), kind=DP), aimag(w1(mxw,iz)),  &
                                             real(w2(mxw,iz), kind=DP), aimag(w2(mxw,iz))
              end do
            end do
          end if
    
          connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(gmy)
              do iz = -global_nz, global_nz-1
                write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                             real(w1(mxw,iz), kind=DP), aimag(w1(mxw,iz)),  &
                                             real(w2(mxw,iz), kind=DP), aimag(w2(mxw,iz))
              end do
            end do
    
        end if
  
      close( omominz )

    end do

END SUBROUTINE Alinz_parity_freq


!SUBROUTINE Alinz_parity_old_forkx0( mx, gmy, loop )
!!-------------------------------------------------------------------------------
!!
!!     Output Al in (z) at mx, gmy, loop
!!                                                   (S. Maeyama, 11 Nov. 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_Al_gettime, rb_Al_mxmyloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
!  use diag_geom, only : kx, gky, gzz, beta, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy, loop
!
!  real(kind=DP) :: time
!  complex(kind=DP), dimension(:), allocatable :: Al
!  complex(kind=DP) :: phi0, Al0, Al_even, Al_odd
!  character(len=4) :: cmx, cmy
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, izw
!
!    call rb_Al_gettime( loop, time )
!    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
!    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )
!
!    write( cmx, fmt="(i4.4)" ) mx
!    write( cmy, fmt="(i4.4)" ) gmy
!    write( cloop, fmt="(i8.8)" ) loop
!    open( omominz, file="./data/Alinz_parity_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
!      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!      write( omominz, "(99a17)" ) "#              zz","Re[Al]","Im[Al]"
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( Al(-global_nz:global_nz-1) )
!
!        call rb_Al_mxmyloop( mx, gmy, loop, Al )
!        if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
!        if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0
!        do iz = -global_nz+1, global_nz-1
!          Al_even = 0.5_DP * (Al(iz) + Al(-iz))
!           Al_odd = 0.5_DP * (Al(iz) - Al(-iz))
!          write( omominz, "(99g17.7e3)" ) gzz(iz), real(Al_even, kind=DP), aimag(Al_even), &
!                                                   real(Al_odd, kind=DP), aimag(Al_odd)
!        end do
!
!        deallocate( Al )
!  
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( Al(-zmin:zmax-1) )
!
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) )
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!              Al(izw) = ck(gmy)**iconnect * Al(izw)
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) ) 
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!              Al(izw) = conjg( ck(gmy)**iconnect ) * Al(izw)
!            end do
!          end do
!
!        if ( trim(flag_normalize) == "phi0" ) Al(:) = (Al(:)/sqrt(beta)) / phi0
!        if ( trim(flag_normalize) == "Al0" )  Al(:) = Al(:) / Al0
!        do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!          Al_even = 0.5_DP * (Al(iz) + Al(-iz))
!           Al_odd = 0.5_DP * (Al(iz) - Al(-iz))
!          write( omominz, "(99g17.7e3)" ) dz * real(iz), real(Al_even, kind=DP), aimag(Al_even),  &
!                                                         real(Al_odd, kind=DP), aimag(Al_odd)
!        end do
!
!        deallocate( Al )
!  
!      end if
!
!    close( omominz )
!
!END SUBROUTINE Alinz_parity_old_forkx0
!
!
!SUBROUTINE Alinz_parity_old_forkx0_freq( mx, gmy )
!!-------------------------------------------------------------------------------
!!
!!     Output Al in (z) at mx, gmy
!!                                                   (S. Maeyama, 11 Nov. 2014)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_Al_gettime, rb_Al_mxmyloop,  &
!                      loop_Al_sta, loop_Al_end
!  use diag_geom, only : kx, gky, gzz, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy
!
!  real(kind=DP) :: time, wtime
!  complex(kind=DP), dimension(:), allocatable :: Al, wal
!  complex(kind=DP) :: Al_even, Al_odd, wal_even, wal_odd, omg_even, omg_odd
!  character(len=4) :: cmx, cmy
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, izw, loop
!
!    loop = loop_Al_sta(snum)
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( Al(-global_nz:global_nz-1) )
!        allocate( wal(-global_nz:global_nz-1) )
!
!        call rb_Al_gettime( loop, time )
!        call rb_Al_mxmyloop( mx, gmy, loop, Al )
!        wtime = time
!        wal(:) = Al(:)
!
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( Al(-zmin:zmax-1) )
!        allocate( wal(-zmin:zmax-1) )
!
!        call rb_Al_gettime( loop, time )
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) )
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!              Al(izw) = ck(gmy)**iconnect * Al(izw)
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) ) 
!            do iz = -global_nz, global_nz-1
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!              Al(izw) = conjg( ck(gmy)**iconnect ) * Al(izw)
!            end do
!          end do
!        wtime = time
!        wal(:) = Al(:)
!
!      end if
!
!    do loop = loop_Al_sta(snum)+1, loop_Al_end(enum)
!      call rb_Al_gettime( loop, time )
!
!      write( cmx, fmt="(i4.4)" ) mx
!      write( cmy, fmt="(i4.4)" ) gmy
!      write( cloop, fmt="(i8.8)" ) loop
!      open( omominz, file="./data/Alinz_parity_freq_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
!        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!        write( omominz, "(99a17)" ) "#              zz","Re[Al]","Im[Al]"
!  
!        if ( dj(gmy) == 0 ) then
!    
!          call rb_Al_mxmyloop( mx, gmy, loop, Al )
!          do iz = -global_nz+1, global_nz-1
!            Al_even = 0.5_DP * (Al(iz) + Al(-iz))
!             Al_odd = 0.5_DP * (Al(iz) - Al(-iz))
!            wal_even = 0.5_DP * (wal(iz) + wal(-iz))
!             wal_odd = 0.5_DP * (wal(iz) - wal(-iz))
!            omg_even = log( Al_even / wal_even ) / ( ui * ( wtime - time ) )
!             omg_odd = log( Al_odd  / wal_odd  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) gzz(iz), real(omg_even, kind=DP), aimag(omg_even), &
!                                                     real(omg_odd, kind=DP), aimag(omg_odd)
!          end do
!          wtime = time
!          wal(:) = Al(:)
!  
!        else
!  
!          if ( connect_min .ne. 0 ) then
!            do iconnect = connect_min, 1, -1
!              mxw = mx+iconnect*dj(gmy)
!              izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!              call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) )
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!                Al(izw) = ck(gmy)**iconnect * Al(izw)
!              end do
!            end do
!          end if
!            do iconnect = 0, connect_max
!              mxw = mx-iconnect*dj(gmy)
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!              call rb_Al_mxmyloop( mxw, gmy, loop, Al(izw:izw+2*global_nz) ) 
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!                Al(izw) = conjg( ck(gmy)**iconnect ) * Al(izw)
!              end do
!            end do
!  
!          do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!            Al_even = 0.5_DP * (Al(iz) + Al(-iz))
!             Al_odd = 0.5_DP * (Al(iz) - Al(-iz))
!            wal_even = 0.5_DP * (wal(iz) + wal(-iz))
!             wal_odd = 0.5_DP * (wal(iz) - wal(-iz))
!            omg_even = log( Al_even / wal_even ) / ( ui * ( wtime - time ) )
!             omg_odd = log( Al_odd  / wal_odd  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) - dz * real(iz), real(omg_even, kind=DP), aimag(omg_even),  &
!                                                             real(omg_odd, kind=DP), aimag(omg_odd)
!          end do
!          wtime = time
!          wal(:) = Al(:)
!  
!        end if
!
!      close( omominz )
!
!    end do
!  
!    deallocate( Al )
!    deallocate( wal )
!
!
!END SUBROUTINE Alinz_parity_old_forkx0_freq


SUBROUTINE mominz_connect( mx, gmy, is, loop )
!-------------------------------------------------------------------------------
!
!     Output mom in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 9 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_mom_mxmyimomisloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
  use diag_geom, only : kx, gky, gzz, beta, dj, ck

  integer, intent(in) :: mx, gmy, is, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1,0:nmom-1) :: mom
  complex(kind=DP) :: phi0, Al0
  character(len=4) :: cmx, cmy
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz, imom

    call rb_mom_gettime( loop, time )
    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cis, fmt="(i1.1)" ) is
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/mominz_connect_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[mom]","Im[mom]"

      if ( dj(gmy) == 0 ) then
  
          do imom = 0, nmom-1
            call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom) )
          end do
          if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
          if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz), (real(mom(iz,imom), kind=DP), &
                                                     aimag(mom(iz,imom)), imom=0,nmom-1)
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            do imom = 0, nmom-1
              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(:,imom) )
            end do
            if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
            if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                           (real(ck(gmy)**iconnect * mom(iz,imom), kind=DP), &
                                           aimag(ck(gmy)**iconnect * mom(iz,imom)), imom=0,nmom-1)
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            do imom = 0, nmom-1
              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(:,imom) )
            end do
            if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
            if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                           (real(conjg( ck(gmy)**iconnect ) * mom(iz,imom), kind=DP), &
                                           aimag(conjg( ck(gmy)**iconnect ) * mom(iz,imom)), imom=0,nmom-1)
            end do
          end do
  
      end if

    close( omominz )

END SUBROUTINE mominz_connect


!SUBROUTINE mominz_parity_old_forkx0( mx, gmy, is, loop )
!!-------------------------------------------------------------------------------
!!
!!     Output mom in (z) at mx, gmy, imom, is, loop
!!                                                   (S. Maeyama, 9 Oct. 2015)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_mom_gettime, rb_mom_mxmyimomisloop, rb_phi_mxmyizloop, rb_Al_mxmyizloop
!  use diag_geom, only : kx, gky, gzz, beta, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy, is, loop
!
!  real(kind=DP) :: time
!  complex(kind=DP), dimension(:,:), allocatable :: mom
!  complex(kind=DP) :: phi0, Al0, mom_even(0:nmom-1), mom_odd(0:nmom-1)
!  character(len=4) :: cmx, cmy
!  character(len=1) :: cis
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, izw, imom
!
!    call rb_mom_gettime( loop, time )
!    call rb_phi_mxmyizloop( mx, gmy, 0, loop, phi0 )
!    call rb_Al_mxmyizloop( mx, gmy, 0, loop, Al0 )
!
!    write( cmx, fmt="(i4.4)" ) mx
!    write( cmy, fmt="(i4.4)" ) gmy
!    write( cis, fmt="(i1.1)" ) is
!    write( cloop, fmt="(i8.8)" ) loop
!    open( omominz, file="./data/mominz_parity_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".dat" )
!      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!      write( omominz, "(99a17)" ) "#              zz","Re[mom]","Im[mom]"
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( mom(-global_nz:global_nz-1,0:nmom-1) )
!
!        do imom = 0, nmom-1
!          call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom) )
!        end do
!        if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
!        if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))
!        do iz = -global_nz+1, global_nz-1
!          mom_even(:) = 0.5_DP * (mom(iz,:) + mom(-iz,:))
!           mom_odd(:) = 0.5_DP * (mom(iz,:) - mom(-iz,:))
!          write( omominz, "(99g17.7e3)" ) gzz(iz), (real(mom_even(imom), kind=DP), &
!                                                   aimag(mom_even(imom)), imom=0,nmom-1), &
!                                                   (real(mom_odd(imom), kind=DP), &
!                                                   aimag(mom_odd(imom)), imom=0,nmom-1)
!        end do
!
!        deallocate( mom )
!  
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( mom(-zmin:zmax-1,0:nmom-1) )
!
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            do imom = 0, nmom-1
!              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) )
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!                mom(izw,imom) = ck(gmy)**iconnect * mom(izw,imom)
!              end do
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            do imom = 0, nmom-1
!              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) ) 
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!                mom(izw,imom) = conjg( ck(gmy)**iconnect ) * mom(izw,imom)
!              end do
!            end do
!          end do
!
!        if ( trim(flag_normalize) == "phi0" ) mom(:,:) = mom(:,:) / phi0
!        if ( trim(flag_normalize) == "Al0" )  mom(:,:) = mom(:,:) / (Al0/sqrt(beta))
!        do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!          mom_even(:) = 0.5_DP * (mom(iz,:) + mom(-iz,:))
!           mom_odd(:) = 0.5_DP * (mom(iz,:) - mom(-iz,:))
!          write( omominz, "(99g17.7e3)" ) - dz * real(iz), ((real(mom_even(imom), kind=DP), &
!                                                           aimag(mom_even(imom))), imom=0,nmom-1), &
!                                                           ((real(mom_odd(imom), kind=DP), &
!                                                           aimag(mom_odd(imom))), imom=0,nmom-1)
!        end do
!
!        deallocate( mom )
!  
!      end if
!
!    close( omominz )
!
!END SUBROUTINE mominz_parity_old_forkx0
!
!
!SUBROUTINE mominz_parity_old_forkx0_freq( mx, gmy, is )
!!-------------------------------------------------------------------------------
!!
!!     Output mom in (z) at mx, gmy, imom, is
!!                                                   (S. Maeyama, 9 Oct. 2015)
!!
!!-------------------------------------------------------------------------------
!  use diag_rb, only : rb_mom_gettime, rb_mom_mxmyimomisloop,  &
!                      loop_mom_sta, loop_mom_end
!  use diag_geom, only : kx, gky, gzz, dj, ck, dz
!
!  integer, intent(in) :: mx, gmy, is
!
!  real(kind=DP) :: time, wtime
!  complex(kind=DP), dimension(:,:), allocatable :: mom, wmm
!  complex(kind=DP), dimension(0:nmom-1) :: mom_even, mom_odd, wmm_even, wmm_odd, omg_even, omg_odd
!  character(len=4) :: cmx, cmy
!  character(len=1) :: cis
!  character(len=8) :: cloop
!  integer :: iconnect, connect_min, connect_max, mxw, zmin, zmax
!  integer :: iz, imom, izw, loop
!
!    loop = loop_mom_sta(snum)
!
!      if ( dj(gmy) == 0 ) then
!  
!        connect_min = 0
!        connect_max = 0
!        zmin = global_nz
!        zmax = global_nz
!        allocate( mom(-global_nz:global_nz-1,0:nmom-1) )
!        allocate( wmm(-global_nz:global_nz-1,0:nmom-1) )
!
!        call rb_mom_gettime( loop, time )
!        do imom = 0, nmom-1
!          call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom) )
!        end do
!        wtime = time
!        wmm(:,:) = mom(:,:)
!
!      else
!
!        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
!        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
!        zmin = 2 * global_nz * connect_min + global_nz
!        zmax = 2 * global_nz * connect_max + global_nz
!        allocate( mom(-zmin:zmax-1,0:nmom-1) )
!        allocate( wmm(-zmin:zmax-1,0:nmom-1) )
!
!        call rb_mom_gettime( loop, time )
!        if ( connect_min .ne. 0 ) then
!          do iconnect = connect_min, 1, -1
!            mxw = mx+iconnect*dj(gmy)
!            izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!            do imom = 0, nmom-1
!              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) )
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!                mom(izw,imom) = ck(gmy)**iconnect * mom(izw,imom)
!              end do
!            end do
!          end do
!        end if
!          do iconnect = 0, connect_max
!            mxw = mx-iconnect*dj(gmy)
!            izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!            do imom = 0, nmom-1
!              call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) ) 
!              do iz = -global_nz, global_nz-1
!                izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!                mom(izw,imom) = conjg( ck(gmy)**iconnect ) * mom(izw,imom)
!              end do
!            end do
!          end do
!        wtime = time
!        wmm(:,:) = mom(:,:)
!
!      end if
!
!    do loop = loop_mom_sta(snum)+1, loop_mom_end(enum)
!      call rb_mom_gettime( loop, time )
!
!      write( cmx, fmt="(i4.4)" ) mx
!      write( cmy, fmt="(i4.4)" ) gmy
!      write( cis, fmt="(i1.1)" ) is
!      write( cloop, fmt="(i8.8)" ) loop
!      open( omominz, file="./data/mominz_parity_freq_mx"//cmx//"my"//cmy//"s"//cis//"_t"//cloop//".dat" )
!        write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
!        write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
!        write( omominz, "(99a17)" ) "#              zz","Re[mom]","Im[mom]"
!  
!        if ( dj(gmy) == 0 ) then
!    
!          do imom = 0, nmom-1
!            call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom) )
!          end do
!          do iz = -global_nz+1, global_nz-1
!            mom_even(:) = 0.5_DP * (mom(iz,:) + mom(-iz,:))
!             mom_odd(:) = 0.5_DP * (mom(iz,:) - mom(-iz,:))
!            wmm_even(:) = 0.5_DP * (wmm(iz,:) + wmm(-iz,:))
!             wmm_odd(:) = 0.5_DP * (wmm(iz,:) - wmm(-iz,:))
!            omg_even(:) = log( mom_even(:) / wmm_even(:) ) / ( ui * ( wtime - time ) )
!             omg_odd(:) = log( mom_odd(:)  / wmm_odd(:)  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) gzz(iz), (real(omg_even(imom), kind=DP), &
!                                                     aimag(omg_even(imom)), imom=0,nmom-1), &
!                                                     (real(omg_odd(imom), kind=DP), &
!                                                     aimag(omg_odd(imom)), imom=0,nmom-1)
!          end do
!          wtime = time
!          wmm(:,:) = mom(:,:)
!  
!        else
!  
!          if ( connect_min .ne. 0 ) then
!            do iconnect = connect_min, 1, -1
!              mxw = mx+iconnect*dj(gmy)
!              izw = -zmin + 2 * global_nz * (connect_min - iconnect)
!              do imom = 0, nmom-1
!                call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) )
!                do iz = -global_nz, global_nz-1
!                  izw = - zmin + 2 * global_nz * (connect_min - iconnect) + global_nz + iz
!                  mom(izw,imom) = ck(gmy)**iconnect * mom(izw,imom)
!                end do
!              end do
!            end do
!          end if
!            do iconnect = 0, connect_max
!              mxw = mx-iconnect*dj(gmy)
!              izw = - zmin + 2 * global_nz * (connect_min + iconnect)
!              do imom = 0, nmom-1
!                call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(izw:izw+2*global_nz,imom) )
!                do iz = -global_nz, global_nz-1
!                  izw = - zmin + 2 * global_nz * (connect_min + iconnect) + global_nz + iz
!                  mom(izw,imom) = conjg( ck(gmy)**iconnect ) * mom(izw,imom)
!                end do
!              end do
!            end do
!  
!          do iz = -min(zmin,zmax-1), min(zmin,zmax-1)
!            mom_even(:) = 0.5_DP * (mom(iz,:) + mom(-iz,:))
!             mom_odd(:) = 0.5_DP * (mom(iz,:) - mom(-iz,:))
!            wmm_even(:) = 0.5_DP * (wmm(iz,:) + wmm(-iz,:))
!             wmm_odd(:) = 0.5_DP * (wmm(iz,:) - wmm(-iz,:))
!            omg_even(:) = log( mom_even(:) / wmm_even(:) ) / ( ui * ( wtime - time ) )
!             omg_odd(:) = log( mom_odd(:)  / wmm_odd(:)  ) / ( ui * ( wtime - time ) )
!            write( omominz, "(99g17.7e3)" ) - dz * real(iz), (real(omg_even(imom), kind=DP), &
!                                                             aimag(omg_even(imom)), imom=0,nmom-1), &
!                                                             (real(omg_odd(imom), kind=DP), &
!                                                             aimag(omg_odd(imom)), imom=0,nmom-1)
!          end do
!          wtime = time
!          wmm(:,:) = mom(:,:)
!  
!        end if
!
!      close( omominz )
!
!    end do
!  
!    deallocate( mom )
!    deallocate( wmm )
!
!
!END SUBROUTINE mominz_parity_old_forkx0_freq


SUBROUTINE eneinz_connect( mx, gmy, loop )
!-------------------------------------------------------------------------------
!
!     Output mom in (z) at mx, gmy, loop
!                                                   (S. Maeyama, 9 Oct. 2015)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_mom_gettime, rb_phi_mxmyloop, &
                      rb_Al_mxmyloop, rb_mom_mxmyimomisloop
  use diag_geom, only : kx, gky, gzz, dj, Anum, Znum, fcs, tau, &
                        fct_e_energy, fct_m_energy

  integer, intent(in) :: mx, gmy, loop

  real(kind=DP) :: time
  complex(kind=DP), dimension(-global_nz:global_nz-1) :: phi, Al
  complex(kind=DP), dimension(-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: mom
  real(kind=DP), dimension(-global_nz:global_nz-1,0:nmom-1,0:ns-1) :: ene
  character(len=4) :: cmx, cmy
  character(len=8) :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: iz, imom, is
                                         !%%% for debug %%%
                                         !real(kind=DP) :: fct, wek, wmk, ent(0:ns-1)
                                         !wek = 0._DP
                                         !wmk = 0._DP
                                         !ent(:) = 0._DP

    call rb_mom_gettime( loop, time )

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( cloop, fmt="(i8.8)" ) loop
    open( omominz, file="./data/eneinz_connect_mx"//cmx//"my"//cmy//"_t"//cloop//".dat" )
      write( omominz, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( omominz, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( omominz, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( omominz, "(99a17)" ) "#              zz","Re[mom]","Im[mom]"

      if ( dj(gmy) == 0 ) then
  
          call rb_phi_mxmyloop( mx, gmy, loop, phi )
          call rb_Al_mxmyloop( mx, gmy, loop, Al )
          do is = 0, ns-1
            do imom = 0, nmom-1
              call rb_mom_mxmyimomisloop( mx, gmy, imom, is, loop, mom(:,imom,is) )
            end do
            do iz = -global_nz, global_nz-1
              ene(iz,0,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * abs(mom(iz,0,is)*Znum(is)/fcs(is))**2
              ene(iz,1,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * abs(mom(iz,1,is)*Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
              ene(iz,2,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * 0.5_DP*abs((2._DP*mom(iz,2,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
              ene(iz,3,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * abs((mom(iz,3,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
              ene(iz,4,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * (1._DP/6._DP)*abs((2._DP*mom(iz,4,is)/tau(is)-3._DP*mom(iz,1,is))  &
                                                      *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
              ene(iz,5,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                           * abs((mom(iz,5,is)/tau(is)-mom(iz,1,is))  &
                                                      *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
            end do
          end do
          do iz = -global_nz, global_nz-1
            write( omominz, "(99g17.7e3)" ) gzz(iz),                             &
               real(0.5_DP*fct_e_energy(mx,gmy,iz) * abs(phi(iz))**2, kind=DP),  &
               real(0.5_DP*fct_m_energy(mx,gmy,iz) * abs( Al(iz))**2, kind=DP),  &
              ((ene(iz,imom,is), imom=0,nmom-1), is=0,ns-1)
          end do

                                         !do iz = -global_nz, global_nz-1
                                         !  fct = grootg(iz) / cfsrf
                                         !  wek = wek + fct * fct_e_energy(mx,gmy,iz) * abs(phi(iz))**2
                                         !  wmk = wmk + fct * fct_m_energy(mx,gmy,iz) * abs( Al(iz))**2
                                         !end do
                                         !do is = 0, ns-1
                                         !  do imom = 0, nmom-1
                                         !    do iz = -global_nz, global_nz-1
                                         !      fct = grootg(iz) / cfsrf
                                         !      ent(is) = ent(is) + fct * 2._DP * ene(iz,imom,is)
                                         !    end do
                                         !  end do
                                         !end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(gmy) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(gmy)
            call rb_phi_mxmyloop( mxw, gmy, loop, phi )
            call rb_Al_mxmyloop( mxw, gmy, loop, Al )
            do is = 0, ns-1
              do imom = 0, nmom-1
                call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(:,imom,is) )
              end do
              do iz = -global_nz, global_nz-1
                ene(iz,0,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs(mom(iz,0,is)*Znum(is)/fcs(is))**2
                ene(iz,1,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs(mom(iz,1,is)*Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
                ene(iz,2,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * 0.5_DP*abs((2._DP*mom(iz,2,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
                ene(iz,3,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs((mom(iz,3,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
                ene(iz,4,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * (1._DP/6._DP)*abs((2._DP*mom(iz,4,is)/tau(is)-3._DP*mom(iz,1,is))  &
                                                        *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
                ene(iz,5,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs((mom(iz,5,is)/tau(is)-mom(iz,1,is))  &
                                                        *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
              end do
            end do
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) - twopi * real(iconnect) + gzz(iz),  &
                 real(0.5_DP*fct_e_energy(mxw,gmy,iz) * abs(phi(iz))**2, kind=DP),  &
                 real(0.5_DP*fct_m_energy(mxw,gmy,iz) * abs( Al(iz))**2, kind=DP),  &
                ((ene(iz,imom,is), imom=0,nmom-1), is=0,ns-1)
            end do

                                         !do iz = -global_nz, global_nz-1
                                         !  fct = grootg(iz) / cfsrf
                                         !  wek = wek + fct * fct_e_energy(mxw,gmy,iz) * abs(phi(iz))**2
                                         !  wmk = wmk + fct * fct_m_energy(mxw,gmy,iz) * abs( Al(iz))**2
                                         !end do
                                         !do is = 0, ns-1
                                         !  do imom = 0, nmom-1
                                         !    do iz = -global_nz, global_nz-1
                                         !      fct = grootg(iz) / cfsrf
                                         !      ent(is) = ent(is) + fct * 2._DP * ene(iz,imom,is)
                                         !    end do
                                         !  end do
                                         !end do

          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(gmy) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(gmy)
            call rb_phi_mxmyloop( mxw, gmy, loop, phi )
            call rb_Al_mxmyloop( mxw, gmy, loop, Al )
            do is = 0, ns-1
              do imom = 0, nmom-1
                call rb_mom_mxmyimomisloop( mxw, gmy, imom, is, loop, mom(:,imom,is) )
              end do
              do iz = -global_nz, global_nz-1
                ene(iz,0,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs(mom(iz,0,is)*Znum(is)/fcs(is))**2
                ene(iz,1,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs(mom(iz,1,is)*Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
                ene(iz,2,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * 0.5_DP*abs((2._DP*mom(iz,2,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
                ene(iz,3,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs((mom(iz,3,is)/tau(is)-mom(iz,0,is))*Znum(is)/fcs(is))**2
                ene(iz,4,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * (1._DP/6._DP)*abs((2._DP*mom(iz,4,is)/tau(is)-3._DP*mom(iz,1,is))  &
                                                        *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
                ene(iz,5,is) = 0.5_DP*fcs(is)/Znum(is)*tau(is)  &
                             * abs((mom(iz,5,is)/tau(is)-mom(iz,1,is))  &
                                                        *Znum(is)/fcs(is)/sqrt(tau(is)/Anum(is)))**2
              end do
            end do
            do iz = -global_nz, global_nz-1
              write( omominz, "(99g17.7e3)" ) + twopi * real(iconnect) + gzz(iz),  &
                 real(0.5_DP*fct_e_energy(mxw,gmy,iz) * abs(phi(iz))**2, kind=DP),  &
                 real(0.5_DP*fct_m_energy(mxw,gmy,iz) * abs( Al(iz))**2, kind=DP),  &
                ((ene(iz,imom,is), imom=0,nmom-1), is=0,ns-1)
            end do

                                         !do iz = -global_nz, global_nz-1
                                         !  fct = grootg(iz) / cfsrf
                                         !  wek = wek + fct * fct_e_energy(mxw,gmy,iz) * abs(phi(iz))**2
                                         !  wmk = wmk + fct * fct_m_energy(mxw,gmy,iz) * abs( Al(iz))**2
                                         !end do
                                         !do is = 0, ns-1
                                         !  do imom = 0, nmom-1
                                         !    do iz = -global_nz, global_nz-1
                                         !      fct = grootg(iz) / cfsrf
                                         !      ent(is) = ent(is) + fct * 2._DP * ene(iz,imom,is)
                                         !    end do
                                         !  end do
                                         !end do

          end do
  
      end if

    close( omominz )
                                         !write(80,"(99g17.7e3)") time, wek, wmk, ent(:)
                                         !%%%%%%%%%%%%%%%%%

END SUBROUTINE eneinz_connect


END MODULE out_mominz
