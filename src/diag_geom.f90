MODULE diag_geom
!-------------------------------------------------------------------------------
!
!     Module for geometry
!                                                   (S. Maeyama, 1 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_header
!- for VMEC -
  use GKV_vmecbzx, only: vmecbzx_boozx_read, vmecbzx_boozx_coeff
!------------
!- for IGS -
  use GKV_igs,    only: igs_read, igs_coeff
!-----------
  implicit none

  public

  real(kind=DP), dimension(0:2*nxw-1)              :: xx
  real(kind=DP), dimension(0:2*nyw-1)              :: yy
  real(kind=DP), dimension(-nx:nx)                 :: kx
  real(kind=DP), dimension(0:global_ny)            :: gky
  real(kind=DP), dimension(-global_nz:global_nz-1) :: gzz, gomg, grootg
  real(kind=DP), dimension(1:2*global_nv)          :: gvl
  real(kind=DP), dimension(0:global_nm)            :: gmu
  real(kind=DP), dimension(-global_nz:global_nz-1,0:global_nm) :: gvp
  real(kind=DP), dimension(-global_nz:global_nz-1) :: dvp
  real(kind=DP),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) :: gfmx
  real(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gksq,  &
                                                            fct_e_energy,  &
                                                            fct_m_energy
  real(kind=DP),  &
    dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: g0, g1

  complex(kind=DP), dimension(0:global_ny)         :: ck
  integer, dimension(0:global_ny)                  :: dj

  real(kind=DP), dimension(0:ns-1) ::   R0_Ln,  &    ! R0/Lns
                                        R0_Lt,  &    ! R0/Lts
                                           nu,  &    ! collision freq.   
                                         Anum,  &    ! mass number
                                         Znum,  &    ! charge number     
                                          fcs,  &    ! charge-density fraction 
                                          sgn,  &    ! signs of charge   
                                          tau,  &    ! T-ratio
                                         dns1
  real(kind=DP) :: dv, cfsrf, lambda_i, beta, q_0, q_bar, tau_ad, vmax
  character(len=8)  :: equib_type
  real(kind=DP) :: r_minor = 84._DP, s_hat
  real(kind=DP) :: eps_r
  real(kind=DP) :: rdeps00, eps_hor, lprd, mprd, lmmq, malpha
  real(kind=DP) :: eps_mor, eps_por, lprdm1, lprdp1, lmmqm1, lmmqp1
  real(kind=DP) :: eps_rnew, rdeps1_0, rdeps1_10, rdeps2_10, rdeps3_10
  real(kind=DP) :: lx, ly, lz, kxmin, kymin, dz, mmax, dm, del_c
  integer       :: n_tht, m_j, ibprime, nx0


 CONTAINS


SUBROUTINE geom_init( inum )
!-------------------------------------------------------------------------------
!
!     Set geometry
!                                                   (S. Maeyama, 1 March 2014)
!
!-------------------------------------------------------------------------------
  use diag_functions, only : besseli0, besseli1
  
  integer, intent(in) :: inum

  character(len=3) :: cnum
  real(kind=DP), dimension(1:3,1:3) :: gg
  real(kind=DP) :: bb
  integer :: mx, my, iz, iv, im, is
  integer :: ny_size, nwk, rankw
  real(kind=DP) :: lz_l, domgdz, domgdx, domgdy, theta, phi_ax
  real(kind=DP) :: s_input, s_0      ! radial label of fluxtube center 
  integer       :: isw, nss, ntheta, nzeta
!- for VMEC -
  namelist /vmecp/ s_input, nss, ntheta, nzeta
!------------
!- for IGS -
  integer       :: mc_type           ! 0:Axisym., 1:Boozer, 2:Hamada
  integer       :: q_type            ! 0:use q and s_hat value in confp, 1:calclated by IGS
  namelist /igsp/ s_input, mc_type, q_type, nss, ntheta
!-----------

  namelist /equib/ equib_type
  namelist /physp/ R0_Ln,  &    ! R0/Lns
                   R0_Lt,  &    ! R0/Lts
                      nu,  &    ! collision freq. 
                    Anum,  &    ! mass number
                    Znum,  &    ! charge number 
                     fcs,  &    ! charge-density fraction 
                     sgn,  &    ! signs of charge 
                     tau,  &    ! T-ratio Ts/T0, T0=reference temp.
                    dns1,  &    ! Initial value
                  tau_ad,  &    ! Temperature ratio in adiabatic model
                lambda_i,  &    ! Debye/rho 
                    beta,  &    ! mu0*ni*Ti/B^2
                 ibprime,  &    ! Grad-P term in magnetic drift
                    vmax,  &    ! Velocity space lenge
                     nx0        ! Initial wave number range

  namelist /nperi/ n_tht, kymin, m_j, del_c
  namelist /confp/ eps_r, eps_rnew,                       &
                   q_0, s_hat,                            &
                   lprd, mprd, eps_hor, eps_mor, eps_por, &
                   rdeps00, rdeps1_0, rdeps1_10,          & 
                   rdeps2_10, rdeps3_10, malpha

  namelist /times/ tend, dtout_fxv, dtout_ptn, dtout_eng, dtout_dtc


! ---- set y range --------------------------
    ny_size = global_ny + 1 
    if( mod(ny_size,nprocw) == 0 )  then
      nwk    = ny_size / nprocw
    else
      nwk    = ny_size / nprocw + 1
    endif
    do rankw = 0, nprocw-1
      !--- global index range ---------------- 
      ist_y_g(rankw)  = nwk*rankw
      iend_y_g(rankw) = min( nwk*(rankw+1)-1, (ny_size-1) )
      nsize_y(rankw)  = iend_y_g(rankw) - ist_y_g(rankw) + 1
      !--- local index range ---------------- 
      ist_y(rankw)    = 0
      iend_y(rankw)   = iend_y_g(rankw) - ist_y_g(rankw)
  
      if( rankw == 0 )   then
         ist1_y(rankw)    = 1
      else 
         ist1_y(rankw)    = 0
      endif
    end do

    write( cnum, fmt="(i3.3)" ) inum

    open( inml, file="../gkvp_namelist."//cnum, status="old", action="read" )

    read(inml,nml=equib)

    read(inml,nml=times)

    read(inml,nml=physp)

    read(inml,nml=nperi)

    if( trim(equib_type) == "slab" ) then

      read(inml,nml=confp)

      lprdm1   = 0._DP
      lprdp1   = 0._DP
      lmmq     = 0._DP
      lmmqm1   = 0._DP
      lmmqp1   = 0._DP
      q_0      = 1._DP ! For now, fixed q_0=1. Changing q_0 can extend parallel z-box size.
      s_hat    = 0._DP ! only shear less slab
      eps_r    = 1._DP
      eps_hor  = 0._DP
      lprd     = 0._DP
      mprd     = 0._DP
      malpha   = 0._DP
      rdeps00  = 0._DP
      eps_mor  = 0._DP
      eps_por  = 0._DP

    else if( trim(equib_type) == "analytic"  .or.  &
             trim(equib_type) == "s-alpha"   .or.  &
             trim(equib_type) == "circ-MHD" ) then

      read(inml,nml=confp)

      lprdm1   = lprd - 1.0_DP
      lprdp1   = lprd + 1.0_DP

      lmmq     = lprd   - mprd * q_0
      lmmqm1   = lprdm1 - mprd * q_0
      lmmqp1   = lprdp1 - mprd * q_0

    else if( trim(equib_type) == "vmec" ) then

      lz_l     = real( n_tht, kind=DP ) * pi / real( nprocz, kind=DP )       ! local z-length

      read(inml,nml=confp)

      read(inml,nml=vmecp)

      call vmecbzx_boozx_read( nss, ntheta, nzeta )

       isw = 0
       iz = 0
       call vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input, iz, 0._DP,  lz_l,   &  ! input 
                         s_0,           q_0,     s_hat,    eps_r,  phi_ax,             &  ! output
                         gomg(iz), grootg(iz),    domgdx,   domgdz,  domgdy,             &
                         gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                      &
                         gg(2,3),   gg(3,3)  )

    else if( trim(equib_type) == "eqdsk" ) then

      lz_l     = real( n_tht, kind=DP ) * pi / real( nprocz, kind=DP )       ! local z-length

      read(inml,nml=confp)
      read(inml,nml=igsp)

      call igs_read( mc_type, nss, ntheta )

      if ( q_type == 1 ) then
        isw = 0
        iz = 0
        call igs_coeff( isw,  mc_type,   nss,    ntheta,  s_input,  0._DP, lz_l, &  ! input 
                        s_0,       q_0,     s_hat,    eps_r,  theta,             &  ! output
                        gomg(iz), grootg(iz),    domgdx,   domgdz, domgdy,       &
                        gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                 &
                        gg(2,3),   gg(3,3)  )
      end if

    else

      write( olog, * ) " # Currently, this equilibrium is not available."
      equib_type = "circ-MHD"
      read(inml,nml=confp)
      !stop

    end if

    !kxmin    = 2._DP * pi * s_hat * kymin / real(m_j, kind=DP)
    if (abs(s_hat) < 1.d-10) then ! When s_hat == ZERO
      m_j = 0
      kxmin = kymin
    else if (m_j == 0) then
      kxmin = kymin
    else
      kxmin    = abs(2._DP * pi * s_hat * kymin / real(m_j, kind=DP))
    end if
    lx       = pi / kxmin
    ly       = pi / kymin
    lz       = real( n_tht, kind=DP ) * pi        ! total z-length

    dz       = lz / real( global_nz, kind=DP )

!      vmax     = 5._DP
    dv       = 2._DP * vmax / real( 2 * nv * nprocv -1, kind=DP )

    mmax     = vmax
    dm       = mmax / real( nprocm * ( nm+1 ) - 1, kind=DP )

    do mx = 0, 2*nxw-1
      xx(mx) = -lx + lx/real(nxw)*mx
    end do

    do my = 0, 2*nyw-1
      yy(my) = -ly + ly/real(nyw)*my
    end do

    do mx = -nx, nx
      kx(mx) = kxmin * real( mx, kind=DP )
    end do

    do my = 0, global_ny
      gky(my) = kymin * real( my, kind=DP )
    end do

    do iz = -global_nz, global_nz-1
      gzz(iz) = dz * real( iz, kind=DP )
      if( trim(equib_type) == "slab" ) then

        gomg(iz) = 1._DP
        grootg(iz) = q_0
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = kx(mx)**2 + gky(my)**2
          end do
        end do

      else if( trim(equib_type) == "analytic" ) then

        gomg(iz) = 1._DP                                              &
                 - eps_r * ( cos( gzz(iz) )                           &
                         + eps_hor * cos( lmmq   * gzz(iz) - malpha ) &
                         + eps_mor * cos( lmmqm1 * gzz(iz) - malpha ) &
                         + eps_por * cos( lmmqp1 * gzz(iz) - malpha ) )
        grootg(iz) = q_0 / gomg(iz)
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = ( kx(mx) + s_hat * gzz(iz) * gky(my) )**2 + gky(my)**2
          end do
        end do

      else if( trim(equib_type) == "s-alpha"  ) then

        gomg(iz)   = 1._DP - eps_r * cos( gzz(iz) )        ! s-alpha with eps-expansion
        !gomg(iz)   = 1._DP / ( 1._DP + eps_r * cos( gzz(iz) ) ) ! for benchmark
        grootg(iz) = q_0 / gomg(iz)
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = ( kx(mx) + s_hat * gzz(iz) * gky(my) )**2 + gky(my)**2
          end do
        end do

      else if( trim(equib_type) == "circ-MHD" ) then

        q_bar = dsqrt( 1._DP - eps_r**2 )*q_0
        theta = 2._DP*atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) &
                                           * tan(gzz(iz)/2._DP) )
        gomg(iz)   = sqrt( q_bar**2 + eps_r**2 ) &
                   / ( 1._DP + eps_r*cos( theta ) ) / q_bar
        grootg(iz) = (q_0**2/q_bar)*( 1._DP+eps_r*cos(theta) )**2

        gg(1,1) = 1._DP
        gg(1,2) = s_hat*gzz(iz) - eps_r*sin(gzz(iz))/(1._DP-eps_r**2)
        gg(2,2) = (s_hat*gzz(iz))**2 - 2._DP*s_hat*gzz(iz)*eps_r*sin(gzz(iz))/(1._DP-eps_r**2) &
                    + (q_bar**2+eps_r**2)/((1._DP+eps_r*cos(theta))**2)/(q_0**2)               &
                    + (eps_r*sin(gzz(iz)))**2/(1._DP-eps_r**2)**2
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)          &
                           + 2._DP*kx(mx)*gky(my)*gg(1,2) &
                           + (gky(my)**2)*gg(2,2)
          end do
        end do

      else if( trim(equib_type) == "vmec" ) then

        isw = 1
        call vmecbzx_boozx_coeff( isw,  nss,  ntheta,  nzeta,  s_input, iz, gzz(iz),  lz_l,   &  ! input 
                                  s_0,       q_0,     s_hat,    eps_r, phi_ax,          &  ! output
                               gomg(iz), grootg(iz),    domgdx,   domgdz, domgdy,          &
                               gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),                  &
                               gg(2,3),   gg(3,3)  )
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)          &
                           + 2._DP*kx(mx)*gky(my)*gg(1,2) &
                           + (gky(my)**2)*gg(2,2)
          end do
        end do

      else if( trim(equib_type) == "eqdsk" ) then

        isw = 1
        call igs_coeff( isw,  mc_type,   nss,    ntheta,  s_input, gzz(iz), lz_l, &  ! input 
                              s_0,       q_0,     s_hat,    eps_r,   theta,       &  ! output
                              gomg(iz), grootg(iz),   domgdx,  domgdz, domgdy,    &
                              gg(1,1),   gg(1,2),   gg(1,3),  gg(2,2),            &
                              gg(2,3),   gg(3,3)  )
        do my = 0, global_ny
          do mx = -nx, nx
            gksq(mx,my,iz) = (kx(mx)**2)*gg(1,1)          &
                           + 2._DP*kx(mx)*gky(my)*gg(1,2) &
                           + (gky(my)**2)*gg(2,2)
          end do
        end do

      end if
    end do

    cfsrf = 0._DP
    do iz = -global_nz, global_nz-1
      cfsrf = cfsrf + grootg(iz)
    end do

    do iv = 1, 2*global_nv
      gvl(iv) = - vmax + dv * real( iv-1, kind=DP )
    end do

    do im = 0, global_nm
      gmu(im) = 0.5_DP * ( dm * real( im, kind=DP ) )**2
    end do
    do im = 0, global_nm
      do iz = -global_nz, global_nz-1
        gvp(iz,im) = sqrt( 2._DP * gmu(im) * gomg(iz) )
      end do
    end do
    do iz = -global_nz, global_nz-1
      dvp(iz) = gvp(iz,1)
    end do

    do my = 0, global_ny
      ck(my)   = exp( ui * 2._DP * pi * del_c &
                         * real( n_tht * my, kind=DP ) )
      dj(my)   = - m_j * n_tht * my
                              !  del_c = q_0*n_alp-int(q_0-n_alp)
                              !  m_j   = 2*n_alp*q_d
    end do

    do im = 0, global_nm
      do iv = 1, 2*global_nv
        do iz = -global_nz, global_nz-1
          gfmx(iz,iv,im) = exp( - 0.5_DP * gvl(iv)**2 - gomg(iz) * gmu(im) ) &
                         / sqrt( twopi**3 )
        end do
      end do
    end do

    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            bb = gksq(mx,my,iz)*tau(is)*Anum(is)/(Znum(is)**2*gomg(iz)**2)
            if ( bb < 150._DP ) then
              g0(mx,my,iz,is) = besseli0(bb) * exp(-bb)
            else
              g0(mx,my,iz,is) = ( 1.d0                                    &
                                + 0.25d0             / ( 2.d0 * bb )      &
                                + 9.d0 / 32.d0       / ( 2.d0 * bb )**2   &
                                + 75.d0 / 128.d0     / ( 2.d0 * bb )**3   &
                                + 3675.d0 / 2048.d0  / ( 2.d0 * bb )**4   &
                                + 59535.d0 / 8192.d0 / ( 2.d0 * bb )**5 ) &
                                / sqrt( 2.d0 * pi * bb )
            end if
            if ( bb < 150._DP ) then
              g1(mx,my,iz,is) = besseli1(bb) * exp(-bb)
            else
              g1(mx,my,iz,is) = ( 1.d0                                    &
                                - 0.75d0             / ( 2.d0 * bb )      &
                                - 15.d0 / 32.d0      / ( 2.d0 * bb )**2   &
                                - 105.d0 / 128.d0    / ( 2.d0 * bb )**3   &
                                - 4725.d0 / 2048.d0  / ( 2.d0 * bb )**4   &
                                - 72765.d0 / 8192.d0 / ( 2.d0 * bb )**5 ) &
                                / sqrt( 2.d0 * pi * bb )
            end if
          end do
        end do
      end do
    end do

    do iz = -global_nz, global_nz-1
      do my = 0, global_ny
        do mx = -nx, nx
          fct_e_energy(mx,my,iz) = lambda_i * gksq(mx,my,iz)
        end do
      end do
    end do
    do is = 0, ns-1
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            fct_e_energy(mx,my,iz) = fct_e_energy(mx,my,iz)  &
                         + Znum(is) * fcs(is) / tau(is) * ( 1._DP - g0(mx,my,iz,is) )
          end do
        end do
        fct_e_energy(0,0,iz) = 0._DP
      end do
    end do

    if ( beta > eps ) then
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            fct_m_energy(mx,my,iz) = gksq(mx,my,iz) / beta
          end do
        end do
        fct_m_energy(0,0,iz) = 0._DP
      end do
    else
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            fct_m_energy(mx,my,iz) = 0._DP
          end do
        end do
      end do
    end if

    close(inml)


END SUBROUTINE geom_init


END MODULE diag_geom
