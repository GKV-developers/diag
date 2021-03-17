MODULE out_ffinvm
!-------------------------------------------------------------------------------
!
!     Output ff in (v,m)
!                                                   (S. Maeyama, 5 May 2022)
!
!-------------------------------------------------------------------------------
  use diag_header
  implicit none

  private

  public fkinvm_fxv, fkinvm_cnt, fluxinvm_fxv, fluxinvm_cnt


 CONTAINS


SUBROUTINE fkinvm_fxv( mx, gmy, rankz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_fxv_gettime, rb_fxv_mxmyrankzisloop, &
                      rb_phi_gettime, rb_phi_mxmyizloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu, gvp

  integer, intent(in) :: mx, gmy, rankz, is, loop

  real(kind=DP) :: time, time_phi
  complex(kind=DP) :: phi
  complex(kind=DP), dimension(1:2*global_nv,0:global_nm) :: fk
  character(len=4) :: cmx, cmy, crz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: giz, iv, im, loop_phi

    call rb_fxv_gettime( loop, time )
    call rb_fxv_mxmyrankzisloop( mx, gmy, rankz, is, loop, fk )

!- For the analysis of phase -
    loop_phi = nint(time / dtout_ptn)
    giz = -global_nz + rankz*(2*nz)
    call rb_phi_gettime( loop_phi, time_phi )
    call rb_phi_mxmyizloop( mx, gmy, giz, loop_phi, phi )
    if (abs(time_phi - time) > min(dtout_ptn,dtout_fxv)) then
      write(*,*) "# Error: wrong time in fluxinvm_fxv."
      write(*,*) "# time(fxv)=",time,",time(phi)=",time_phi
    end if
!---

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( crz, fmt="(i4.4)" ) rankz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    giz = -global_nz + rankz*(2*nz)
    open( offinvm, file="./data/fkinvm_mx"//cmx//"my"//cmy//"rankz"//crz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(a,g17.7e3)"       ) "# For the analysis of phase at time_phi=",time_phi
      write( offinvm, "(a,g17.7e3,a,g17.7e3)" ) "#   Re[phi]=",real(phi,kind=DP),", Im[phi]=",aimag(phi)
      write( offinvm, "(99a17)" ) "#              vl","mu","vp","Re[ff]","Im[ff]"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), &
                                          real(fk(iv,im), kind=DP), aimag(fk(iv,im))
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fkinvm_fxv


SUBROUTINE fluxinvm_fxv( rankz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output flux in (vl,mu) at a given z, while x,y are averaged.
!                                                   (S. Maeyama, 10 Feb. 2021)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_fxv_gettime, rb_fxv_rankzivimisloop, &
                      rb_phi_gettime, rb_phi_izloop
  use diag_geom, only : gky, gzz, gvl, gmu, gvp, Anum, tau, Znum, gomg, gksq
  use diag_functions, only : besselj0

  integer, intent(in) :: rankz, is, loop

  real(kind=DP) :: time, time_phi, kmo, wr
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: fk, phi
  real(kind=DP), dimension(1:2*global_nv,0:global_nm) :: flux
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:global_nm) :: j0
  character(len=4) :: crz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, giz, iv, im, loop_phi

    call rb_fxv_gettime( loop, time )

    loop_phi = nint(time / dtout_ptn)
    giz = -global_nz + rankz*(2*nz)
    call rb_phi_gettime( loop_phi, time_phi )
    call rb_phi_izloop( giz, loop_phi, phi )
    if (abs(time_phi - time) > min(dtout_ptn,dtout_fxv)) then
      write(*,*) "# Error: wrong time in fluxinvm_fxv."
      write(*,*) "# time(fxv)=",time,",time(phi)=",time_phi
    end if
    
!$OMP parallel do default(none) &
!$OMP shared(gksq,gmu,gomg,tau,Anum,Znum,j0,is,giz) &
!$OMP private(mx,my,im,kmo)
    do im = 0, global_nm
      do my = 0, global_ny
        do mx = -nx, nx
          kmo = sqrt( 2._DP * gksq(mx,my,giz) * gmu(im) / gomg(giz) ) &
              * dsqrt( tau(is)*Anum(is) ) / Znum(is)
          j0(mx,my,im) = besselj0(kmo) ! *  2._DP * pi * vp(giz,im)
        end do
      end do
    end do

    do im = 0, global_nm
      do iv = 1, 2*global_nv
        call rb_fxv_rankzivimisloop( rankz, iv, im, is, loop, fk )
        wr = 0._DP
        do my = 0, global_ny
          do mx = -nx, nx
            wr = wr + real( - ui * gky(my) * j0(mx,my,im) * phi(mx,my)  &
                            * conjg( fk(mx,my) ), kind=DP )
          end do
        end do
        flux(iv,im) = wr
      end do
    end do

    write( crz, fmt="(i4.4)" ) rankz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinvm, file="./data/fluxinvm_rankz"//crz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(99a17)" ) "#              vl","mu","vp","flux"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), flux(iv,im)
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fluxinvm_fxv


SUBROUTINE fkinvm_cnt( mx, gmy, giz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output fk (a mode of ff) in (vl,mu)
!                                                   (S. Maeyama, 5 May 2020)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_mxmyizisloop, &
                      rb_phi_gettime, rb_phi_mxmyizloop
  use diag_geom, only : kx, gky, gzz, gvl, gmu, gvp

  integer, intent(in) :: mx, gmy, giz, is, loop

  real(kind=DP) :: time, time_phi
  complex(kind=DP) :: phi
  complex(kind=DP), dimension(1:2*global_nv,0:global_nm) :: fk
  character(len=4) :: cmx, cmy, ciz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: iv, im, loop_phi

    call rb_cnt_gettime( loop, time)
    call rb_cnt_mxmyizisloop( mx, gmy, giz, is, loop, fk )

!- For the analysis of phase -
    loop_phi = nint(time / dtout_ptn)
    call rb_phi_gettime( loop_phi, time_phi )
    call rb_phi_mxmyizloop( mx, gmy, giz, loop_phi, phi )
    if (abs(time_phi - time) > min(dtout_ptn,dtout_fxv)) then
      write(*,*) "# Error: wrong time in fluxinvm_fxv."
      write(*,*) "# time(fxv)=",time,",time(phi)=",time_phi
    end if
!---

    write( cmx, fmt="(i4.4)" ) mx
    write( cmy, fmt="(i4.4)" ) gmy
    write( ciz, fmt="(i4.4)" ) giz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinvm, file="./data/fkinvm_mx"//cmx//"my"//cmy//"z"//ciz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#   mx=",mx,   ",   kx=",kx(mx)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  gmy=",gmy,  ",  gky=",gky(gmy)
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(a,g17.7e3)"       ) "# For the analysis of phase at time_phi=",time_phi
      write( offinvm, "(a,g17.7e3,a,g17.7e3)" ) "#   Re[phi]=",real(phi,kind=DP),", Im[phi]=",aimag(phi)
      write( offinvm, "(99a17)" ) "#              vl","mu","vp","Re[ff]","Im[ff]"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), &
                                          real(fk(iv,im), kind=DP), aimag(fk(iv,im))
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fkinvm_cnt


SUBROUTINE fluxinvm_cnt( giz, is, loop )
!-------------------------------------------------------------------------------
!
!     Output flux in (vl,mu) at a given z, while x,y are averaged.
!                                                   (S. Maeyama, 10 Feb. 2021)
!
!-------------------------------------------------------------------------------
  use diag_rb, only : rb_cnt_gettime, rb_cnt_izivimisloop, &
                      rb_phi_gettime, rb_phi_izloop
  use diag_geom, only : gky, gzz, gvl, gmu, gvp, Anum, tau, Znum, gomg, gksq
  use diag_functions, only : besselj0

  integer, intent(in) :: giz, is, loop

  real(kind=DP) :: time, time_phi, kmo, wr
  complex(kind=DP), dimension(-nx:nx,0:global_ny) :: fk, phi
  real(kind=DP), dimension(1:2*global_nv,0:global_nm) :: flux
  real(kind=DP), dimension(-nx:nx,0:global_ny,0:global_nm) :: j0
  character(len=4) :: ciz
  character(len=1) :: cis
  character(len=8) :: cloop
  integer :: mx, my, iv, im, loop_phi

    call rb_cnt_gettime( loop, time )

    loop_phi = nint(time / dtout_ptn)
    call rb_phi_gettime( loop_phi, time_phi )
    call rb_phi_izloop( giz, loop_phi, phi )
    if (abs(time_phi - time) > min(dtout_ptn,dtout_fxv)) then
      write(*,*) "# Error: wrong time in fluxinvm_fxv."
      write(*,*) "# time(fxv)=",time,",time(phi)=",time_phi
    end if
    
!$OMP parallel do default(none) &
!$OMP shared(gksq,gmu,gomg,tau,Anum,Znum,j0,is,giz) &
!$OMP private(mx,my,im,kmo)
    do im = 0, global_nm
      do my = 0, global_ny
        do mx = -nx, nx
          kmo = sqrt( 2._DP * gksq(mx,my,giz) * gmu(im) / gomg(giz) ) &
              * dsqrt( tau(is)*Anum(is) ) / Znum(is)
          j0(mx,my,im) = besselj0(kmo) ! *  2._DP * pi * vp(giz,im)
        end do
      end do
    end do

    do im = 0, global_nm
      do iv = 1, 2*global_nv
        call rb_cnt_izivimisloop( giz, iv, im, is, loop, fk )
        wr = 0._DP
        do my = 0, global_ny
          do mx = -nx, nx
            wr = wr + real( - ui * gky(my) * j0(mx,my,im) * phi(mx,my)  &
                            * conjg( fk(mx,my) ), kind=DP )
          end do
        end do
        flux(iv,im) = wr
      end do
    end do

    write( ciz, fmt="(i4.4)" ) giz
    write( cis, fmt="(i1.1)" ) is
    write( cloop, '(i8.8)' ) loop

    open( offinvm, file="./data/fluxinvm_z"//ciz//"s"//cis//"_t"//cloop//".dat" )
      write( offinvm, "(a,i17,a,g17.7e3)" ) "# loop=",loop, ", time=",time
      write( offinvm, "(a,i17,a,g17.7e3)" ) "#  giz=",giz,  ",  gzz=",gzz(giz)
      write( offinvm, "(99a17)" ) "#              vl","mu","vp","flux"
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( offinvm, '(99g17.7e3)' ) gvl(iv), gmu(im), gvp(giz,im), flux(iv,im)
        end do
        write( offinvm, * )
      end do
    close( offinvm )

END SUBROUTINE fluxinvm_cnt


END MODULE out_ffinvm
