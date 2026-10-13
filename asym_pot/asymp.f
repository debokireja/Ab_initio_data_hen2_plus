c     Long-range tail for He+ + N2 (fixed N2 r), R = Rmatch+1 .. Rend
c     Input : *ci.dat   (R, E1)  -- first line is a header
c     Output: fort.91   (R, V)          -- same energy scale as 15d.dat
       implicit double precision (a-h,o-z)
       parameter (nmax=500)
       dimension rd(nmax),ed(nmax)
       dimension aq(5),aa0(5),aa2(5)
c----- Table 2.4 (C0..C4), a.u.
       data aq /-1.125d0, 0.989d0, 0.100d0,-0.113d0,-0.110d0/
       data aa0/11.561d0, 6.033d0, 0.459d0, 1.660d0,-1.894d0/
       data aa2/ 4.791d0, 6.892d0, 0.479d0, 6.679d0,-5.942d0/
c----- user settings ------------------------------------------------
       gama   = 45.d0      ! angle (deg) of the data in *ci.dat
       sreq   = 2.0811d0   ! N2 equilibrium r (bohr)
       sfix   = 2.0811d0   ! fixed N2 r used in the ab initio data
       rmatch = 12.0d0      ! take the last point of ab intio R where tail is joined to ab initio data
       rend   = 30.d0
       eref   = 0.d0       ! subtract from energies if you want a
c                          ! zero-at-infinity scale (e.g. E(inf))
c--------------------------------------------------------------------
       pi = dacos(-1.d0)
       ct = dcos(gama*pi/180.d0)
       p2 = 0.5d0*(3.d0*ct*ct-1.d0)
       x  = sfix-sreq
c      r-dependent properties (only C0 survives when sfix = sreq)
       al0=0.d0
       al2=0.d0
       qm =0.d0
       do j=1,5
        al0=al0+aa0(j)*x**(j-1)
        al2=al2+aa2(j)*x**(j-1)
        qm =qm +aq(j) *x**(j-1)
       enddo
c----- read ab initio data ----------------------------
       open(10,file='45ci.dat',status='old')
       nd=0
 10    read(10,*,end=20) rr,ee
       nd=nd+1
       rd(nd)=rr
       print*,rr,ee
       ed(nd)=ee-eref
       goto 10
 20    close(10)
c----- ab initio energy at the matching point -----------------------
       im=0
       do i=1,nd
        if (dabs(rd(i)-rmatch).lt.1.d-6) im=i
       enddo
       if (im.eq.0) stop 'rmatch not found in data d.dat'
c----- shift so that the tail is continuous at rmatch ---------------
       df = vmod(rmatch)-ed(im)
       write(*,*) 'shift df =',df
c----- diagnostic: model vs data slope between last two points ------
       if (im.gt.1) then
        dd = ed(im)-ed(im-1)
        dm = vmod(rd(im))-vmod(rd(im-1))
        write(*,'(a,f6.1,a,f6.1)') ' R =',rd(im-1),' ->',rd(im)
        write(*,*) 'ab initio dE =',dd,'  model dE =',dm
       endif
c----- write tail ---------------------------------------------------
       r = rmatch
 30    r = r+1.d0
       if (r.gt.rend+1.d-6) goto 40
       write(91,50) r, vmod(r)-df
 50    format(f7.2,3x,f15.8)
       goto 30
 40    stop
       contains
       double precision function vmod(rr)
       double precision rr
c      charge-quadrupole + charge-induced-dipole (alpha2 = alpha_par -
c      alpha_perp, so alpha(theta) = alpha0 + (2/3)*alpha2*P2)
       vmod = qm*p2/rr**3
     .      - (al0+2.d0*al2*p2/3.d0)/(2.d0*rr**4)
       end function
       end
