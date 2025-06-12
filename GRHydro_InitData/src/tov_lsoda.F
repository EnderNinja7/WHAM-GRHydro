#ifdef HAVE_ODEPACK
c===========================================================
c     History: sode.f
c
c     Driver routine which integrates ODEs defining 
c     TOV star
c 
c    Based on bstar.f by Matthew Choptuik
c===========================================================
      subroutine     tov_lsoda ( rho, p, m, alf,
     &                          r,  nr, tol , n_succ, ssrc) 

      implicit       none


      integer       vsrc,  vsxynt
      real*8         rho2p,         p2rho

      integer        nr
      real*8         rho(nr),p(nr), m(nr), alf(nr)
      real*8         r(nr)
      real*8         tol
      integer        n_succ, ssrc

      character*12    codenm
      parameter    ( codenm = 'tov_lsoda' )

      integer        iargc,         indlnb,         i4arg
      real*8         r8arg

      real*8         r8_never
      parameter    ( r8_never = -1.0d-60 )

c-----------------------------------------------------------
c     Order of system.
c-----------------------------------------------------------
      integer        neq
      parameter    ( neq = 3 )

c-----------------------------------------------------------
c     Storage for solution at requested output radii.
c-----------------------------------------------------------
      integer        maxout
      parameter    ( maxout = 10 000 )

      real*8         y0(neq)
      real*8         vxout(maxout), vyout(maxout,neq),
     &               work(maxout)
      integer        nxout,         ixout,       nxout_succ

      integer        ieq

      logical        ltrace
      parameter    ( ltrace = .false. )

      integer        maxdump
      parameter    ( maxdump = 50 )

      logical        lsodatrace
      parameter    ( lsodatrace = .false. )


c-----------------------------------------------------------
c     LSODA Variables.
c-----------------------------------------------------------
      external       tov_fcn,        tov_jac
      real*8         y(neq)
      real*8         tbgn,       tend
      integer        itol
      real*8         rtol,       atol
      integer        itask,      istate,     iopt
      integer        lrw

c      parameter    ( lrw = 22 + neq * max(16 ,neq + 9) )
      parameter    ( lrw = 22 + neq * (16) )
      real*8         rwork(lrw)

      integer        liw
      parameter    ( liw = 20 + neq )
      integer        iwork(liw)
      integer        jt

      integer        i

c-----------------------------------------------------------
c     Common communication with routine 'tov_fcn' in 'tov_fcn.f' ...
c-----------------------------------------------------------
c   rho0:  Value of rho at r=0
c-----------------------------------------------------------
      real*8  rho0
      common / com_tov / rho0

c-----------------------------------------------------------
c     Parse arguments.  
c-----------------------------------------------------------
      rho0    = rho(1)

      nxout   = nr
      do ixout = 1, nxout
         vxout(ixout) = r(ixout)
      enddo

c-----------------------------------------------------------
c     Initialize those "boundary" conditions which are 
c     fixed by regularity, or which are arbitrary.
c-----------------------------------------------------------
      y0(1) = rho2p(rho0)
      y0(2) = 0.0d0
      y0(3) = 1.0


c-----------------------------------------------------------
c     Echo command line arguments if "local tracing" 
c     enabled ...
c-----------------------------------------------------------
      if( ltrace ) then
         write(0,*) 'rho0: ', rho0
         write(0,*) 'p(1): ', y0(1)
         write(0,*) 'm(1): ', y0(2)
         write(0,*) 'alf(1): ', y0(3)
      end if

c-----------------------------------------------------------
c     Get output radii 
c-----------------------------------------------------------
      if( ltrace ) then
         if( nxout .le. maxdump ) then
            call dvdump(vxout,nxout,codenm//': output radii',0)
         else 
            call dvdmp1(vxout,work,
     &                  int(1.0d0 * nxout / maxdump + 0.5d0),
     &                  nxout,codenm//': selected output radii',0)
         end if
         write(0,*)
         write(0,*) codenm//': Initial time: ', vxout(1)
         write(0,*)
      end if

c-----------------------------------------------------------
c     Set LSODA parameters ...
c-----------------------------------------------------------
      itol   = 1
      rtol   = tol
      atol   = tol
      itask  = 1 
      iopt   = 0
      jt     = 2

c-----------------------------------------------------------
c     Initialize the solution ...
c-----------------------------------------------------------
      do ieq = 1 , neq
         y(ieq)       = y0(ieq)
         vyout(1,ieq) = y0(ieq)
      end do

c-----------------------------------------------------------
c     Do the integration ...
c-----------------------------------------------------------
      ssrc = 0
      do ixout = 2 , nxout
         istate = 1
c-----------------------------------------------------------
c        Need these temporaries since lsoda overwrites 
c        tend ... 
c-----------------------------------------------------------
100      continue
         tbgn = vxout(ixout-1)
         tend = vxout(ixout)
         call lsoda(tov_fcn,neq,y,tbgn,tend,
     &           itol,rtol,atol,itask,    
     &           istate,iopt,rwork,lrw,iwork,liw,tov_jac,jt)
         if( lsodatrace ) then
            write(0,1000) codenm, ixout, vxout(ixout), 
     &      vxout(ixout+1)
1000        format(' ',a,': Step ',i4,'  t = ',1pe10.3,
     &      ' .. ',1pe10.3)
            write(0,*) codenm//': lsoda reurns ', istate
         end if

c         if (istate .eq. -2) then
c            write(0,*) 'istate = -2, trying again...'
c            itol   = 1
c            rtol   = tol
c            atol   = tol
c            itask  = 1
c            iopt   = 0
c            jt     = 2
c            istate = 3
c            goto 100
c         endif

         if( istate .lt. -1 ) then
            write(0,1500) codenm, istate, ixout, nxout, 
     &                    vxout(ixout-1), vxout(ixout)
1500        format(/' ',a,': Error return ',i2,
     &              ' from integrator LSODA.'/
     &              '       At output time ',i5,' of ',i5/
     &              '       Interval ',1pe11.3,' .. ',
     &             1pe11.3/)
             nxout_succ = ixout - 1
             ssrc = 1
            go to 500

         end if

         do ieq = 1 , neq
            vyout(ixout,ieq) = y(ieq)
         end do

      end do
      nxout_succ = nxout

500   continue

c-----------------------------------------------------------
c     OUTPUT
c-----------------------------------------------------------
      do ixout = 1 , nxout_succ
         p(ixout)    = vyout(ixout, 1)
         m(ixout)    = vyout(ixout, 2)
         alf(ixout)    = vyout(ixout, 3)
         rho(ixout)    = p2rho(p(ixout))
      end do

      n_succ = nxout_succ

      return

      end





c===========================================================
c     Some double precision vector utility routines. 
c
c     Originally from ~matt/vutil/dveclib.f
c===========================================================
 
c-----------------------------------------------------------
c     Dumps vector labelled with LABEL on UNIT.
c-----------------------------------------------------------
      subroutine dvdump(v,n,label,unit)
         implicit       none

         real*8         v(*)
         character*(*)  label
         integer        i, n, st, unit
 
         if( n .lt. 1 ) go to 130
            write(unit,100) label
 100        format(/' <<< ',A,' >>>'/)
            st = 1
 110        continue
               write(unit,120) ( v(i) , i = st , min(st+3,n))
 120           FORMAT(' ',4(1PE19.10))
               st = st + 4
            if( st .le. n ) go to 110
 
 130     continue
 
         return
      end
 
c-----------------------------------------------------------
c     Extension of dvdump which dumps every 'inc'th element.
c-----------------------------------------------------------
      subroutine dvdmp1(v,w,inc,n,label,unit)
         implicit       none

         real*8         v(*),    w(*)
         character*(*)  label
         integer        inc,     n,      unit
 
         call dvinj(v,w,inc,n)
         call dvdump(w,1+(n-1)/inc,label,unit)

         return

      end
 
c---------------------------------------------------------
c     Injects every incth element of v1 into v2
c---------------------------------------------------------
      subroutine dvinj(v1,v2,inc,n)
         implicit    none
   
         real*8      v1(*),     v2(*)
         integer     i,         inc,     j,      n
 
         j = 1
         do i = 1 , n , inc
            v2(j) = v1(i)
            j = j + 1
         end do

 
         return
      end

#else
      subroutine tov_lsoda()
      return 
      end
 
#endif
