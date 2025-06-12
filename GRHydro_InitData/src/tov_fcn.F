#ifdef HAVE_ODEPACK
c------------------------------------------------------------------------
c    equation of state, converts p to rho
c      currently using a polytrope
c------------------------------------------------------------------------
      real*8  function p2rho(p)
      implicit none
      real*8 p

      p2rho = p ** (3/4.0d0)
  
      return 
      end

c------------------------------------------------------------------------
c    equation of state, converts rho to p
c      currently using a polytrope
c------------------------------------------------------------------------
      real*8  function rho2p(rho)
      implicit none
      real*8 rho

      rho2p = rho ** (4/3.0d0)
 
      return
      end


c-----------------------------------------------------------
c     Integrates static equations of motion for spherically 
c     symmetric, general-relativistic TOV star:
c
c     dp/dr = - (p + rho)*(m+4*pi*r**3*p)/ r / (r-2*m)
c     dm/dr = 4*pi*r**2 * rho
c     d alpha/dr = alpha*(m+4*pi*r**3*rho )/ r / (r-2*m)
c
c  (see, e.g. Walds textbook pp126-127.)
c  where g_tt = - alpha**2, and g_rr = 1/(1-2m/r)
c
c  Note that (my rho) = (the-rest-of-GRHydros rho)*(1+epsilon)
c
c-----------------------------------------------------------
c     neq = 3
c-----------------------------------------------------------

      subroutine tov_fcn(neq,r,y,yprime)
         implicit    none

         real*8   rho0
         common / com_tov / rho0

         real*8      p2rho

         integer     neq
         real*8      r,     y(neq),    yprime(neq)

         real*8      pi
         parameter   ( pi = 3.1415926 5358 9793 d0 )

         real*8      p, m, alf, rho, fourpr2

         p = y(1)
         m = y(2)
         alf  = y(3)

         if( r .eq. 0.0d0 ) then
            yprime(1) = 0.0d0
            yprime(2) = 0.0d0
            yprime(3) = 0.0d0
         else
            fourpr2 = 4*pi*r*r
c  BUG: should be using call to CactusEOS           
            rho = p2rho(p)

            yprime(1) = - (p + rho)*(m+fourpr2*r*p)/ r / (r-2*m)

            yprime(2) = fourpr2 * rho

            yprime(3) = alf*(m+fourpr2*r*rho)/ r / (r-2*m)
          endif

         return
      end 

c-----------------------------------------------------------
c     Implements Jacobian  (optional) ...
c   Does absolutely nothing
c-----------------------------------------------------------
      subroutine tov_jac
         implicit    none

         real*8   rho0
         common / com_tov / rho0

         return
      end
#else
      subroutine tov_fcn()
      return 
      end
#endif
