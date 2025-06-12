 /*@@
   @file      GRHydro_Bondi.c 
   @date      Wed Jan 13 13:00:49 EST 2010
   @author    Scott C. Noble
   @desc 
   Hydro initial data for the relativistic Bondi solution about 
   a single Schwarzschild black hole. 
   @enddesc 
 @@*/

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************
   Calculates the Bondi solution, or the spherically symmetric hydrostationary 
   solution to a fluid on a static fixed background spacetime. We assume that one can
   calculate a radius "r" from the grid and that with respect to this radial coordinate, 
   the solution satisfies 

   d (\rho u^r) / dr    = 0 

   Assumes that the equation of state is  P = K \rho^\Gamma   and K  is set by 
   the location of the sonic point.     


 -- Implicitly assumes that there is no spin in the geometry as there is no Bondi 
    solution for spinning black holes.  If a spin is specified, a spherically symmetric 
    is still assumed but the 4-velocity is set consistently with the spinning spacetime. 

***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// set to 1 to tracing output
#define LTRACE   1


/* Mnemonics: */

#define NDIM (4)  /* Number of spacetime dimensions */

/* mnemonics for dimensional indices */
#define TT      (0)     
#define RR      (1)
#define TH      (2)
#define PH      (3)

/* mnemonics for dimensional indices */
#define XX      (1)
#define YY      (2)
#define ZZ      (3)

/* mnemonics for coordinate system choice */
#define COORD_BOYERLINDQUIST (0)
#define COORD_KERRSCHILD     (1)
#define COORD_ISOTROPIC      (2)

/*  Macros:  */
#define DLOOP1     for(i=0  ;i<NDIM ;i++)
#define DLOOP2     for(i=0  ;i<NDIM ;i++) for(j=0  ;j<NDIM ;j++)

#define sincos sincos1   // Avoid name conflict with libc function
static inline void sincos(CCTK_REAL theta,
                          CCTK_REAL *restrict sth, CCTK_REAL *restrict cth)
{
  *sth = sin(theta);
  *cth = cos(theta);
}

//Newton-Raphson parameters:
#define NEWT_DIM_B        (1      )
#define MAX_NEWT_ITER_B   (30     )    /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL_B        (1.0e-15)    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL_B    (1.0e-10)    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER_B (2      )
#define SMALL_BONDI (1.e-20)

//some math helpers
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))


static CCTK_REAL  Mdot, rs, vs_sq, vs, cs_sq, cs, rhos, hs, K, Qdot, gamma_eos, r_sol;
static CCTK_REAL  M, a, asq, Msq; 
static CCTK_REAL  x_cen, y_cen, z_cen;
static unsigned short int coord_type; 


/******************************************************************************
  bl_to_ks_con():
  ----------
       -- transforms a contravariant vector in BL coordinates to KS coordinates;
 ******************************************************************************/
static void bl_to_ks_con(CCTK_REAL *x, CCTK_REAL blcon[], CCTK_REAL kscon[] )
{
  CCTK_REAL delta_m1, dtKS_drBL, dphiKS_drBL;
  
  delta_m1 = 1./( x[RR] * ( x[RR] - 2*M ) + asq );
  dtKS_drBL   = 2*M*x[RR] * delta_m1; 
  dphiKS_drBL = a * delta_m1;

  kscon[TT] = blcon[TT]  +  blcon[RR] * dtKS_drBL;
  kscon[RR] = blcon[RR];
  kscon[TH] = blcon[TH];
  kscon[PH] = blcon[PH]  +  blcon[RR] * dphiKS_drBL;
  
  return;
}

/******************************************************************************
  dxc_dxs_ks_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^\hat{a}_a  defined:

      x^\hat{a}[Cartesian]  = \Lambda^\hat{a}_a  x^a[Spherical] 

            where dxc_dxs[i][j] = \Lambda^i_j 

        for  Kerr-Schild coordinates.

 ******************************************************************************/
static void dxc_dxs_ks_calc(CCTK_REAL *x_cart, CCTK_REAL *x_spher, CCTK_REAL dxc_dxs[NDIM][NDIM] )
{
  int i, j; 
  CCTK_REAL r, th, ph; 
  CCTK_REAL sth,cth,sph,cph;

  for( i = 0 ; i < NDIM ; i++ ) {
    for( j = 0 ; j < NDIM ; j++ ) {
      dxc_dxs[j][i] = 0. ;
    } 
  }
  
  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] = cph * sth             ;  // dx/dr
  dxc_dxs[XX][TH] = (r*cph - a*sph) * cth ;  // dx/dtheta
  dxc_dxs[XX][PH] = -(x_cart[YY] - y_cen) ;  // dx/dphi
  dxc_dxs[YY][RR] = sph * sth             ;  // dy/dr
  dxc_dxs[YY][TH] = (r*sph + a*cph) * cth ;  // dy/dtheta
  dxc_dxs[YY][PH] = x_cart[XX] - x_cen    ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r*sth                ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi
  
  return;

}

/******************************************************************************
  dxc_dxs_bl_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^\hat{a}_a  defined:

      x^\hat{a}[Cartesian]  = \Lambda^\hat{a}_a  x^a[Spherical] 

            where dxc_dxs[i][j] = \Lambda^i_j 

        for  Boyer-Lindquist coordinates.

 ******************************************************************************/
static void dxc_dxs_bl_calc(CCTK_REAL *x_cart, CCTK_REAL *x_spher, CCTK_REAL dxc_dxs[NDIM][NDIM] )
{
  int i, j; 
  CCTK_REAL r, th, ph, rterm, dr; 
  CCTK_REAL sth,cth,sph,cph;

  for( i = 0 ; i < NDIM ; i++ ) {
    for( j = 0 ; j < NDIM ; j++ ) {
      dxc_dxs[j][i] = 0. ;
    } 
  }
  
  r  = x_spher[RR]; 
  th = x_spher[TH];
  ph = x_spher[PH]; 

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  rterm = sqrt( r*r + asq );
  dr = r / rterm; 
  
  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] = cph * sth * dr        ;  // dx/dr
  dxc_dxs[XX][TH] = rterm * cph * cth     ;  // dx/dtheta
  dxc_dxs[XX][PH] = -(x_cart[YY] - y_cen) ;  // dx/dphi
  dxc_dxs[YY][RR] = sph * sth * dr        ;  // dy/dr
  dxc_dxs[YY][TH] = rterm * sph * cth     ;  // dy/dtheta
  dxc_dxs[YY][PH] = x_cart[XX] - x_cen    ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r*sth                ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi
  
  return;

}

/******************************************************************************
  dxc_dxs_iso_calc():
  ----------------------
       -- calculates the transformation matrix  Lambda^\hat{a}_a  defined:

      x^\hat{a}[Cartesian]  = \Lambda^\hat{a}_a  x^a[Spherical] 

            where dxc_dxs[i][j] = \Lambda^i_j 

        for  "Isotropic" coordinates.

 ******************************************************************************/
static void dxc_dxs_iso_calc(CCTK_REAL *x_cart, CCTK_REAL *x_spher, CCTK_REAL dxc_dxs[NDIM][NDIM] )
{
  int i, j; 
  CCTK_REAL th, ph,r_iso;
  CCTK_REAL sth,cth,sph,cph;

  for( i = 0 ; i < NDIM ; i++ ) {
    for( j = 0 ; j < NDIM ; j++ ) {
      dxc_dxs[j][i] = 0. ;
    } 
  }

  /* BL spherical coordinates : */
  th = x_spher[TH];
  ph = x_spher[PH]; 

  r_iso = x_spher[TT];

  sincos( th, &sth , &cth ); 
  sincos( ph, &sph , &cph ); 

  dxc_dxs[TT][TT] = 1.                    ;  // dt/dt
  dxc_dxs[XX][RR] = cph * sth             ;  // dx/dr
  dxc_dxs[XX][TH] = r_iso * cph * cth     ;  // dx/dtheta
  dxc_dxs[XX][PH] = -(x_cart[YY] - y_cen) ;  // dx/dphi
  dxc_dxs[YY][RR] = sph * sth             ;  // dy/dr
  dxc_dxs[YY][TH] = r_iso * sph * cth     ;  // dy/dtheta
  dxc_dxs[YY][PH] = x_cart[XX] - x_cen    ;  // dy/dphi
  dxc_dxs[ZZ][RR] = cth                   ;  // dz/dr
  dxc_dxs[ZZ][TH] = -r_iso*sth            ;  // dz/dtheta
  //  dxc_dxs[ZZ][PH] = 0.                    ;  // dz/dphi

  return;

}

/******************************************************************************
  ks_cart_to_ks_spher_pos():
  ----------
   -- transforms the position in Cartesian KS coordinates to Spherical KS coordinates;
 ******************************************************************************/
static void ks_cart_to_ks_spher_pos(CCTK_REAL *x_cart, CCTK_REAL *x_spher)
{
  CCTK_REAL xx,yy,zz,r,t3;

  xx = x_cart[XX] - x_cen; 
  yy = x_cart[YY] - y_cen; 
  zz = x_cart[ZZ] - z_cen; 

  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  r = sqrt(   t3 + sqrt( t3*t3 + asq*zz*zz )  );

  x_spher[TT] = x_cart[TT];
  x_spher[RR] = r; 
  x_spher[TH] = acos(zz / r); 
  t3 = atan2( (yy*r - xx*a) , (xx*r + yy*a) ); 
  if( t3 < 0. ) {  t3 += 2.*M_PI; } 
  x_spher[PH] = t3;

  return;
}

/******************************************************************************
  bl_cart_to_bl_spher_pos():
  ----------
   -- transforms the position in Cartesian BL coordinates to Spherical BL coordinates;
 ******************************************************************************/
static void bl_cart_to_bl_spher_pos(CCTK_REAL *x_cart, CCTK_REAL *x_spher)
{
  CCTK_REAL xx,yy,zz,r,t3;

  xx = x_cart[XX] - x_cen; 
  yy = x_cart[YY] - y_cen; 
  zz = x_cart[ZZ] - z_cen; 

  t3 = 0.5*(xx*xx + yy*yy + zz*zz - asq); 
  r = sqrt(   t3 + sqrt( t3*t3 + asq*zz*zz )  );

  x_spher[TT] = x_cart[TT];
  x_spher[RR] = r; 
  x_spher[TH] = acos(zz / r); 
  t3 = atan2( yy , xx ); 
  if( t3 < 0. ) {  t3 += 2.*M_PI; } 
  x_spher[PH] = t3;

  return;
}

/******************************************************************************
  iso_cart_to_bl_spher_pos():
  ----------
   -- transforms the position in Cartesian "Isotropic" coordinates to Spherical BL coordinates;
   -- r^2  =  x^2 + y^2 + z^2 in these coordinates, so they are slightly distorted from flatspace
 ******************************************************************************/
static void iso_cart_to_bl_spher_pos(CCTK_REAL *x_cart, CCTK_REAL *x_spher)
{
  CCTK_REAL xx,yy,zz,riso,r,t3;

  xx = x_cart[XX] - x_cen; 
  yy = x_cart[YY] - y_cen; 
  zz = x_cart[ZZ] - z_cen; 
 
  riso = sqrt(xx*xx + yy*yy + zz*zz);
  r = 0.25 * ( 2.*riso + M + a ) * ( 2.*riso + M - a ) / riso ; 

  x_spher[TT] = riso;
  x_spher[RR] = r; 
  x_spher[TH] = acos(zz / riso); 
  t3 = atan2( yy , xx ); 
  if( t3 < 0. ) {  t3 += 2.*M_PI; } 
  x_spher[PH] = t3;

  return;
}

/****************************************************************************
  setutcon():
 ------------
     -- find the contravariant time-component of a time-like vector 
        pointing forward in time;
****************************************************************************/
static void setutcon(CCTK_REAL *vcon, CCTK_REAL gcov[][NDIM])
{
  CCTK_REAL d,b,c;
  
  d=gcov[TT][TT];

  b = gcov[TT][1]*vcon[1] + gcov[TT][2]*vcon[2] + gcov[TT][3]*vcon[3] ;

  c = gcov[1][1] * vcon[1] * vcon[1] 
    + gcov[2][2] * vcon[2] * vcon[2] 
    + gcov[3][3] * vcon[3] * vcon[3] 
    + 2.*(   gcov[1][2] * vcon[1] * vcon[2] 
	   + gcov[1][3] * vcon[1] * vcon[3] 
           + gcov[2][3] * vcon[2] * vcon[3] );

  c += 1. ;  /* vector is timelike */

  vcon[0]=(-b-sqrt(b*b-d*c))/(d);   /* sign for pointing forward in time */

  return;
}

/****************************************************************************
  bl_gcov_func():
  ---------------
     -- Covariant Kerr metric in Boyer-Lindquist coordinates. 
****************************************************************************/
static void bl_gcov_func( CCTK_REAL *x, CCTK_REAL gcov[NDIM][NDIM])
{
  int i,j ;
  CCTK_REAL sth,cth,s2,r2,DD,mu ;
  CCTK_REAL r,th;

  for( i = 0 ; i < NDIM ; i++ ) {
    for( j = 0 ; j < NDIM ; j++ ) {
      gcov[i][j] = 0. ;
    } 
  }

  r  = x[RR];
  th = x[TH];
  sincos( th, &sth , &cth );   

  s2 = sth*sth ;
  r2 = r*r ;
  DD = 1. - 2.*M/r + asq/r2 ;
  mu = 1. + asq*cth*cth/r2 ;
	
  gcov[TT][TT] = -(1. - 2.*M/(r*mu)) ;
  gcov[TT][3] = -2.*M*a*s2/(r*mu) ;
  gcov[3][TT] = gcov[TT][3] ;
  gcov[1][1] = mu/DD ;
  gcov[2][2] = r2*mu ;
  gcov[3][3] = r2*s2*(1. + asq/r2 + 2.*M*asq*s2/(r2*r*mu)) ;

  return;
}


/***************************************************************************

 set_bondi_parameters():
 ---------------------
   -- finds the values of the hydro. quantities  at the sonic point, which
         serves as a reference point for the conservation equations given 
         in Shapiro and Teukolsky equations (G.21,G.22).  

   -- The sonic point values are then used in find_bondi_solution() to determine
       the hydro. quantities at an arbitrary radius;

   -- the "boundary conditions" that uniquely determine the Bondi solution
       are the radius of the sonic point, "rs",  and the mass accretion rate, 
       Mdot;

   -- the Bondi solution here in is the isentropic, spherically symmetric, 
       perfect fluid solution to Einstein's equations.  That is, we only 
       assume an r-dependence, there's a in-going radial velocity only, 
       and the EOS are :  P = (G-1)*rho*eps   and   P = K rho^G  
       where  K = const.  and  G is the adiabatic constant "gam".  

***************************************************************************/
static void set_bondi_parameters( CCTK_REAL M_in, CCTK_REAL Mdot_in,  CCTK_REAL rs_in, CCTK_REAL gam )
{

  CCTK_REAL  gtemp; 
  static int first_call = 1;

  /* Set the solution-determining parameters:    */
  M    = M_in;
  Mdot = Mdot_in;
  Msq  = M*M;
  rs   = rs_in;
  gamma_eos = gam;

  
  /* Calculate the hydro. quantities: */
  cs_sq  =  M / ( 2.*rs - 3.*M ) ; 

  if( cs_sq > (gam - 1.) ) { 
    cs_sq = gam - 1.;
    rs = 0.5 * M * ( 3. + 1./cs_sq ) ;
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"set_bondi_parameters(): bad value of rs, need to increase it !! \n"); 
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"set_bondi_parameters(): Need to change rs !! \n"); 
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"set_bondi_parameters(): rs[old] = %28.15e    rs[new] = %28.16e  \n\n",rs_in,rs); 
  }

  cs     =  sqrt(cs_sq);
  vs_sq  =  M / ( 2. * rs ) ; 
  vs     =  sqrt(vs_sq); 
  rhos   =  Mdot / ( 4. * M_PI * vs * rs * rs ) ; 
  gtemp  =  gam - 1.;
  hs     =  1. / ( 1. - cs_sq / (gam - 1.) );
  K      = hs * cs_sq * pow( rhos, (-gtemp) ) / gam ; 
  Qdot   = hs * hs * ( 1. - 3. * vs_sq ) ;
  gamma_eos = gam;

  if( first_call ) {
    CCTK_VInfo(CCTK_THORNSTRING,"#######################################################");         
    CCTK_VInfo(CCTK_THORNSTRING,"Bondi Solution Parameters1: ");					  
    CCTK_VInfo(CCTK_THORNSTRING,"-------------------------");					  
    CCTK_VInfo(CCTK_THORNSTRING,"M  = %28.20e     Mdot = %28.20e     rs   = %28.20e",M,Mdot,rs);	  
    CCTK_VInfo(CCTK_THORNSTRING,"vs = %28.20e     cs   = %28.20e     rhos = %28.20e",vs,cs,rhos);	  
    CCTK_VInfo(CCTK_THORNSTRING,"hs = %28.20e     K    = %28.20e     Qdot = %28.20e",hs,K,Qdot);	  
    CCTK_VInfo(CCTK_THORNSTRING,"gam= %28.20e     r_sol= %28.20e",gamma_eos, r_sol);		  
    CCTK_VInfo(CCTK_THORNSTRING,"#######################################################");	  
    first_call = 0;
  }

  return;
  
}

/**********************************************************************/
/************************************************************

  gnr_bondi():
  -----------
    -- should be just like the routine general_newton_raphson() in utoprim*.c 
       except the "physicality" condition is different;

    -- performs Newton-Rapshon method on an arbitrary system
        though tailored to calculate  rho  second Bondi Conservation eq. 
        by ensuring that  rho > 0    (look near "METHOD specific:")

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int gnr_bondi( CCTK_REAL x[], int n, 
	       void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
			      CCTK_REAL [][NEWT_DIM_B], CCTK_REAL *, 
			      CCTK_REAL *, int) )
{
  CCTK_REAL f, df, dx[NEWT_DIM_B], resid[NEWT_DIM_B], 
    jac[NEWT_DIM_B][NEWT_DIM_B];
  CCTK_REAL errx;
  int    n_iter, id, i_extra, doing_extra;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    errx = 0.;

    /* don't use line search : */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at relative error in indep. variable: */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = (x[0] == 0.) ? (SMALL_BONDI)  :  fabs(x[0]);
    

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL_B) && (doing_extra == 0) && (EXTRA_NEWT_ITER_B > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL_B)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER_B) || (n_iter >= (MAX_NEWT_ITER_B-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (isfinite(f)==0) || (isfinite(df)==0)  ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL_B){
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"newt: errx = %28.20e \n", errx); 
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL_B) && (fabs(errx) > NEWT_TOL_B) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL_B ){
    return(0);
  }

  return(0);

}

/******************************************************************************/
/******************************************************************************

 bondi_resid():
 --------------
     -- routine to calculate the residual and jacobian used by 
           the Newton-Raphson routine general_newton_raphson(), which is 
           used to find X2 from theta;

***********************************************************************************/
static void bondi_resid(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], CCTK_REAL jac[][NEWT_DIM_B], 
               CCTK_REAL *f, CCTK_REAL *df, int n )
{
  CCTK_REAL v, vp, h, hp, term;

  hp  =  K * gamma_eos * pow( x[0], (gamma_eos - 2.) );   //   dh/drho 
  h   =  1. +  hp * x[0] / ( gamma_eos - 1. ); 
  v   =  Mdot / ( 4. * M_PI * r_sol * r_sol * x[0] );    
  vp  =  -v / x[0];   //  dv/drho
  term = 1. - 2.*M/r_sol + v*v;
  resid[0]  =  -Qdot  +  h * h * term;
  jac[0][0] =  2. * h *( hp*term + h*v*vp );
  dx[0] = -resid[0] / jac[0][0];
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;

}


/***************************************************************************/
/***************************************************************************

 find_bondi_solution():
 ---------------------

   -- essentially just calls gnr_bondi() to find the solution for
      the density at a given radius, given the parameters calculated from 
      set_bondi_parameters();

   -- after the density is found, the density (rho), internal energy densit (u), magnitude 
        of the radial component of the 4-velocity (v) are returned to the calling 
        routine; 

   -- requires r = radius at which we want solution

   -- note that v is a magnitude, so the user will have to set u^r = -v ;

   -- if there is an error in finding the solution, it returns the error 
       status from the root-finding routine.  See documentation of the 
       gnr_bondi() for further details;

***************************************************************************/
static int find_bondi_solution( CCTK_REAL r, CCTK_REAL *rho, CCTK_REAL *u, CCTK_REAL *v )
{

  int  retval=0;
  const int  ntries = 10000;
  int  itry;

  CCTK_REAL rhotmp;
  CCTK_REAL dr,ur;

  
  /************************************************************************/
  /* Find the initial guess for the newton iterations:                    */
  /* Take the sonic point values if we have no better guess (when rho<0)  */
  /************************************************************************/

  if( *rho < 0. ) {  
    if( r > 0.9*rs && r < 1.1*rs ) { 
      *rho = rhos;
    }
    else { 
      //  rhotmp = (sqrt(Qdot) - 1.) * (gamma_eos - 1.) / ( gamma_eos * K );
      if(r < rs) {  ur = pow(r,-0.5)     ;  }
      else       {  ur = 0.5*pow(r,-1.5) ;  }
      *rho = Mdot / (4.*M_PI * r * r * ur); 
    }
  } 

  // set global variables needed by residual function:
  r_sol = r ; 


  // Use Newton's method to find rho:
  retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     


  // If we were unsure about the guess and solver fails, then creep from known solution to desired point:
  if( retval ) { 

    dr = (r - rs)/(1.*(ntries-1));
    
    *rho = rhos;   // start with sonic point value and near sonic point

    // go gradually away from sonic point toward location where we want the solution:
    r_sol = rs ;
    for( itry = 1; itry < ntries; itry++ ) { 
      r_sol  += dr;

      retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     

      if( retval ) { 
	CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"find_bondi_solution: Incr. guess failed, decrease dfactor, retval = %d, r = %g, itry = %d \n", retval, r, itry);
	return(10);
      }
    }
    
    // No try where we want the solution:
    r_sol = r ; 
    retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid );  

    if( retval ) { 
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"find_bondi_solution: Final Incr. guess failed, decrease dfactor??, retval = %d \n", retval);
      return(11);
    }

  }

  rhotmp = *rho;

  // Calculate other quantities:
  *u = K * pow( rhotmp, gamma_eos ) / (gamma_eos - 1.);

  *v = Mdot / ( 4. * M_PI * r * r * rhotmp );

  if( *u <= 0. ) { retval = -1; } 

  return( retval ) ;

}  

/***********************************************************************************/
/***********************************************************************************
  calc_b_bondi():
  ---------
   -- calculates the contravariant magnetic vector Bcons from the amplitude of the 
        radial component of the contravariant magnetic field vector Boyer-Lindquist coordinates.  

***********************************************************************************/
#if 0
static void  calc_b_bondi( CCTK_REAL B0, CCTK_REAL vtmp, CCTK_REAL x[NDIM], CCTK_REAL x_spher[NDIM],  CCTK_REAL bcon[NDIM] )
{
  int i,j;
  CCTK_REAL  bcon_bl[NDIM] = {0.}; // F^{ab} = (b^a u^b - b^b u^a)
  CCTK_REAL  bcon_ks[NDIM];
  CCTK_REAL  gcov[NDIM][NDIM];
  CCTK_REAL  r_iso;
  CCTK_REAL dxc_dxs[NDIM][NDIM];
  CCTK_REAL f_bl = 1. - 2*M/x[RR];
  CCTK_REAL sqrt_detg_bl = SQR(x[RR])*sin(x[TH])/sqrt(f_bl);
  CCTK_REAL lapse_bl = sqrt(f_bl);

  CCTK_REAL  ucon_bl[NDIM] = {0.};
  CCTK_REAL  w_lorentz_bl; // Lorentz factor in Boyer-Lindquist coords

  assert(a == 0.);

  ucon_bl[RR] = -vtmp;

  /* Find time component of 4-velocity: */
  bl_gcov_func( x_spher, gcov );
  setutcon( ucon_bl, gcov );
  w_lorentz_bl = lapse_bl*ucon_bl[TT];

  /* Find time component of 4-velocity: */
  bcon_bl[RR] = B0*SQR(M)*w_lorentz_bl/sqrt_detg_bl;                                // -- " --
  bcon_bl[TT] = -gcov[RR][RR]/gcov[TT][TT] * bcon_bl[RR] * ucon_bl[RR]/ucon_bl[TT]; // Equ A16 of Villiers&Hawley (corrected)

  switch( coord_type ) {
  
  case COORD_BOYERLINDQUIST :
    dxc_dxs_bl_calc(x, x_spher, dxc_dxs );
    DLOOP1 { bcon[i] = 0.; } 
    DLOOP2 { bcon[i] += dxc_dxs[i][j] * bcon_bl[j] ; }
    break;


  case COORD_KERRSCHILD    :
    bl_to_ks_con(x_spher,bcon_bl,bcon_ks); 
    dxc_dxs_ks_calc(x, x_spher, dxc_dxs );
    DLOOP1 { bcon[i] = 0.; } 
    DLOOP2 { bcon[i] += dxc_dxs[i][j] * bcon_ks[j] ; }
    break;


  case COORD_ISOTROPIC     :
    r_iso = x_spher[TT]; 
    bcon_bl[RR] /=  1. + 0.25 * (asq-Msq) / (r_iso*r_iso)  ;    /* BL to Isotropic coordinate transformation */
    /* I believe we can use BL's cartesian transformation, while using BL's radius : */
    dxc_dxs_iso_calc(x, x_spher, dxc_dxs );
    DLOOP1 { bcon[i] = 0.; } 
    DLOOP2 { bcon[i] += dxc_dxs[i][j] * bcon_bl[j] ; }
    break; 

  }
  
  return;
}
#endif
static void  calc_b_bondi( CCTK_REAL B0, CCTK_REAL vtmp, CCTK_REAL x[NDIM], CCTK_REAL x_spher[NDIM],  CCTK_REAL ucon[NDIM], CCTK_REAL bcon[NDIM])
{
  int i,j;
  CCTK_REAL  bcon_spher[NDIM] = {0.}; // F^{ab} = (b^a u^b - b^b u^a)
  CCTK_REAL  ucon_spher[NDIM] = {0.};
  CCTK_REAL dxc_dxs[NDIM][NDIM];

  assert(a == 0.);
  assert(coord_type == COORD_KERRSCHILD);

  // Schwarzschild and spherical Kerr-Schild (Anil's Eddington coords) ^r components are identical
  ucon_spher[RR] = -vtmp;
  ucon_spher[TT] = (x_spher[RR] + (x_spher[RR]+2*M) * SQR(ucon_spher[RR])) / ( sqrt( (x_spher[RR]-2*M+x_spher[RR]*SQR(ucon_spher[RR]))*x_spher[RR] ) - 2*M*ucon_spher[RR] );

  // Schwarzschild and Kerr-Schild ^t components differ since t_KS = t_SW + 2*M*ln(r/(2*M)-1)
  // NB: we remove the sin(theta) factor from detg_BL  since it factors out of the divergence operator 
  bcon_spher[RR] = B0*SQR(M)/SQR(x_spher[RR]) *  // Equ A16 of Villiers&Hawley plus explict compuation of detg_BL and W_BL
                   sqrt(1-2*M/x_spher[RR] + SQR(ucon_spher[RR]));
  bcon_spher[TT] = B0*SQR(M)/SQR(x_spher[RR]) * 
                   (4*SQR(M) - SQR(ucon_spher[RR])*(SQR(x_spher[RR]+2*M*x_spher[RR]))) / 
                   (2*M*sqrt(SQR(x_spher[RR])-2*M*x_spher[RR]+SQR(x_spher[RR]*ucon_spher[RR])) - 
                    ucon_spher[RR]*SQR(x_spher[RR]));

  dxc_dxs_ks_calc(x, x_spher, dxc_dxs );

  DLOOP1 { bcon[i] = 0.; } 
  DLOOP2 { bcon[i] += dxc_dxs[i][j] * bcon_spher[j] ; }

  DLOOP1 { ucon[i] = 0.; } 
  DLOOP2 { ucon[i] += dxc_dxs[i][j] * ucon_spher[j] ; }

  return;
}

/***********************************************************************************/
/***********************************************************************************
  calc_vel_bondi():
  ---------
   -- calculates the 4-velocity from the amplitude of the 
        radial component of the 4-velocity in Boyer-Lindquist coordinates.  

***********************************************************************************/
#if 0
static void  calc_vel_bondi( CCTK_REAL vtmp, CCTK_REAL x[NDIM], CCTK_REAL x_spher[NDIM],  CCTK_REAL ucon[NDIM] )
{
  int i,j;
  CCTK_REAL  ucon_bl[NDIM] = {0.};
  CCTK_REAL  ucon_ks[NDIM];
  CCTK_REAL  gcov[NDIM][NDIM];
  CCTK_REAL  r_iso;
  CCTK_REAL dxc_dxs[NDIM][NDIM];

  ucon_bl[RR] = -vtmp;

  /* Find time component of 4-velocity: */
  bl_gcov_func( x_spher, gcov );      
  setutcon( ucon_bl, gcov );


  switch( coord_type ) {
  
  case COORD_BOYERLINDQUIST :
    dxc_dxs_bl_calc(x, x_spher, dxc_dxs );
    DLOOP1 { ucon[i] = 0.; } 
    DLOOP2 { ucon[i] += dxc_dxs[i][j] * ucon_bl[j] ; }
    break;


  case COORD_KERRSCHILD    :
    bl_to_ks_con(x_spher,ucon_bl,ucon_ks); 
    dxc_dxs_ks_calc(x, x_spher, dxc_dxs );
    DLOOP1 { ucon[i] = 0.; } 
    DLOOP2 { ucon[i] += dxc_dxs[i][j] * ucon_ks[j] ; }
    break;


  case COORD_ISOTROPIC     :
    r_iso = x_spher[TT]; 
    ucon_bl[RR] /=  1. + 0.25 * (asq-Msq) / (r_iso*r_iso)  ;    /* BL to Isotropic coordinate transformation */
    /* I believe we can use BL's cartesian transformation, while using BL's radius : */
    dxc_dxs_iso_calc(x, x_spher, dxc_dxs );
    DLOOP1 { ucon[i] = 0.; } 
    DLOOP2 { ucon[i] += dxc_dxs[i][j] * ucon_bl[j] ; }
    break; 

  default:
    assert(0 && "Internal error");
    return; /* NOTREACHED */

  }
  
  return;
}
#endif

/***********************************************************************************/
/***********************************************************************************
  GRHydro_BondiM(): 
  ---------
   -- driver routine for the Bondi solution;

***********************************************************************************/
static void GRHydro_BondiM_Internal(CCTK_ARGUMENTS, 
  CCTK_REAL range_min, CCTK_REAL range_max, 
  const int range_imin[], const int range_imax[]);

void GRHydro_BondiM(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  const int null[3] = {0,0,0};
  GRHydro_BondiM_Internal(CCTK_PASS_CTOC, 1e100, -1e100, cctk_lsh, null);
}

void GRHydro_BondiM_Range(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  const int null[3] = {0,0,0};
  GRHydro_BondiM_Internal(CCTK_PASS_CTOC, bondi_freeze_inner_radius,
                          bondi_freeze_outer_radius,
                          cctk_lsh, null);
}

void GRHydro_BondiM_Boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  const int imax[3] = {cctk_lsh[0]-cctk_nghostzones[0],
                       cctk_lsh[1]-cctk_nghostzones[1],
                       cctk_lsh[2]-cctk_nghostzones[2]};
  GRHydro_BondiM_Internal(CCTK_PASS_CTOC, 1e100, -1e100, 
                          cctk_nghostzones, imax);
}

static void GRHydro_BondiM_Internal(CCTK_ARGUMENTS, 
  CCTK_REAL range_min, CCTK_REAL range_max,
  const int range_imin[], const int range_imax[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL det, inv_alp;
  CCTK_REAL rhotmp, utmp, vtmp, rspher, xpos[NDIM], x_spher[NDIM], ucon[NDIM], bcon[NDIM];
  int retval;
  CCTK_REAL  *r_bondi, *logr_bondi, *rho_bondi, *u_bondi, *v_bondi;
  CCTK_REAL  dlogr,logrmin;
  CCTK_REAL  rmin_bondi,rmax_bondi;
  int bvec_method;
  enum {BVEC_DIRECT, BVEC_TRANSFORM};

  int nbondi; 
  
  const int size = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
  int i, imin, j;

#    define velx(i_)  vel[i_ + 0*size]
#    define vely(i_)  vel[i_ + 1*size]
#    define velz(i_)  vel[i_ + 2*size]
#    define   sx(i_) scon[i_ + 0*size]
#    define   sy(i_) scon[i_ + 1*size]
#    define   sz(i_) scon[i_ + 2*size]
#    define Bconsx(i_) Bcons[i_ + 0*size]
#    define Bconsy(i_) Bcons[i_ + 1*size]
#    define Bconsz(i_) Bcons[i_ + 2*size]
#    define Bvecx(i_) Bvec[i_ + 0*size]
#    define Bvecy(i_) Bvec[i_ + 1*size]
#    define Bvecz(i_) Bvec[i_ + 2*size]

  if( CCTK_EQUALS(bondi_coordinates,"Boyer-Lindquist") ) { 
    coord_type = COORD_BOYERLINDQUIST ;
  } else if( CCTK_EQUALS(bondi_coordinates,"Kerr-Schild") ) { 
    coord_type = COORD_KERRSCHILD  ;
  } else if( CCTK_EQUALS(bondi_coordinates,"Isotropic") ) { 
    coord_type = COORD_ISOTROPIC ;
  } else {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Unknown coordinate type '%s'", bondi_coordinates);
  }

  if( CCTK_EQUALS(bondi_Bvec_method, "direct")) {
    bvec_method = BVEC_DIRECT;
  } else if( CCTK_EQUALS(bondi_Bvec_method, "transform")) {
    bvec_method = BVEC_TRANSFORM;
  } else {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
               "Unknown bvec setup method '%s'", bondi_Bvec_method);
  }



  /* xyz location of the black hole : */
  i = 0; 
  x_cen = bh_bondi_pos_x[i] ; 
  y_cen = bh_bondi_pos_y[i] ; 
  z_cen = bh_bondi_pos_z[i] ; 

  M      = bondi_central_mass[i];
  a      = M*bondi_central_spin[i];
  asq    = a * a; 
  rmin_bondi = M * bondi_rmin[i];
  rmax_bondi = M * bondi_rmax[i];

  nbondi = n_bondi_pts[i];

  set_bondi_parameters( M, mdot_sonicpt_bondi, r_sonicpt_bondi, gl_gamma);


#if(LTRACE)
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\n##################################################\n");                       
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"  Geometry PARAMETERS \n------------------------------------\n");		  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  M           =  %28.18e \n",M       ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  a           =  %28.18e \n",a       ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  x_cen       =  %28.18e \n",x_cen   ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  y_cen       =  %28.18e \n",y_cen   ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  z_cen       =  %28.18e \n",z_cen   ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  rmin_bondi  =  %28.18e \n",rmin_bondi   ); 				  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  rmax_bondi  =  %28.18e \n",rmax_bondi   ); 				  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  nbondi      =  %d      \n",nbondi   ); 					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\n------------------------------------\n");					  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"   SOLUTION PARAMETERS \n------------------------------------\n");		  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  bondi_coordinates          =  %s      \n",bondi_coordinates ); 		  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  mdot_sonicpt_bondi         =  %28.18e \n",mdot_sonicpt_bondi); 		  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  r_sonicpt_bondi            =  %28.18e \n",r_sonicpt_bondi   ); 		  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  gl_gamma    =  %28.18e \n",gl_gamma     ); 	  
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\n------------------------------------\n");					  
#endif

  if( r_sonicpt_bondi < rmin_bondi ) { 
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"sonic point lies below solution domain!!");
  }

  if( r_sonicpt_bondi > rmax_bondi ) { 
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"sonic point lies beyond solution domain!!");
  }


  /*********************************************************************************
     ARRAY ALLOCATIONS :
  *********************************************************************************/
  /* global solution : */
  if( (r_bondi = (CCTK_REAL *) calloc(nbondi, sizeof(CCTK_REAL))) == NULL ) { 
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc r_bondi \n");
    return;
  }
  if( (logr_bondi = (CCTK_REAL *) calloc(nbondi, sizeof(CCTK_REAL))) == NULL ) { 
    free(r_bondi);    
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc logr_bondi \n");
    return;
  }
  if( (rho_bondi = (CCTK_REAL *) calloc(nbondi, sizeof(CCTK_REAL))) == NULL ) { 
    free(r_bondi);    free(logr_bondi);
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc rho_bondi \n");
    return;
  }
  if( (u_bondi = (CCTK_REAL *) calloc(nbondi, sizeof(CCTK_REAL))) == NULL ) { 
    free(r_bondi);    free(logr_bondi);    free(rho_bondi);
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc u_bondi \n");
    return;
  }
  if( (v_bondi = (CCTK_REAL *) calloc(nbondi, sizeof(CCTK_REAL))) == NULL ) { 
    free(r_bondi);    free(logr_bondi);    free(rho_bondi);    free(u_bondi);
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc v_bondi \n");
    return;
  }
  

  /*********************************************************************************
     SOLUTION DOMAIN :
  *********************************************************************************/
  logrmin = log10(rmin_bondi);
  dlogr = (log10(rmax_bondi) - logrmin)/(1.*(nbondi-1));

  for(i=0; i < nbondi; i++) {   
    logr_bondi[i] = logrmin + dlogr*i;   
  }

  for(i=0; i < nbondi; i++) {   
    r_bondi[i] = pow(10.,logr_bondi[i]);
  }

  rhotmp = 1.e200; 
  imin = 0; 

  /* find the position in the array where the sonic point lies */
  for(i=0; i < nbondi; i++) {
    utmp = fabs(r_bondi[i] - r_sonicpt_bondi); 
    if( utmp < rhotmp ) { 
      rhotmp = utmp; 
      imin = i ; 
    }
  }

  /*********************************************************************************
     DETERMINE BONDI SOLUTION :
  *********************************************************************************/

  /* start at the sonic point (where we know the solution) and spread out from there using the 
     adjacent point as the guess for the next furthest point : */
  rhotmp = -1.;  // start with guess

  for(i=imin; i < nbondi; i++) {
    rspher = r_bondi[i]; 
    retval = find_bondi_solution( rspher, &rhotmp, &utmp, &vtmp );
    if( retval ) { 
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"Problem1 with find_bondi_solution() at  i = %d,  r = %28.16e \n",i,rspher);
      rhotmp = -1. ;      vtmp = 0.;       /* trigger the floor and set to staticity  */
    }
    if(rhotmp < initial_rho_abs_min)  {
      rhotmp = initial_rho_abs_min;
      utmp = K * pow( rhotmp, gl_gamma ) / (gl_gamma - 1.);
    }
    rho_bondi[i] = rhotmp;    u_bondi[i]   = utmp;     v_bondi[i]   = vtmp; 
  }

  rhotmp = -1.;  // start with guess

  for(i=imin-1; i >= 0; i--) {
    rspher = r_bondi[i]; 
    retval = find_bondi_solution( rspher, &rhotmp, &utmp, &vtmp );
    if( retval ) { 
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"Problem2 with find_bondi_solution() at  i = %d,  r = %28.16e \n",i,rspher);
      rhotmp = -1. ;      vtmp = 0.;       /* trigger the floor and set to staticity  */
    }
    if(rhotmp < initial_rho_abs_min)  {
      rhotmp = initial_rho_abs_min;
      utmp = K * pow( rhotmp, gl_gamma ) / (gl_gamma - 1.);
    }
    rho_bondi[i] = rhotmp;    u_bondi[i]   = utmp;     v_bondi[i]   = vtmp; 
  }

#if(LTRACE) 
  for(i=0; i < nbondi; i++) {  
    CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"radial-bondisoln %12d %26.16e %26.16e %26.16e %26.16e \n",i,r_bondi[i],rho_bondi[i],u_bondi[i],v_bondi[i]);  
  }
#endif

  /*********************************************************************************
     LOAD GRID FUNCTIONS WITH THE SOLUTION
  *********************************************************************************/
  xpos[TT] = 0.;
  
  for(int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
  for(int jj = 0 ; jj < cctk_lsh[1] ; jj++) {
  for(int ii = 0 ; ii < cctk_lsh[0] ; ii++) {

    i = CCTK_GFINDEX3D(cctkGH, ii,jj,kk);

    xpos[XX] = x[i] ;
    xpos[YY] = y[i] ;
    xpos[ZZ] = z[i] ;

    if (bondi_radial_offset > 0.) {
      double rspher_orig = sqrt(SQR(xpos[XX])+SQR(xpos[YY])+SQR(xpos[ZZ]));
      double rspher_new = rspher_orig + bondi_radial_offset;
      if(rspher_orig < SMALL_BONDI) { 
        xpos[XX] = rspher_new;
        xpos[YY] = 0.;
        xpos[ZZ] = 0.;
      } else {
        for(int n = XX ; n <= ZZ ; n++) {
          xpos[n] *= rspher_new/rspher_orig;
        }
      }
    }


    switch( coord_type ) {
  
    case COORD_BOYERLINDQUIST :
      bl_cart_to_bl_spher_pos( xpos, x_spher);
      break;

    case COORD_KERRSCHILD    :
      ks_cart_to_ks_spher_pos( xpos, x_spher);
      break;

    case COORD_ISOTROPIC     :
      iso_cart_to_bl_spher_pos(xpos, x_spher); 
      break; 

    default:
      assert(0 && "Internal error");
      return; /* NOTREACHED */

    }

    rspher = x_spher[RR]; 

    // the conditions here are a bit strange. The logic goes like this:
    // range_min,range_max (and range_imin, range_imax) are radii (and indices)
    // bounding the volume where we want to evolve the fields. So during
    // evolution we only reset to ID data if we are NOT inside this volume. For
    // ID we pass special values for min/max such that these conditions never
    // trigger. 
    // TODO: do this more nicely. The fact that a comment is needed should tell
    // me that this is getting too complex.
    if(rspher > range_min && rspher < range_max)
      continue;
    if(ii >= range_imin[0] && ii < range_imax[0] &&
       jj >= range_imin[1] && jj < range_imax[1] &&
       kk >= range_imin[2] && kk < range_imax[2]) {
      continue;
    }

    /* Find nearest point in the Bondi solution : */
    j = (int)  ( 0.5 +  (log10(rspher) - logrmin) / dlogr ) ; 

    if(  j < 0  ) { 
      CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"Grid point resides outside solution domain:  %26.16e %26.16e %26.16e %26.16e %26.16e\n", rspher, rmin_bondi, xpos[XX],xpos[YY],xpos[ZZ]); 
      j = 0; 
    }
    else if( j > (nbondi - 1) ) {
      CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"GRHydro_Bondi::  Grid point resides outside solution domain:  %26.16e %26.16e %26.16e %26.16e %26.16e\n", rspher, rmin_bondi, xpos[XX],xpos[YY],xpos[ZZ]); 
      j = nbondi - 1; 
    }

    rhotmp = rho_bondi[j];
    retval = find_bondi_solution( rspher, &rhotmp, &utmp, &vtmp );
    if( retval ) { 
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"Problem3 with find_bondi_solution() at  i = %d,  r = %28.16e \n",i,rspher);
      rhotmp = -1. ;      vtmp = 0.;       /* trigger the floor and set to staticity  */
    }

    if(rhotmp < initial_rho_abs_min)  {
      rhotmp = initial_rho_abs_min;
      utmp = K * pow( rhotmp, gl_gamma ) / (gl_gamma - 1.);
    }

    rho[i] = rhotmp;
    eps[i] = utmp/rhotmp;
    calc_b_bondi(bondi_bmag, vtmp, xpos, x_spher, ucon, bcon);

    // verify normaliztion of ucon
    {
      CCTK_REAL gtt = -SQR(alp[i]) +
                      gxx[i]*SQR(betax[i]) + gyy[i]*SQR(betay[i]) +
                      gzz[i]*SQR(betaz[i]) + 
                      2*gxy[i]*betax[i]*betay[i] + 2*gxz[i]*betax[i]*betaz[i] +
                      2*gyz[i]*betay[i]*betaz[i];
      CCTK_REAL gtx = gxx[i]*betax[i]+gxy[i]*betay[i]+gxz[i]*betaz[i];
      CCTK_REAL gty = gxy[i]*betax[i]+gyy[i]*betay[i]+gyz[i]*betaz[i];
      CCTK_REAL gtz = gxz[i]*betax[i]+gyz[i]*betay[i]+gzz[i]*betaz[i];
      CCTK_REAL umag = gtt*SQR(ucon[TT]) +
                       gxx[i]*SQR(ucon[XX]) + gyy[i]*SQR(ucon[YY]) + 
                       gzz[i]*SQR(ucon[ZZ]) +
                       2*(gtx*ucon[XX]*ucon[TT] + gty*ucon[YY]*ucon[TT] +
                          gtz*ucon[ZZ]*ucon[TT]) +
                       2*(gxy[i]*ucon[XX]*ucon[YY] + gxz[i]*ucon[XX]*ucon[ZZ] +
                          gyz[i]*ucon[YY]*ucon[ZZ]);
      CCTK_REAL abssum = fabs(gtt*SQR(ucon[TT])) + 
                       fabs(gxx[i]*SQR(ucon[XX])) + fabs(gyy[i]*SQR(ucon[YY])) + 
                       fabs(gzz[i]*SQR(ucon[ZZ])) +
                       2*(fabs(gtx*ucon[XX]*ucon[TT]) + fabs(gty*ucon[YY]*ucon[TT]) +
                          fabs(gtz*ucon[ZZ]*ucon[TT])) +
                       2*(fabs(gxy[i]*ucon[XX]*ucon[YY]) + fabs(gxz[i]*ucon[XX]*ucon[ZZ]) +
                          fabs(gyz[i]*ucon[YY]*ucon[ZZ]));
      if(fabs(umag-(-1.0)) > 1e-13*abssum) {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Suspicious four velocity (%.16e,%.16e,%.16e,%.16e) with "
                   "normalization %.16e at  i = %d x = (%g,%g,%g), r = %28.16e, "
                   "3-metric = (%.16e,%.16e,%.16e,%.16e,%.16e,%.16e) "
                   "gtt = %.16e, beta_i = (%.16e,%.16e,%.16e)",
                   ucon[TT],ucon[XX],ucon[YY],ucon[ZZ],umag,i,
                   x[i],y[i],z[i],rspher,
                   gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],gtt,gtx,gty,gtz);
      }
    }

    inv_alp = 1./alp[i];

    velx(i) = (ucon[XX]/ucon[TT] + betax[i]) * inv_alp;
    vely(i) = (ucon[YY]/ucon[TT] + betay[i]) * inv_alp;
    velz(i) = (ucon[ZZ]/ucon[TT] + betaz[i]) * inv_alp;

    SpatialDet(gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],&det);

    if(bvec_method == BVEC_TRANSFORM) {
      const CCTK_REAL w = ucon[TT] * alp[i];
      Bvecx(i) = w*bcon[XX] - alp[i]*bcon[TT]*ucon[XX];
      Bvecy(i) = w*bcon[YY] - alp[i]*bcon[TT]*ucon[YY];
      Bvecz(i) = w*bcon[ZZ] - alp[i]*bcon[TT]*ucon[ZZ];
    } else {
      Bvecx(i) = bondi_bmag*SQR(M)*x[i]/sqrt(det)/CUBE(r[i]);
      Bvecy(i) = bondi_bmag*SQR(M)*y[i]/sqrt(det)/CUBE(r[i]);
      Bvecz(i) = bondi_bmag*SQR(M)*z[i]/sqrt(det)/CUBE(r[i]);
    }

    // damp everything down to atmosphere inside of 1M
    if(rspher < M) {
      const CCTK_REAL smooth = 0.5*(1+tanh(tan(M_PI*(rspher/M-0.5))));
      rho[i] *= smooth;  
      velx(i) *= smooth;
      vely(i) *= smooth;
      velz(i) *= smooth;
      eps[i] *= smooth;
      Bvecx(i) *= smooth;
      Bvecy(i) *= smooth;
      Bvecz(i) *= smooth;
    }

    Prim2ConGenM(*GRHydro_eos_handle,gxx[i],gxy[i],
                 gxz[i],gyy[i],gyz[i],gzz[i],
                 det, &dens[i],&sx(i),&sy(i),&sz(i),
                 &tau[i],
                 &Bconsx(i),&Bconsy(i),&Bconsz(i),
                 rho[i],
                 velx(i),vely(i),velz(i),
                 eps[i],&press[i],
                 Bvecx(i),Bvecy(i),Bvecz(i),
                 &w_lorentz[i]);

    GRHydro_C2P_failed[i] = 0.;
    if(psidc)
      psidc[i] = 0.;
                 
  } // ii
  } // jj
  } // kk

  free(r_bondi);      free(logr_bondi);      free(rho_bondi);     free(u_bondi);   free(v_bondi);


#    undef velx
#    undef vely
#    undef velz
#    undef sx
#    undef sy
#    undef sz
#    undef Bconsx
#    undef Bconsy
#    undef Bconsz
#    undef Bvecx
#    undef Bvecy
#    undef Bvecz


  return;
}




#undef LTRACE 
#undef NEWT_DIM_B        
#undef MAX_NEWT_ITER_B   
#undef NEWT_TOL_B        
#undef MIN_NEWT_TOL_B    
#undef EXTRA_NEWT_ITER_B 
#undef SMALL_BONDI 
#undef NDIM 
#undef TT
#undef RR
#undef TH
#undef PH
#undef XX
#undef YY
#undef ZZ
#undef COORD_BOYERLINDQUIST
#undef COORD_KERRSCHILD    
#undef COORD_ISOTROPIC     
#undef DLOOP1
#undef DLOOP2
