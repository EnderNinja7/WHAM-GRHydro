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


#if !defined(__INTEL_COMPILER)
# define LOCAL_SINCOS
# define sincos( theta_ , sth_ , cth_ )  {  *(sth_) = sin((theta_)) ;  *(cth_) = cos((theta_)) ;  }
#endif


//Newton-Raphson parameters:
#define NEWT_DIM_B        (1      )
#define MAX_NEWT_ITER_B   (30     )    /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL_B        (1.0e-15)    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL_B    (1.0e-10)    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER_B (2      )
#define SMALL_BONDI (1.e-20)


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
       and the EOS are :  P = (G-1)*rho   and   P = K rho^G  
       where  K = const.  and  G is the adiabatic constant "gam".  

***************************************************************************/
static void set_bondi_parameters( CCTK_REAL M_in, CCTK_REAL Mdot_in,  CCTK_REAL rs_in, CCTK_REAL gam )
{

  CCTK_REAL  gtemp; 

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

#if( LTRACE ) 
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\n#######################################################\n");         
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"Bondi Solution Parameters1: \n");
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"------------------------- \n\n");
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"M  = %28.20e     Mdot = %28.20e     rs   = %28.20e  \n",M,Mdot,rs);
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"vs = %28.20e     cs   = %28.20e     rhos = %28.20e  \n",vs,cs,rhos);
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"hs = %28.20e     K    = %28.20e     Qdot = %28.20e   \n",hs,K,Qdot);
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"gam= %28.20e     r_sol= %28.20e       \n",gamma_eos, r_sol);
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"rs   = : %28.20e \n", rs) ;
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"urs  = : %28.20e \n", vs) ;
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"rhos = : %28.20e \n", rhos) ;
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"K    = : %28.20e \n", K) ;
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"#######################################################\n\n");
#endif

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

  CCTK_REAL rhotmp, rho_guess;
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
      //  rho_guess = pow( rhotmp , (1./(gamma_eos - 1.)) );
      if(r < rs) {  ur = pow(r,-0.5)     ;  }
      else       {  ur = 0.5*pow(r,-1.5) ;  }
      *rho = Mdot / (4.*M_PI * r * r * ur); 
    }
  }

  // safe guess value for multiple tries
  rho_guess = *rho;
  

  // set global variables needed by residual function:
  r_sol = r ; 


  // Use Newton's method to find rho:
  retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     


  // first try guess if failure 
  if( retval ) { 
    *rho = rho_guess;
    retval = gnr_bondi( rho, NEWT_DIM_B, bondi_resid);     
  }

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
  calc_vel_bondi():
  ---------
   -- calculates the 4-velocity from the amplitude of the 
        radial component of the 4-velocity in Boyer-Lindquist coordinates.  

***********************************************************************************/
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

/***********************************************************************************/
/***********************************************************************************
  GRHydro_Bondi(): 
  ---------
   -- driver routine for the Bondi solution;

***********************************************************************************/
void GRHydro_Bondi(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL det;
  CCTK_REAL rhotmp, utmp, vtmp, rspher, xpos[NDIM], x_spher[NDIM], ucon[NDIM];
  int retval;
  CCTK_REAL  *r_bondi, *logr_bondi, *rho_bondi, *u_bondi, *v_bondi;
  CCTK_INT   *bad_bondi;
  CCTK_REAL  dlogr,logrmin;
  CCTK_REAL  rmin_bondi,rmax_bondi;
//  CCTK_INT GRHydro_reflevel=0;

  CCTK_INT nbondi; 
  CCTK_INT ileft, iright, nans_exist;
  CCTK_REAL sigma;
  
  
  const int size = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  int i, imin, j;

#    define velx(i_)  vel[i_ + 0*size]
#    define vely(i_)  vel[i_ + 1*size]
#    define velz(i_)  vel[i_ + 2*size]
#    define   sx(i_) scon[i_ + 0*size]
#    define   sy(i_) scon[i_ + 1*size]
#    define   sz(i_) scon[i_ + 2*size]

  if( CCTK_EQUALS(bondi_coordinates,"Boyer-Lindquist") ) { 
    coord_type = COORD_BOYERLINDQUIST ;
  }
  if( CCTK_EQUALS(bondi_coordinates,"Kerr-Schild") ) { 
    coord_type = COORD_KERRSCHILD  ;
  }
  if( CCTK_EQUALS(bondi_coordinates,"Isotropic") ) { 
    coord_type = COORD_ISOTROPIC ;
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
  CCTK_VWarn(CCTK_WARN_DEBUG, __LINE__, __FILE__, CCTK_THORNSTRING,"\t  nbondi      =  %d      \n",(int)nbondi  ); 					  
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
  if( (bad_bondi = (CCTK_INT *) calloc(nbondi, sizeof(CCTK_INT))) == NULL ) { 
    free(r_bondi);    free(logr_bondi);    free(rho_bondi);    free(u_bondi); free(v_bondi);
    CCTK_WARN(CCTK_WARN_ABORT,"Cannot alloc bad_bondi \n");
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


  /* Verify Bondi solution, extrapolate/interpolate over Nans */ 
  nans_exist = 0 ; 

  for(i=0; i < nbondi; i++) {   bad_bondi[i] = 0 ; } 
  for(i=0; i < nbondi; i++) {     
    if( (!isfinite(rho_bondi[i])) || (!isfinite(u_bondi[i])) || (!isfinite(v_bondi[i])) || (rho_bondi[i]==0.) || (v_bondi[i]==0.) || (u_bondi[i]==0.) ) {
      bad_bondi[i] = 1 ; 
      nans_exist++; 
    }
  }

  if( nans_exist > (nbondi-3) ) { 
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,"Too many bad points %12d out of  %12d  \n",(int)nans_exist,(int)nbondi); 
  }

  /* We find the 2-point stencil by favoring interpolation over extrapolation, no matter the distance between the points */
  if( nans_exist ) { 
    for(i=0; i < nbondi; i++) { 
      if( bad_bondi[i] ) { 
	j = i-1;
	while( (j > 0) && (bad_bondi[j]) ) { j--; }

	if( (j==0) && bad_bondi[j] ) { 
	  /* No good points to the left, need two points on the right: */

	  while( (j < (nbondi-1)) && (bad_bondi[j]) ) { j++; }
	  if( (j==(nbondi-1)) && bad_bondi[j] ) { 
	    CCTK_WARN(CCTK_WARN_ABORT,"No available points with which to interpolate!\n");
            ileft = 0; /* NOTREACED */
            iright = 0; /* NOTREACHED */
	  }
	  else { 
	    ileft = j ; 

	    /* Continue searching to the right: */
	    j++; 
	    while( (j < (nbondi-1)) && (bad_bondi[j]) ) { j++; }
	    if( (j==(nbondi-1)) && bad_bondi[j] ) { 
	      CCTK_WARN(CCTK_WARN_ABORT,"Need another point to the right but cannot find one! \n");
              iright = 0; /* NOTREACHED */
	    }
	    else { 
	      iright = j ; 
	    }
	  }
	}
	else { 
	  /* Found a left point, need a right: */
	  ileft = j ; 

	  j = i+1;
	  while( (j < (nbondi-1)) && (bad_bondi[j]) ) { j++; }
	  if( (j==(nbondi-1)) && bad_bondi[j] ) { 
	    /* No good points to the right, need to find another point on the left: */
	    iright = ileft; 
	    
	    j = ileft-1; 
	    while( (j > 0) && (bad_bondi[j]) ) { j--; }
	    if( (j==0) && bad_bondi[j] ) { 
	      CCTK_WARN(CCTK_WARN_ABORT,"Need another point to the left but cannot find one! \n");
              ileft = 0; /* NOTREACED */
	    }
	    else {
	      ileft = j;
	    }
	  }
	  else { 
	    iright = j;
	  }
	}
	/* Now interpolate over the bad point with the good points found: */ 
	
	sigma = (r_bondi[i]-r_bondi[ileft])/(r_bondi[iright]-r_bondi[ileft]); 
	rho_bondi[i] = rho_bondi[iright] * sigma   +  rho_bondi[ileft] * (1.-sigma) ;
	  u_bondi[i] =   u_bondi[iright] * sigma   +    u_bondi[ileft] * (1.-sigma) ;
	  v_bondi[i] =   v_bondi[iright] * sigma   +    v_bondi[ileft] * (1.-sigma) ;
      }
    }
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
  
  for(i=0; i < size; i++) {

    xpos[XX] = x[i] ;
    xpos[YY] = y[i] ;
    xpos[ZZ] = z[i] ;


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

    if( rspher < rmin_bondi )  rspher = rmin_bondi;
    

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
    calc_vel_bondi(vtmp, xpos, x_spher, ucon);

    det = 1./alp[i];   /* temp var */

    velx(i) = (ucon[XX]/ucon[TT] + betax[i]) * det;
    vely(i) = (ucon[YY]/ucon[TT] + betay[i]) * det;
    velz(i) = (ucon[ZZ]/ucon[TT] + betaz[i]) * det;


    SpatialDet(gxx[i],gxy[i],gxz[i],gyy[i],gyz[i],gzz[i],&det);

      Prim2ConGen(*GRHydro_eos_handle,gxx[i],gxy[i],
		  gxz[i],gyy[i],gyz[i],gzz[i],
		  det, &dens[i],&sx(i),&sy(i),&sz(i),
		  &tau[i],rho[i],
		  velx(i),vely(i),velz(i),
		  eps[i],&press[i],&w_lorentz[i]);

      if( (!isfinite(dens[i])) || (!isfinite(sx(i))) || (!isfinite(sy(i))) || (sz(i)==0.) || (tau[i]==0.) || (press[i]==0.) || (w_lorentz[i]==0.) ) {
	CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,"bad point at %12d  :  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  %26.16e  \n",i,rspher, rho[i], eps[i], press[i], velx(i), vely(i), velz(i), dens[i], sx(i), sy(i), sz(i) ); 
	dens[i] = rho_bondi[0];
	tau[i]  =   u_bondi[0];
	press[i] = (gamma_eos-1.)*u_bondi[0];
	w_lorentz[i]  = 1.;
	sx(i) = 0.;	sy(i) = 0.;	sz(i) = 0.;
      }

  }

  free(r_bondi);      free(logr_bondi);      free(rho_bondi);     free(u_bondi);   free(v_bondi);


#    undef velx
#    undef vely
#    undef velz
#    undef sx
#    undef sy
#    undef sz


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

#ifdef LOCAL_SINCOS
# undef sincos
# undef LOCAL_SINCOS
#endif
