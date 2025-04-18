#include <cmath>
#include <iostream>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "GRHydro_Reconstruct_drv_impl.hh"
#include "GRHydro_WENOReconstruct.hh"

using namespace std;

template <typename T> static inline T SQR (T const & x) { return x*x; }

/**
   WENO5 reconstruction operator.
   Supports standard WENO5 (with and without adaptive epsilon), and WENO-z.
*/
template <bool do_wenoz, bool do_adaptive_epsilon>
template <int dir>
void GRHydro_WENOReconstruct1d_cxx<do_wenoz,do_adaptive_epsilon>::
apply(const int nx, const CCTK_REAL* const restrict a,
      CCTK_REAL* const restrict aminus, CCTK_REAL* const restrict aplus,
      const cGH* const cctkGH, const int j, const int k)
{
   DECLARE_CCTK_PARAMETERS;
   
#define A(i_) (a[ijk[i_]])
#define Aplus(i_) (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])
   
   for (int i=GRHydro_stencil-1; i < nx-GRHydro_stencil+1; ++i)
   {
      const int ijk[5] = {
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i-2, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i-2, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i-2), 
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i-1, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i-1, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i-1),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i  , j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i  , k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i  ),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i+1, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i+1, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i+1),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i+2, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i+2, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i+2)
                         };
                       
   
         
   
      static_assert (! (do_wenoz && do_adaptive_epsilon), "Adaptive_epsilon not supported for WENO-Z");

      // from https://arxiv.org/abs/0704.2608 (WHAM paper) but generic for WENO
      // these are from equ. (18) and (19)  (not A18) when substituting into
      // each other and sorting by the index 0 <= (i-r+j) <= 4 instead of 0 <=
      // j <= 4 (and 0 <= r <= 2) and the matricies for c_ij are from appendix
      // A2
      //
      // Maple code for the c_{rj} in (18):
      //
      // First is the matrix (*interpolate* from center to left side of cell):
      //
      // C_cl := matrix([[15/8, -5/4, 3/8], [3/8, 3/4, -1/8], [-1/8, 3/4, 3/8]]);
      //
      // Next is the matrix (de-average to center of cell):
      //
      // C_ac := matrix([[23/24, 1/12, -1/24],[-1/24, 13/12, -1/24],[-1/24, 1/12, 23/24]]);
      //
      // we want to de-average and interpolate so:
      //
      // C_al := multiply(C_cl, C_ac);
      //
      // which gives:
      //
      //        [ 11/6  -7/6   1/3 ]
      // C_al = [  1/3   5/6  -1/6 ]
      //        [ -1/6   5/6   1/3 ]
      //
      // Order ends up being inverted (different counting convention).
      //
      static const CCTK_REAL 
                      weno_coeffs[3][5] = { { 2.0/6.0, -7.0/6.0, 11.0/6.0, 0,        0 }, 
                                            { 0,       -1.0/6.0, 5.0/6.0,  2.0/6.0,  0 },
                                            { 0,        0,       2.0/6.0,  5.0/6.0, -1.0/6.0 } };
      if (do_wenoz)
      {
      
         const CCTK_REAL beta1 = 13.0/12.0*SQR(A(0)-2.0*A(1)+A(2)) + 1.0/4.0*SQR(A(0)-4.0*A(1)+3.0*A(2));
         const CCTK_REAL beta2 = 13.0/12.0*SQR(A(1)-2.0*A(2)+A(3)) + 1.0/4.0*SQR(A(1)-A(3));
         const CCTK_REAL beta3 = 13.0/12.0*SQR(A(2)-2.0*A(3)+A(4)) + 1.0/4.0*SQR(3.0*A(2)-4.0*A(3)+A(4));
            
            
         //    Compute weights according to weno-z alorithm
         const CCTK_REAL wbarplus1 = 1.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta1));
         const CCTK_REAL wbarplus2 = 3.0/5.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta2));
         const CCTK_REAL wbarplus3 = 3.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta3));

         const CCTK_REAL wbarminus1 = 3.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta1));
         const CCTK_REAL wbarminus2 = 3.0/5.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta2));
         const CCTK_REAL wbarminus3 = 1.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta3));
         
         const CCTK_REAL iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3);
         
         const CCTK_REAL wplus1 = wbarplus1 * iwbarplussum;
         const CCTK_REAL wplus2 = wbarplus2 * iwbarplussum;
         const CCTK_REAL wplus3 = wbarplus3 * iwbarplussum;
         
         const CCTK_REAL iwbarminussum = 1.0 / (wbarminus1 + wbarminus2 + wbarminus3);
         
         const CCTK_REAL wminus1 = wbarminus1 * iwbarminussum;
         const CCTK_REAL wminus2 = wbarminus2 * iwbarminussum;
         const CCTK_REAL wminus3 = wbarminus3 * iwbarminussum;
         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         Aminus(2) = 0;
         for (int n=0; n < 5; ++n) {
               Aplus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);
               Aminus(2) += (wminus1 * weno_coeffs[2][4-n]
                           + wminus2 * weno_coeffs[1][4-n]
                           + wminus3 * weno_coeffs[0][4-n]) * A(n);
         }
      } else {
         
         // equ. A18 of https://arxiv.org/abs/0704.2608 (WHAM paper) when exapnded out
         // beta (omega_r in the paper) is the same for center to interface,
         // average to center and center to average reconstruction (appendix A3
         // before equ A18)
         static const CCTK_REAL beta_shu[3][6] = { { 4.0/3.0,  -19.0/3.0, 25.0/3.0, 11.0/3.0, -31.0/3.0, 10.0/3.0 },
                                      { 4.0/3.0,  -13.0/3.0, 13.0/3.0, 5.0/3.0,  -13.0/3.0, 4.0/3.0 },
                                      { 10.0/3.0, -31.0/3.0, 25.0/3.0, 11.0/3.0, -19.0/3.0, 4.0/3.0 } };

         // Compute smoothness indicators
         // This is from Tchekhovskoy et al 2007 (WHAM code paper).
         CCTK_REAL beta1  = beta_shu[0][0]*SQR(A(0))
                  + beta_shu[0][1]*A(0)*A(1)
                  + beta_shu[0][2]*SQR(A(1))
                  + beta_shu[0][3]*A(0)*A(2)
                  + beta_shu[0][4]*A(1)*A(2)
                  + beta_shu[0][5]*SQR(A(2));
         
         CCTK_REAL beta2  = beta_shu[1][0]*SQR(A(1))
                  + beta_shu[1][1]*A(1)*A(2)
                  + beta_shu[1][2]*SQR(A(2))
                  + beta_shu[1][3]*A(1)*A(3)
                  + beta_shu[1][4]*A(2)*A(3)
                  + beta_shu[1][5]*SQR(A(3));
         
         CCTK_REAL beta3  = beta_shu[2][0]*SQR(A(2))
                  + beta_shu[2][1]*A(2)*A(3)
                  + beta_shu[2][2]*SQR(A(3))
                  + beta_shu[2][3]*A(2)*A(4)
                  + beta_shu[2][4]*A(3)*A(4)
                  + beta_shu[2][5]*SQR(A(4));
         
         
         if (do_adaptive_epsilon) {
            const CCTK_REAL vnorm = (SQR(A(0)) + SQR(A(1)) + SQR(A(2)) + SQR(A(3)) + SQR(A(4)));
               
            beta1 += 100.0*weno_eps*(vnorm + 1.0);
            beta2 += 100.0*weno_eps*(vnorm + 1.0);
            beta3 += 100.0*weno_eps*(vnorm + 1.0);
               
            const CCTK_REAL ibetanorm = 1.0 / (beta1 + beta2 + beta3);
               
            beta1 *= ibetanorm;
            beta2 *= ibetanorm;
            beta3 *= ibetanorm;
         }
         
         const CCTK_REAL wbarplus1 = 1.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarplus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarplus3 = 5.0/16.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3);
         
         const CCTK_REAL wplus1 = wbarplus1 * iwbarplussum;
         const CCTK_REAL wplus2 = wbarplus2 * iwbarplussum;
         const CCTK_REAL wplus3 = wbarplus3 * iwbarplussum;

         const CCTK_REAL wbarminus1 = 5.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarminus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarminus3 = 1.0/16.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarminussum = 1.0 / (wbarminus1 + wbarminus2 + wbarminus3);
         
         const CCTK_REAL wminus1 = wbarminus1 * iwbarminussum;
         const CCTK_REAL wminus2 = wbarminus2 * iwbarminussum;
         const CCTK_REAL wminus3 = wbarminus3 * iwbarminussum;
                                         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         Aminus(2) = 0;
         for (int n=0; n < 5; ++n) {
               Aplus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);
               Aminus(2) += (wminus1 * weno_coeffs[2][4-n]
                           + wminus2 * weno_coeffs[1][4-n]
                           + wminus3 * weno_coeffs[0][4-n]) * A(n);
         }
      }
   }
}

// instantiate all copies we need to avoid possible duplicates when multiple
// files need this and each instantiates their own copy
template class GRHydro_WENOReconstruct1d_cxx<false,false>;
template class GRHydro_WENOReconstruct1d_cxx<false,true >;
template class GRHydro_WENOReconstruct1d_cxx<true ,false>;

// the typedefs are required since the preprocessor chokes on the comma in the
// template specialization. Alternative:
// #define COMMA ,
// or #define / #undef a la Carpet's typecase.hh
typedef GRHydro_WENOReconstruct1d_cxx<false,false>
  GRHydro_WENOReconstruct1d_cxx00;
typedef GRHydro_WENOReconstruct1d_cxx<false,true >
  GRHydro_WENOReconstruct1d_cxx01;
typedef GRHydro_WENOReconstruct1d_cxx<true ,false>
  GRHydro_WENOReconstruct1d_cxx10;

INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_WENOReconstruct1d_cxx00)
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_WENOReconstruct1d_cxx01)
INSTANTIATE_RECONSTRUCTION_OPERATOR(GRHydro_WENOReconstruct1d_cxx10)
