#include "bout/fv_ops.hxx"

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <derivs.hxx>
#include <output.hxx>
#include <utils.hxx>

#include <cmath>

///////////////////////////////////////////////////////////////////
// Finite volume methods

BoutReal BOUTMIN(const BoutReal &a, const BoutReal &b, const BoutReal &c, const BoutReal &d) {
  BoutReal r1 = (a < b) ? a : b;
  BoutReal r2 = (c < d) ? c : d;
  return (r1 < r2) ? r1 : r2;
}

struct Stencil1D {
  // Cell centre values
  BoutReal c, m, p, mm, pp;
  
  // Left and right cell face values
  BoutReal L, R;
};

// First order upwind for testing
void Upwind(Stencil1D &n, const BoutReal h) {
  n.L = n.R = n.c;
}

// Fromm method
void Fromm(Stencil1D &n, const BoutReal h) {
  n.L = n.c - 0.25*(n.p - n.m);
  n.R = n.c + 0.25*(n.p - n.m);
}

void XPPM(Stencil1D &n, const BoutReal h) {
  // 4th-order PPM interpolation in X
  
  const BoutReal C = 1.25; // Limiter parameter
  
  BoutReal h2 = h*h;
  
  n.R = (7./12)*(n.c + n.p) - (1./12)*(n.m + n.pp);
  n.L = (7./12)*(n.c + n.m) - (1./12)*(n.mm + n.p);
  
  // Apply limiters
  if( (n.c - n.R)*(n.p - n.R) > 0.0 ) {
    // Calculate approximations to second derivative
    
    BoutReal D2 = (3./h2)*(n.c - 2*n.R + n.p);
    BoutReal D2L = (1./h2)*(n.m - 2*n.c + n.p);
    BoutReal D2R = (1./h2)*(n.c - 2.*n.p + n.pp);
    
    BoutReal D2lim; // Value to be used in limiter
    
    // Check if they all have the same sign
    if( (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
      // Same sign
	    
      D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), fabs(D2) );
    }else {
      // Different sign
      D2lim = 0.0;
    }
    
    n.R = 0.5*(n.c + n.p) - (h2/6)*D2lim;
  }
  
  if( (n.m - n.L)*(n.c - n.L) > 0.0 ) {
    // Calculate approximations to second derivative
    
    BoutReal D2 = (3./h2)*(n.m - 2*n.L + n.c);
    BoutReal D2L = (1./h2)*(n.mm - 2*n.m + n.c);
    BoutReal D2R = (1./h2)*(n.m - 2.*n.c + n.p);
    
    BoutReal D2lim; // Value to be used in limiter
    
    // Check if they all have the same sign
    if( (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
      // Same sign
      
      D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), fabs(D2) );
    }else {
      // Different sign
      D2lim = 0.0;
    }
    
    n.L = 0.5*(n.m + n.c) - (h2/6)*D2lim;
  }
  
  if( ( (n.R - n.c)*(n.c - n.L) <= 0.0 ) || ( (n.m - n.c)*(n.c - n.p) <= 0.0 ) ) {
    // At a local maximum or minimum
    
    BoutReal D2 = (6./h2)*(n.L - 2.*n.c + n.R);
    
    if(fabs(D2) < 1e-10) {
      n.R = n.L = n.c;
    }else {
      BoutReal D2C = (1./h2)*(n.m - 2.*n.c + n.p);
      BoutReal D2L = (1./h2)*(n.mm - 2*n.m + n.c);
      BoutReal D2R = (1./h2)*(n.c - 2.*n.p + n.pp);
    
      BoutReal D2lim;
      // Check if they all have the same sign
      if( (D2*D2C > 0.0) && (D2*D2L > 0.0) && (D2*D2R > 0.0) ) {
	// Same sign
	
	D2lim = SIGN(D2) * BOUTMIN( C*fabs(D2L), C*fabs(D2R), C*fabs(D2C), fabs(D2) );
	n.R = n.c + (n.R - n.c)*D2lim / D2;
	n.L = n.c + (n.L - n.c)*D2lim / D2;
      }else {
	// Different signs
	n.R = n.L = n.c;
      }
    }
  }
}

// Div (n * b x Grad(f)/B)
const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n, const Field3D &f) {
  Field3D result = 0;
  
  //////////////////////////////////////////
  // X-Z advection.
  // 
  //             Z
  //             |
  // 
  //    fmp --- vU --- fpp
  //     |      nU      |
  //     |               |
  //    vL nL        nR vR    -> X
  //     |               |
  //     |      nD       |
  //    fmm --- vD --- fpm
  //
  
  for(int i=mesh->xstart;i<=mesh->xend;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++) {
	int kp = (k+1) % (mesh->ngz-1);
	int kpp = (kp+1) % (mesh->ngz-1);
	int km = (k-1+mesh->ngz-1) % (mesh->ngz-1);
	int kmm = (km-1+mesh->ngz-1) % (mesh->ngz-1);
	
	// 1) Interpolate stream function f onto corners fmp, fpp, fpm
	
	BoutReal fmm = 0.25*(f(i,j,k) + f(i-1,j,k) + f(i,j,km) + f(i-1,j,km));
	BoutReal fmp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i-1,j,k) + f(i-1,j,kp)); // 2nd order accurate
	BoutReal fpp = 0.25*(f(i,j,k) + f(i,j,kp) + f(i+1,j,k) + f(i+1,j,kp));
	BoutReal fpm = 0.25*(f(i,j,k) + f(i+1,j,k) + f(i,j,km) + f(i+1,j,km));
	
	// 2) Calculate velocities on cell faces
	
	BoutReal vU = mesh->J(i,j)*(fmp - fpp)/mesh->dx(i,j); // -J*df/dx
	BoutReal vD = mesh->J(i,j)*(fmm - fpm)/mesh->dx(i,j); // -J*df/dx
	
	BoutReal vR = 0.5*(mesh->J(i,j)+mesh->J(i+1,j))*(fpp - fpm)/mesh->dz; // J*df/dz 
	BoutReal vL = 0.5*(mesh->J(i,j)+mesh->J(i-1,j))*(fmp - fmm)/mesh->dz; // J*df/dz 
	
	// 3) Calculate n on the cell faces. The sign of the
	//    velocity determines which side is used.
	
	// X direction
	Stencil1D s;
	s.c  = n(i,  j,k);
	s.m  = n(i-1,j,k);
	s.mm = n(i-2,j,k);
	s.p  = n(i+1,j,k);
	s.pp = n(i+2,j,k);
	
	//Upwind(s, mesh->dx(i,j));
	//Fromm(s, mesh->dx(i,j));
	XPPM(s, mesh->dx(i,j)); 
	
	if(vR > 0.0) {
	  // Flux out into next cell
	  BoutReal flux = vR * s.R;
	  result(i,j,k)   += flux / (mesh->dx(i,j) * mesh->J(i,j));
	  result(i+1,j,k) -= flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
	}
	
	if(vL < 0.0) {
	  BoutReal flux = vL * s.L;
	  result(i,j,k)   -= flux / (mesh->dx(i,j) * mesh->J(i,j));
	  result(i-1,j,k) += flux / (mesh->dx(i+1,j) * mesh->J(i+1,j));
	}

	/// NOTE: Need to communicate fluxes

	// Z direction
	s.m  = n(i,j,km);
	s.mm = n(i,j,kmm);
	s.p  = n(i,j,kp);
	s.pp = n(i,j,kpp);
	
	//Upwind(s, mesh->dz);
	//Fromm(s, mesh->dz);
	XPPM(s, mesh->dz);

	if(vU > 0.0) {
	  BoutReal flux = vU * s.R / mesh->dz;
	  result(i,j,k)   += flux;
	  result(i,j,kp)  -= flux;
	}
	if(vD < 0.0) {
	  BoutReal flux = vD * s.L / mesh->dz;
	  result(i,j,k)   -= flux;
	  result(i,j,km)  += flux;
	}
	
      }
  
  return result;
}

