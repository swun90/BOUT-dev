/*
  Finite-volume discretisation methods. Flux-conservative form
  
  NOTE: EXPERIMENTAL
 */

#ifndef __FV_OPS_H__
#define __FV_OPS_H__

#include "../field3d.hxx"
#include "../vector2d.hxx"

const Field3D Div_n_bxGrad_f_B_XPPM(const Field3D &n, const Field3D &f);

#endif // __FV_OPS_H__
