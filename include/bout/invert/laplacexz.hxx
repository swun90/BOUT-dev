/**************************************************************************
 * Perpendicular Laplacian inversion in X-Z
 *
 * Equation solved is:
 *
 * Div( A * Grad_perp(x) ) + B*x = b
 *
 *
 **************************************************************************
 * Copyright 2015 B.D.Dudson
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __LAPLACEXZ_H__
#define __LAPLACEXZ_H__

#include <options.hxx>
#include <field3d.hxx>
#include <bout/mesh.hxx>

class LaplaceXZ {
public:
  LaplaceXZ(Mesh *m, Options *options) {}
  virtual ~LaplaceXZ() {}

  virtual void setCoefs(const Field2D &A, const Field2D &B) = 0;
  virtual void setCoefs(const Field3D &A, const Field3D &B) { setCoefs(A.DC(), B.DC()); }

  virtual Field3D solve(const Field3D &b, const Field3D &x0) = 0;

  static LaplaceXZ* create(Mesh *m, Options *opt = NULL);
private:

};

#endif // __LAPLACEXZ_H__
