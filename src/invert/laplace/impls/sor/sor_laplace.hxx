/**************************************************************************
 * Perpendicular Laplacian inversion using SOR method
 * 
 *
 **************************************************************************
 * Copyright 2014 B.D.Dudson
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

class LaplaceSOR;

#ifndef __LAP_SOR_H__
#define __LAP_SOR_H__

#include <invert_laplace.hxx>
#include <options.hxx>
#include <boutexception.hxx>

class LaplaceSOR : public Laplacian {
public:
  LaplaceSOR(Options *opt = NULL);
  ~LaplaceSOR();
  
  void setCoefA(const Field3D &val) { A = val; }
  void setCoefC(const Field3D &val) { C = val; }
  void setCoefD(const Field3D &val) { D = val; }
  
  void setCoefA(const Field2D &val) { A = val; }
  void setCoefC(const Field2D &val) { C = val; }
  void setCoefD(const Field2D &val) { D = val; }

  void setCoefEx(const Field2D &val) { throw BoutException("LaplaceSOR does not have Ex coefficient"); }
  void setCoefEz(const Field2D &val) { throw BoutException("LaplaceSOR does not have Ez coefficient"); }
  
  const FieldPerp solve(const FieldPerp &b) {return solve(b,b);}
  const FieldPerp solve(const FieldPerp &b, const FieldPerp &x0);
  
private:
  Field3D A, C, D;
  
  
  
};

#endif
