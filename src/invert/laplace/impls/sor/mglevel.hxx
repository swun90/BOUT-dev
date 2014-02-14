
class MGLevel;
class MGField;

#ifndef __MGLEVEL_H__
#define __MGLEVEL_H__

#include <bout/mesh.hxx>

class MGLevel {
public:
  MGLevel(Mesh *m, int level);
  MGLevel(const MGLevel &m, int level);
  
  const BoutReal& dx(int jx, int jz) const {return _dx[jx][jz];}
  const BoutReal& dz() const {return _dz;}
  
  int xSize() const {return nx;}
  int ySize() const {return ny;}
  int zSize() const {return nz;}

  void communicate(MGField &f);
private:
  int nx, ny, nz;
  BoutReal **_dx, _dz;
  
  void average(BoutReal **in, int xsize, int zsize,
               BoutReal **out, int xfac, int zfac);
  
  void sum(BoutReal **in, int xsize, int zsize,
           BoutReal **out, int xfac, int zfac);
  
  int findFactor(int n, int factor);
};

class MGField {
public:
  MGField(MGLevel *m);
  ~MGField();
  
  BoutReal& operator()(int jx, int jy, int jz) {return data[jx][jz];}
  const BoutReal& operator()(int jx, int jy, int jz) const {return data[jx][jz];}
  
  // Assignment. Handles coarsening or refinement
  MGField & operator=(const MGField &rhs);
  
  void communicate() {fieldmesh->communicate(*this);}
  
  MGLevel* mesh() const {return fieldmesh;}
  
private:
  MGField(); // Prevents construction without mesh
  
  MGLevel *fieldmesh;
  BoutReal **data;
  
};

#endif // __MGLEVEL_H__
