
#include "mglevel.hxx"

#include <utils.hxx>

MGLevel::MGLevel(Mesh *m, int level) {
  // Start from a mesh, and reduce each dimension by (up to) 2^level
  
  nx = findFactor(m->xend - m->xstart + 1, level);
  ny = m->ngy;
  nz = findFactor(m->ngz-1, level);
  
  int xfactor = (m->xend - m->xstart + 1) / nx;
  int zfactor = (m->ngz-1) / nz;
  
  _dz = m->dz * zfactor;
  
}

void MGLevel::communicate(MGField &f) {
  // Check that this is the correct mesh
  if(f.mesh() != this)
    throw BoutException("Incorrect MGLevel used to communicate MGField");
  
  // Not implemented yet
}

int findFactor(int n, int level) {
  for(int i=0;i<level;i++) {
    if(n % 2 != 0) {
      // n not divisible by 2 -> Don't divide
      break;
    }
    
    n /= 2;
  }
  return n;
}


////////////////////////////////////////////////////////


MGField::MGField(MGLevel *m) : fieldmesh(m) {
  // Allocate memory
  
  data = rmatrix(m->xSize(), m->zSize());
  
  
}

MGField::~MGField() {
  // Free memory
  free_rmatrix(data);
}

