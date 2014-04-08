/*
 * Passive advection test
 *
 */

#include <bout/physicsmodel.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>

#include "bout/fv_ops.hxx"

class Advect : public PhysicsModel {
public:
  int init(bool restarting) {
    
    initial_profile("psi", psi);
    SAVE_ONCE(psi);
    mesh->communicate(psi);
    
    solver->add(n, "n");
    
    return 0;
  }
  
  int rhs(BoutReal t) {
    mesh->communicate(n);
    
    ddt(n) = -Div_n_bxGrad_f_B_XPPM(n, psi);
    //ddt(n) = bracket(n, psi, BRACKET_ARAKAWA);
    
    //ddt(n) += 10.*(SQ(mesh->dx)*D2DX2(n) + SQ(mesh->dz)*D2DZ2(n));
    //ddt(n) -= 10.*(SQ(SQ(mesh->dx))*D4DX4(n) + SQ(SQ(mesh->dz))*D4DZ4(n));
    
    return 0;
  }
  
private:
  Field3D n, psi;
};

BOUTMAIN(Advect);
