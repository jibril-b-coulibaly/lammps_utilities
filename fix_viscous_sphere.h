/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   adaptation of fix viscous for atom_style sphere, applied on the angular velocity
   written by Jibril B. Coulibaly @ Northwestern University, 04/12/2019
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(viscous/sphere,FixViscousSphere)

#else

#ifndef LMP_FIX_VISCOUS_SPHERE_H
#define LMP_FIX_VISCOUS_SPHERE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixViscousSphere : public Fix {
 public:
  FixViscousSphere(class LAMMPS *, int, char **);
  virtual ~FixViscousSphere();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 protected:
  double *gamma;
  int ilevel_respa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
