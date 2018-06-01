#include "Atom.h"

Atom::Atom ( int i, Point3d f_c )
{
    index = i;
    frac_coords = f_c;
}

Atom::Atom(){}
Atom::~Atom(){}
