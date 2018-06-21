// Calculating the cartesian coordinates given the fractional coordinates and the transformation matrix.

#pragma once

#include "Atom.h"

void Frac_To_Cart_Coords ( double ** matrix, Atom& atom )
{
    atom.cart_coords.x = matrix[0][0] * atom.frac_coords.x + matrix[0][1] * atom.frac_coords.y + matrix[0][2] * atom.frac_coords.z;
    atom.cart_coords.y = matrix[1][0] * atom.frac_coords.x + matrix[1][1] * atom.frac_coords.y + matrix[1][2] * atom.frac_coords.z;
    atom.cart_coords.z = matrix[2][0] * atom.frac_coords.x + matrix[2][1] * atom.frac_coords.y + matrix[2][2] * atom.frac_coords.z;
}
