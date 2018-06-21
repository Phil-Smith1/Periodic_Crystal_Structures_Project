#pragma once

#include <gemmi/cif.hpp>

#include <string>
#include <vector>

#include "Frac_To_Cart_Coords.h"

using namespace std;
namespace cif = gemmi::cif;

void Read_Atom_Coords ( cif::Block * block, double ** matrix, vector<Atom>& cloud )
{
    vector<string> labels;
    labels.push_back( "_atom_site_fract_x" );
    labels.push_back( "_atom_site_fract_y" );
    labels.push_back( "_atom_site_fract_z" );
    
    cif::Table table = block->find( labels );
    
    size_t num_atoms = table.length();
    
    cloud.reserve( num_atoms );
    
    for (int counter = 0; counter < num_atoms; ++counter)
    {
        cif::Table::Row row = table.operator[]( counter );
        
        double x_value = stod( row.operator[]( 0 ) ), y_value = stod( row.operator[]( 1 ) ), z_value = stod( row.operator[]( 2 ) );
        
        cloud.push_back( Atom( counter, Point3d( x_value, y_value, z_value ) ) );
        
        Frac_To_Cart_Coords( matrix, cloud[counter] );
    }
}
