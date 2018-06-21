#pragma once

#include <gemmi/cif.hpp>

#include <string>
#include <vector>

using namespace std;
namespace cif = gemmi::cif;

void Read_Cell_Shape ( cif::Block * block, vector<pair<string, double>>& cell_shape )
{
    vector<string> labels;
    labels.reserve( 6 );
    
    labels.push_back( "_cell_length_a" );
    labels.push_back( "_cell_length_b" );
    labels.push_back( "_cell_length_c" );
    labels.push_back( "_cell_angle_alpha" );
    labels.push_back( "_cell_angle_beta" );
    labels.push_back( "_cell_angle_gamma" );
    
    cif::Table table = block->find( labels );
    cif::Table::Row row = table.operator[]( 0 );
    
    cell_shape.reserve( 6 );
    
    for (int counter = 0; counter < 6; ++counter)
    {
        double value = stod( row.operator[]( counter ) );
        cell_shape.push_back( pair<string, double>( labels[counter], value ) );
    }
}
