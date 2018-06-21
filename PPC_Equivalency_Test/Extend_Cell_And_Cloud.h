#pragma once

#include "Extend_Cell.h"
#include "Transformation_Matrix.h"
#include "Extend_Cloud.h"

void Extend_Cell_And_Cloud ( double dist, int& min_cell_length, vector<pair<string, double>>& cell_shape, int *cell_length_scale_factor, double ** matrix, vector<Atom>const& cloud, vector<Atom>& extended_cloud )
{
    Extend_Cell( dist, min_cell_length, cell_shape, cell_length_scale_factor );
    
    Extend_Cloud( cell_length_scale_factor, matrix, cloud, extended_cloud );
}

void Extend_Cell_And_Cloud_2 ( int *cell_length_scale_factor, double **matrix, vector<Atom>const& cloud, vector<Atom>& extended_cloud )
{
    ++cell_length_scale_factor[0];
    
    Extend_Cloud( cell_length_scale_factor, matrix, cloud, extended_cloud );
}
