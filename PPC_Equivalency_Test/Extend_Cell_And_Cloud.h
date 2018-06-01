#pragma once

#include "Extend_Cell.h"
#include "Transformation_Matrix.h"
#include "Extend_Cloud.h"

void Extend_Cell_And_Cloud ( double dist, int& min_cell_length, vector<pair<string, double>>& cell_shape, double **matrix, vector<Atom>const& cloud, vector<Atom>& extended_cloud )
{
    int cell_length_scale_factor[3] = { 1, 1, 1 };
    
    Extend_Cell( dist, min_cell_length, cell_shape, cell_length_scale_factor );
    
    Transformation_Matrix( cell_shape, matrix );
    
    Extend_Cloud( cell_length_scale_factor, matrix, cloud, extended_cloud );
}
