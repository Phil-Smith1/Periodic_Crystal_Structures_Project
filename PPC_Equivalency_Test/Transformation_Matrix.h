#pragma once

#include "Cell_Volume.h"

void Transformation_Matrix ( vector<pair<string, double>>const& cell_shape, double **matrix )
{
    matrix[0][0] = cell_shape[0].second;
    matrix[0][1] = cell_shape[1].second * cos( pi * cell_shape[5].second / (double)180 );
    matrix[0][2] = cell_shape[2].second * cos( pi * cell_shape[4].second / (double)180 );
    matrix[1][0] = 0;
    matrix[1][1] = cell_shape[1].second * sin( pi * cell_shape[5].second / (double)180 );
    matrix[1][2] = cell_shape[2].second * (cos( pi * cell_shape[3].second / (double)180 ) - cos( pi * cell_shape[4].second / (double)180 ) * cos( pi * cell_shape[5].second / (double)180 )) / (double)sin( pi * cell_shape[5].second / (double)180 );
    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = Cell_Volume( cell_shape ) / (double)(cell_shape[0].second * cell_shape[1].second * sin( pi * cell_shape[5].second / (double)180 ));
}
