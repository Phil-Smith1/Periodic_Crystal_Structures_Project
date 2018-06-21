#pragma once

#include <string>
#include <vector>

using namespace std;

#include "Min_Cell_Length.h"

void Extend_Cell ( double min_dist, int& min_cell_length, vector<pair<string, double>>& cell_shape, int * cell_length_scale_factor )
{
    while (cell_shape[min_cell_length].second * cell_length_scale_factor[min_cell_length] < min_dist)
    {
        while (cell_shape[min_cell_length].second * cell_length_scale_factor[min_cell_length] < min_dist)
        {
            ++cell_length_scale_factor[min_cell_length];
        }
        
        min_cell_length = Min_Cell_Length(cell_shape[0].second * cell_length_scale_factor[0], cell_shape[1].second * cell_length_scale_factor[1], cell_shape[2].second * cell_length_scale_factor[2]);
    }
}
