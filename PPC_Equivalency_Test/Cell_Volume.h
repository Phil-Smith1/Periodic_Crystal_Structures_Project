// Calculating the volume of a cell from the cell shape.

#pragma once

#include <cmath>
#include <string>
#include <vector>

using namespace std;

#ifndef PI
#define PI 1
const double pi = 3.14159265359;
#endif

double Cell_Volume ( vector<pair<string, double>>const& cell_shape )
{
    double volume = cell_shape[0].second * cell_shape[1].second * cell_shape[2].second * sqrt( 1 - pow( cos( pi * cell_shape[3].second / (double)180 ), 2 ) - pow( cos( pi * cell_shape[4].second / (double)180 ), 2 ) - pow( cos( pi * cell_shape[5].second / (double)180 ), 2 ) + 2 * cell_shape[1].second * cell_shape[2].second * cos( pi * cell_shape[3].second / (double)180 ) * cos( pi * cell_shape[4].second / (double)180 ) * cos( pi * cell_shape[5].second / (double)180 ) );
    
    return volume;
}
