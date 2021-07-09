// Calculates the minimum distances within a cell.

#pragma once

#include <vector>

#include "Sort_Min_Distances.h"

using namespace std;

void Min_Dist_Within_Cell ( vector<Atom>const& cloud, size_t unit_cloud_size, int num_distances, vector<pair<double, pair<Atom, Atom>>>& min_dist )
{
    min_dist.clear();
    min_dist.reserve( num_distances );
    
    for (int counter = 0; counter < num_distances; ++counter)
    {
        min_dist.push_back( pair<double, pair<Atom, Atom>>( 1e10, pair<Atom, Atom>( cloud[0], cloud[0] ) ) );
    }
    
    size_t cloud_size = cloud.size();
    
    for (int counter_1 = 0; counter_1 < unit_cloud_size; ++counter_1)
    {
        for (int counter_2 = counter_1 + 1; counter_2 < cloud_size; ++counter_2)
        {
            double distance = norm( cloud[counter_1].cart_coords - cloud[counter_2].cart_coords );
            
            if (distance < min_dist[0].first + 1e-10)
            {
                min_dist[0] = pair<double, pair<Atom, Atom>>( distance, pair<Atom, Atom>( cloud[counter_1], cloud[counter_2] ) );
                
                sort( min_dist.begin(), min_dist.end(), Sort_Min_Distances );
            }
        }
    }
}
