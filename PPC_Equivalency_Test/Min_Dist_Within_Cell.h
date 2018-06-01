// Calculates the minimum distances within a cell.

#pragma once

#include <vector>

#include "Sort_Min_Distances.h"

using namespace std;

void Min_Dist_Within_Cell ( vector<Atom>const& cloud, vector<pair<double, pair<Atom, Atom>>>& min_dist )
{
    size_t cloud_size = cloud.size();
    
    for (unsigned counter_1 = 0; counter_1 < cloud_size; ++counter_1)
    {
        for (unsigned counter_2 = counter_1 + 1; counter_2 < cloud_size; ++counter_2)
        {
            double distance = norm( cloud[counter_1].cart_coords - cloud[counter_2].cart_coords );
            
            if (distance < min_dist[0].first)
            {
                min_dist[0] = pair<double, pair<Atom, Atom>>( distance, pair<Atom, Atom>( cloud[counter_1], cloud[counter_2] ) );
                
                sort( min_dist.begin(), min_dist.end(), Sort_Min_Distances );
            }
        }
    }
}
