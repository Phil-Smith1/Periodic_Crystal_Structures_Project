#pragma once

#include <string>
#include <vector>

#include "Frac_To_Cart_Coords.h"

using namespace std;

void Min_Dist_Within_Nbhd_Cells ( int * cell_length_scale_factor, double ** matrix, vector<Atom>const& cloud, int num_distances, vector<pair<double, pair<Atom, Atom>>>& min_dist )
{
    vector<Point3d> cell_shift;
    cell_shift.reserve( 13 );
    
    cell_shift.push_back( Point3d( 1, -1, -1 ) );
    cell_shift.push_back( Point3d( 1, 0, -1 ) );
    cell_shift.push_back( Point3d( 1, 1, -1 ) );
    cell_shift.push_back( Point3d( 1, -1, 0 ) );
    cell_shift.push_back( Point3d( 1, 0, 0 ) );
    cell_shift.push_back( Point3d( 1, 1, 0 ) );
    cell_shift.push_back( Point3d( 1, -1, 1 ) );
    cell_shift.push_back( Point3d( 1, 0, 1 ) );
    cell_shift.push_back( Point3d( 1, 1, 1 ) );
    cell_shift.push_back( Point3d( 0, -1, 1 ) );
    cell_shift.push_back( Point3d( 0, 0, 1 ) );
    cell_shift.push_back( Point3d( 0, 1, 1 ) );
    cell_shift.push_back( Point3d( 0, 1, 0 ) );
    
    size_t cloud_size = cloud.size();
    
    for (unsigned counter_1 = 0; counter_1 < 13; ++counter_1)
    {
        vector<Atom> shifted_cloud;
        shifted_cloud.reserve( cloud_size );
        
        cell_shift[counter_1] = Point3d( cell_shift[counter_1].x * cell_length_scale_factor[0], cell_shift[counter_1].y * cell_length_scale_factor[1], cell_shift[counter_1].z * cell_length_scale_factor[2] );
        
        for (auto c : cloud)
        {
            shifted_cloud.push_back( Atom( c.index, c.frac_coords + cell_shift[counter_1] ) );
        }
        
        for (unsigned counter_2 = 0; counter_2 < cloud_size; ++counter_2)
        {
            Frac_To_Cart_Coords( matrix, shifted_cloud[counter_2] );
        }
        
        for (int counter_1 = 0; counter_1 < 2; ++counter_1)
        {
            for (auto s : shifted_cloud)
            {
                double distance = norm( cloud[counter_1].cart_coords - s.cart_coords );
                
                if (distance < min_dist[0].first)
                {
                    min_dist[0] = pair<double, pair<Atom, Atom>>( distance, pair<Atom, Atom>( cloud[counter_1], s ) );
                    
                    sort( min_dist.begin(), min_dist.end(), Sort_Min_Distances );
                }
            }
        }
    }
}
