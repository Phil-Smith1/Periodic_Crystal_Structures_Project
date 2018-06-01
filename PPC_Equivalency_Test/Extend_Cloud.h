#pragma once

#include "Frac_To_Cart_Coords.h"

void Extend_Cloud ( int *cell_length_scale_factor, double **matrix, vector<Atom>const& cloud, vector<Atom>& extended_cloud )
{
    int product_scale_factor = cell_length_scale_factor[0] * cell_length_scale_factor[1] * cell_length_scale_factor[2];
    
    size_t cloud_size = cloud.size();
    
    extended_cloud.reserve( product_scale_factor * cloud_size );
    
    unsigned counter_1 = 0;
    
    for (unsigned counter_2 = 0; counter_2 < cell_length_scale_factor[2]; ++counter_2)
    {
        for (unsigned counter_3 = 0; counter_3 < cell_length_scale_factor[1]; ++counter_3)
        {
            for (unsigned counter_4 = 0; counter_4 < cell_length_scale_factor[0]; ++counter_4)
            {
                for (auto c : cloud)
                {
                    extended_cloud.push_back( Atom( counter_1, c.frac_coords + Point3d( counter_4, counter_3, counter_2 ) ) );
                    Frac_To_Cart_Coords( matrix, extended_cloud[counter_1] );
                    cout << extended_cloud[counter_1].frac_coords << endl;
                    ++counter_1;
                }
            }
        }
    }
}
