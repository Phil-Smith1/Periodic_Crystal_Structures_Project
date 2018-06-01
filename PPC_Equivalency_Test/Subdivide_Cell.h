#pragma once

#include <vector>

#include "Atom.h"

using namespace std;

void Subdivide_Cell ( int factor, vector<Atom>const& atom_data, vector<pair<vector<Atom>, int>>& subcells )
{
    subcells.resize( pow( factor, 3 ) );
    
    for (unsigned counter = 0; counter < subcells.size(); ++counter)
    {
        subcells[counter].second = counter;
    }
    
    for (auto a : atom_data)
    {
        int box[3] = { factor - 1, factor - 1, factor - 1 };
        
        for (unsigned counter = 0; counter + 1 < factor; ++counter)
        {
            if (a.frac_coords.x < (counter + 1) / (double)factor)
            {
                box[0] = counter;
                break;
            }
        }
        
        for (unsigned counter = 0; counter + 1 < factor; ++counter)
        {
            if (a.frac_coords.y < (counter + 1) / (double)factor)
            {
                box[1] = counter;
                break;
            }
        }
        
        for (unsigned counter = 0; counter + 1 < factor; ++counter)
        {
            if (a.frac_coords.z < (counter + 1) / (double)factor)
            {
                box[2] = counter;
                break;
            }
        }
        
        subcells[factor * factor * box[2] + factor * box[1] + box[0]].first.push_back( a );
    }
}
