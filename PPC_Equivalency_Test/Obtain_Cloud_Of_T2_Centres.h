#pragma once

#include <vector>

#include "Frac_To_Cart_Coords.h"

using namespace std;

void Obtain_Cloud_Of_T2_Centres ( vector<Atom>const& atom_cloud, double ** matrix, vector<Atom>& T2_centre_cloud )
{
    int num_molecules = atom_cloud.size() / (double)46;
    
    T2_centre_cloud.reserve( num_molecules );
    
    for (int counter = 0; counter < num_molecules; ++counter)
    {
        Atom atom_1 = atom_cloud[9 + counter * 23], atom_2 = atom_cloud[22 + counter * 23];
        
        T2_centre_cloud.push_back( Atom( counter, (atom_1.frac_coords + atom_2.frac_coords) * 0.5 ) );
        
        Frac_To_Cart_Coords( matrix, T2_centre_cloud[counter] );
    }
}
