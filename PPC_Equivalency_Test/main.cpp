#include <gemmi/cif.hpp>

#include <iostream>

#include "Read_Cell_Shape.h"
#include "Read_Atom_Coords.h"
#include "Min_Dist_Within_Cell.h"
#include "Extend_Cell_And_Cloud.h"
#include "Min_Dist_Within_Nbhd_Cells.h"
#include "Print_Summary.h"

namespace cif = gemmi::cif;

const string cif_file_name = "/Users/philsmith/Documents/Xcode Projects/Cif/T2_1_num_molGeom.cif";

const int num_distances = 2;

int main ( int, char*[] )
{
    // Start time.
    
    clock_t start_time = clock();
    
    // Accessing relevant block.
    
    cif::Document cif_file = cif::read_file( cif_file_name );
    
    cif::Block *block = cif_file.find_block( "-223.695_t2_4_4211" );
    
    // Reading cell shape.
    
    vector<pair<string, double>> cell_shape;
    
    Read_Cell_Shape( block, cell_shape );
    
    int min_cell_length = Min_Cell_Length( cell_shape[0].second, cell_shape[1].second, cell_shape[2].second );
    
    // Calculating the transformation matrix.
    
    double **matrix;
    matrix = new double *[3];
    for (unsigned counter = 0; counter < 3; ++counter)
    {
        matrix[counter] = new double [3];
    }
    
    Transformation_Matrix( cell_shape, matrix );
    
    // Obtaining the cloud.
    
    vector<Atom> cloud;
    
    Read_Atom_Coords( block, matrix, cloud );
    
    vector<Atom> reduced_cloud;
    reduced_cloud.reserve(2);
    
    Atom atom_1 = cloud[9], atom_2 = cloud[22];
    
    reduced_cloud.push_back( Atom( 0, (atom_1.frac_coords + atom_2.frac_coords) * 0.5 ) );
    
    atom_1 = cloud[32];
    atom_2 = cloud[45];
    
    reduced_cloud.push_back( Atom( 1, (atom_1.frac_coords + atom_2.frac_coords) * 0.5 ) );
    
    Frac_To_Cart_Coords( matrix, reduced_cloud[0] );
    Frac_To_Cart_Coords( matrix, reduced_cloud[1] );
    
    // Calculating the minimum distances within the cell.
    
    vector<pair<double, pair<Atom, Atom>>> min_dist;
    min_dist.reserve( num_distances );
    for (unsigned counter = 0; counter < num_distances; ++counter)
    {
        min_dist.push_back( pair<double, pair<Atom, Atom>>( 1e10, pair<Atom, Atom>( cloud[0], cloud[0] ) ) );
    }
    
    Min_Dist_Within_Cell( reduced_cloud, min_dist );
    
    // Extending the cloud if necessary.
    
    vector<Atom> extended_cloud;
    bool extension = false;
    
    if (min_dist[0].first > cell_shape[min_cell_length].second)
    {
        min_dist[0].first = min_dist[1].first + 10;
        
        extension = true;
        
        Extend_Cell_And_Cloud( min_dist[0].first, min_cell_length, cell_shape, matrix, reduced_cloud, extended_cloud );
        
        Min_Dist_Within_Cell( extended_cloud, min_dist );
    }
    
    // Calculating the minimum distances including neighbouring cells.
    
    if (extension) Min_Dist_Within_Nbhd_Cells( cell_shape, matrix, extended_cloud, min_dist );
    
    else Min_Dist_Within_Nbhd_Cells( cell_shape, matrix, reduced_cloud, min_dist );
    
    for (auto iter = min_dist.rbegin(); iter != min_dist.rend(); ++iter)
    {
        cout << iter->first << endl;
    }
    
    cout << endl;
    
    cout << "The minimum distance is " << min_dist.rbegin()->first << "." << endl << endl;
    
    // Printing the summary.
    
    Print_Summary( start_time );
    
    return 0;
}
