#include <gemmi/cif.hpp>

#include <fstream>

#include "Read_Cell_Shape.h"
#include "Read_Atom_Coords.h"
#include "Obtain_Cloud_Of_T2_Centres.h"
#include "Min_Dist_Within_Cell.h"
#include "Extend_Cell_And_Cloud.h"
#include "Min_Dist_Within_Nbhd_Cells.h"
#include "Print_Summary.h"

const int num_distances = 2;

const int iterations = 5688;

const string dataset_directory = "/Users/philsmith/Documents/Xcode Projects/PPC_Equivalency_Test/T2_Dataset/";

bool Increasing_Distance ( pair<int, Point2d>const& p1, pair<int, Point2d>const& p2 )
{
    return p2.second.x > p1.second.x;
}

int main ( int, char*[] )
{
    clock_t start_time = clock(); // Start time.
    
    vector<pair<int, Point2d>> output_cloud;
    
    for (int counter_1 = 1; counter_1 < iterations + 1; ++counter_1)
    {
        string filename = "T2_" + to_string( counter_1 ) + "_num_molGeom.cif";
        string file_path = dataset_directory + filename;
        
        cif::Document cif_file = cif::read_file( file_path ); // Accessing CIF file.
        
        cif::Block * block = &(*(++cif_file.blocks.begin())); // Pointer to relevant block.
        
        vector<pair<string, double>> cell_shape;
        
        Read_Cell_Shape( block, cell_shape ); // Reading cell shape.
        
        int min_cell_length = Min_Cell_Length( cell_shape[0].second, cell_shape[1].second, cell_shape[2].second );
        
        double ** matrix;
        matrix = new double *[3];
        for (int counter_2 = 0; counter_2 < 3; ++counter_2)
        {
            matrix[counter_2] = new double [3];
        }
        
        Transformation_Matrix( cell_shape, matrix ); // Calculating the transformation matrix.
        
        vector<Atom> atom_cloud;
        
        Read_Atom_Coords( block, matrix, atom_cloud ); // Obtaining the atom cloud.
        
        vector<Atom> T2_centre_cloud;
        
        Obtain_Cloud_Of_T2_Centres( atom_cloud, matrix, T2_centre_cloud ); // Obtaining the cloud of molecule centres.
        
        vector<pair<double, pair<Atom, Atom>>> min_dist;
        
        size_t unit_cloud_size = T2_centre_cloud.size();
        
        Min_Dist_Within_Cell( T2_centre_cloud, unit_cloud_size, num_distances, min_dist ); // Calculating the minimum distances within the cell.
        
        vector<Atom> extended_cloud;
        bool extension = false;
        int cell_length_scale_factor[3] = { 1, 1, 1 };
        
        while (min_dist[0].first > 1e9) // Not enough distances.
        {
            extension = true;
            
            Extend_Cell_And_Cloud_2( cell_length_scale_factor, matrix, T2_centre_cloud, extended_cloud );
            
            Min_Dist_Within_Cell( extended_cloud, unit_cloud_size, num_distances, min_dist );
        }
        
        if (min_dist[0].first > cell_shape[min_cell_length].second) // Extending the cloud if necessary.
        {
            extension = true;
            
            Extend_Cell_And_Cloud( min_dist[0].first, min_cell_length, cell_shape, cell_length_scale_factor, matrix, T2_centre_cloud, extended_cloud );
            
            Min_Dist_Within_Cell( extended_cloud, unit_cloud_size, num_distances, min_dist );
        }
        
        // Calculating the minimum distances including neighbouring cells.
        
        if (extension) Min_Dist_Within_Nbhd_Cells( cell_length_scale_factor, matrix, extended_cloud, num_distances, min_dist );
        
        else Min_Dist_Within_Nbhd_Cells( cell_length_scale_factor, matrix, T2_centre_cloud, num_distances, min_dist );
        
        output_cloud.push_back( pair<int, Point2d>( counter_1, Point2d( min_dist[1].first, min_dist[0].first ) ) );
        
        Print_Info( counter_1, min_dist );
    }
    
    sort( output_cloud.begin(), output_cloud.end(), Increasing_Distance );
    
    ofstream ofs( "/Users/philsmith/Documents/Xcode Projects/PPC_Equivalency_Test/Results/ppc" + to_string( iterations ) + ".txt" );
    
    for (int counter = 0; counter < output_cloud.size(); ++counter)
    {
        ofs << "T2_" + to_string( output_cloud[counter].first ) + "_num_molGeom.cif " << output_cloud[counter].second.x << " " << output_cloud[counter].second.y << endl;
    }
    
    Print_Summary( start_time ); // Printing the summary.
    
    return 0;
}
