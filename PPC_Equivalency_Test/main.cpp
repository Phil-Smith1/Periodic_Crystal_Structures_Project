#include <gemmi/cif.hpp>

#include <fstream>

#include "Read_Cell_Shape.h"
#include "Read_Atom_Coords.h"
#include "Obtain_Cloud_Of_T2_Centres.h"
#include "Min_Dist_Within_Cell.h"
#include "Extend_Cell_And_Cloud.h"
#include "Min_Dist_Within_Nbhd_Cells.h"
#include "Print_Summary.h"

const int num_distances = 7764;

const int iterations = 1;

const string dataset_directory = "/Users/philsmith/Documents/Work/Xcode Projects/T2_Dataset/";

bool Increasing_Distance ( pair<int, Point2d>const& p1, pair<int, Point2d>const& p2 )
{
    return p2.second.x > p1.second.x;
}

int main ( int, char*[] )
{
    clock_t start_time = clock(); // Starts the time of the code.
    
    vector<pair<int, Point2d>> output_cloud;
    
    for (int counter_1 = 1; counter_1 < iterations + 1; ++counter_1)
    {
        string filename = "T2_" + to_string( counter_1 ) + "_num_molGeom.cif";
        string file_path = dataset_directory + filename;
        
        cif::Document cif_file = cif::read_file( file_path ); // Accessing CIF file.
        
        cif::Block * block = &(*(++cif_file.blocks.begin())); // Pointer to relevant block.
        
        vector<pair<string, double>> cell_shape;
        
        Read_Cell_Shape( block, cell_shape ); // Reading cell shape.
        
        cell_shape.clear();
        cell_shape.push_back( pair<string, double>( "_cell_length_a", 100 ) );
        cell_shape.push_back( pair<string, double>( "_cell_length_b", 100 ) );
        cell_shape.push_back( pair<string, double>( "_cell_length_c", 100 ) );
        cell_shape.push_back( pair<string, double>( "_cell_angle_alpha", 90 ) );
        cell_shape.push_back( pair<string, double>( "_cell_angle_beta", 90 ) );
        cell_shape.push_back( pair<string, double>( "_cell_angle_gamma", 90 ) );
        
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
        
        double u = -0.03;
        vector<Point3d> frac_base_pts;
        
        frac_base_pts.push_back( Point3d( u, 0, 0.25 ) );
        frac_base_pts.push_back( Point3d( -u, 0.5, 0.25 ) );
        frac_base_pts.push_back( Point3d( 0.5 - u, 0, 0.75 ) );
        frac_base_pts.push_back( Point3d( u + 0.5, 0.5, 0.75 ) );
        frac_base_pts.push_back( Point3d( 0.25, u, 0 ) );
        
        frac_base_pts.push_back( Point3d( 0.25, -u, 0.5 ) );
        
        frac_base_pts.push_back( Point3d( 0.75, 0.5 - u, 0 ) );
        frac_base_pts.push_back( Point3d( 0.75, u + 0.5, 0.5 ) );
        frac_base_pts.push_back( Point3d( 0, 0.25, u ) );
        frac_base_pts.push_back( Point3d( 0.5, 0.25, -u ) );
        frac_base_pts.push_back( Point3d( 0, 0.75, 0.5 - u ) );
        frac_base_pts.push_back( Point3d( 0.5, 0.75, u + 0.5 ) );
        frac_base_pts.push_back( Point3d( -u, 0, 0.75 ) );
        frac_base_pts.push_back( Point3d( u, 0.5, 0.75 ) );
        frac_base_pts.push_back( Point3d( u + 0.5, 0, 0.25 ) );
        frac_base_pts.push_back( Point3d( 0.5 - u, 0.5, 0.25 ) );
        frac_base_pts.push_back( Point3d( 0.75, -u, 0 ) );
        frac_base_pts.push_back( Point3d( 0.75, u, 0.5 ) );
        frac_base_pts.push_back( Point3d( 0.25, u + 0.5, 0 ) );
        frac_base_pts.push_back( Point3d( 0.25, 0.5 - u, 0.5 ) );
        frac_base_pts.push_back( Point3d( 0, 0.75, -u ) );
        frac_base_pts.push_back( Point3d( 0.5, 0.75, u ) );
        frac_base_pts.push_back( Point3d( 0, 0.25, u + 0.5 ) );
        frac_base_pts.push_back( Point3d( 0.5, 0.25, 0.5 - u ) );
        
        T2_centre_cloud.clear();
        
        for (int counter = 0; counter < frac_base_pts.size(); ++counter)
        {
            if (frac_base_pts[counter].x < -1e-10) frac_base_pts[counter].x = 1 + frac_base_pts[counter].x;
            if (frac_base_pts[counter].y < -1e-10) frac_base_pts[counter].y = 1 + frac_base_pts[counter].y;
            if (frac_base_pts[counter].z < -1e-10) frac_base_pts[counter].z = 1 + frac_base_pts[counter].z;
            
            Atom a( counter, frac_base_pts[counter] );
            
            Frac_To_Cart_Coords( matrix, a );
            
            T2_centre_cloud.push_back( a );
        }
        
        vector<pair<double, pair<Atom, Atom>>> min_dist;
        
        size_t unit_cloud_size = T2_centre_cloud.size();
        
        Min_Dist_Within_Cell( T2_centre_cloud, unit_cloud_size, num_distances, min_dist ); // Calculating the minimum distances within the cell.
        
        vector<Atom> extended_cloud;
        bool extension = false;
        int cell_length_scale_factor[3] = { 1, 1, 1 };
        
        /*while (min_dist[0].first > 1e9) // Not enough distances.
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
        
        else Min_Dist_Within_Nbhd_Cells( cell_length_scale_factor, matrix, T2_centre_cloud, num_distances, min_dist );*/
        
        Min_Dist_Within_Nbhd_Cells( cell_length_scale_factor, matrix, T2_centre_cloud, num_distances, min_dist );
        
        sort( min_dist.begin(), min_dist.end(), Sort_Min_Distances_2 );
        
        //output_cloud.push_back( pair<int, Point2d>( counter_1, Point2d( min_dist[1].first, min_dist[0].first ) ) );
        
        ofstream ofs( "/Users/philsmith/Documents/Work/Xcode Projects/PPC_Equivalency_Test/Results/ppc" + to_string( iterations ) + ".txt" );
        
        for (int counter = 0; counter < min_dist.size(); ++counter)
        {
            //ofs << min_dist[counter].first << endl;
            
            Point3d v = min_dist[num_distances - 1 - counter].second.first.cart_coords - min_dist[num_distances - 1 - counter].second.second.cart_coords;
            
            if (abs( v.x ) < tiny_num)
            {
                if (abs( v.y ) < tiny_num)
                {
                    if (v.z < -tiny_num) v.z *= -1;
                    
                    while (v.z > 100 - tiny_num) v.z -= 100;
                }
                
                else
                {
                    if (v.y < -tiny_num)
                    {
                        v.y *= -1;
                        v.z *= -1;
                    }
                    
                    while (v.y > 100 - tiny_num) v.y -= 100;
                    
                    while (v.z < -tiny_num) v.z += 100;
                    while (v.z > 100 - tiny_num) v.z -= 100;
                }
            }
            
            else
            {
                if (v.x < -tiny_num)
                {
                    v.x *= -1;
                    v.y *= -1;
                    v.z *= -1;
                }
                
                while (v.x > 100 - tiny_num) v.x -= 100;
                
                while (v.y < -tiny_num) v.y += 100;
                while (v.y > 100 - tiny_num) v.y -= 100;
                while (v.z < -tiny_num) v.z += 100;
                while (v.z > 100 - tiny_num) v.z -= 100;
            }
            
            ofs << v.x << " " << v.y << " " << v.z << " " << counter << endl;
        }
        
        Print_Info( counter_1, min_dist );
    }
    
    sort( output_cloud.begin(), output_cloud.end(), Increasing_Distance );
    
    /*ofstream ofs( "/Users/philsmith/Documents/Xcode Projects/PPC_Equivalency_Test/Results/ppc" + to_string( iterations ) + ".txt" );
    
    for (int counter = 0; counter < output_cloud.size(); ++counter)
    {
        ofs << "T2_" + to_string( output_cloud[counter].first ) + "_num_molGeom.cif " << output_cloud[counter].second.x << " " << output_cloud[counter].second.y << endl;
    }*/
    
    Print_Summary( start_time ); // Printing the summary.
    
    ifstream ifs_1( "/Users/philsmith/Documents/Work/Xcode Projects/PPC_Equivalency_Test/Results/0.03.txt" );
    ifstream ifs_2( "/Users/philsmith/Documents/Work/Xcode Projects/PPC_Equivalency_Test/Results/-0.03.txt" );
    
    string line_data;
    vector<Point3d> vec1, vec2;
    
    while (getline( ifs_1, line_data ))
    {
        stringstream stream;
        double x, y, z;
        
        stream << line_data;
        
        stream >> x >> y >> z;
        
        vec1.push_back( Point3d( x, y, z ) );
    }
    
    while (getline( ifs_2, line_data ))
    {
        stringstream stream;
        double x, y, z;
        
        stream << line_data;
        
        stream >> x >> y >> z;
        
        vec2.push_back( Point3d( x, y, z ) );
    }
    
    double max_dist = 0;
    int c = 0;
    
    for (int counter = 0; counter < num_distances; ++counter)
    {
        double dist = abs( vec1[counter].x - vec2[counter].x );
        
        if (dist > max_dist)
        {
            c = counter;
            max_dist = dist;
        }
        
        dist = abs( vec1[counter].y - vec2[counter].y );
        
        if (dist > max_dist)
        {
            c = counter;
            max_dist = dist;
        }
        
        dist = abs( vec1[counter].z - vec2[counter].z );
        
        if (dist > max_dist)
        {
            c = counter;
            max_dist = dist;
        }
    }
    
    cout << "Max dist: " << max_dist << endl;
    cout << c << endl;
    
    return 0;
}
