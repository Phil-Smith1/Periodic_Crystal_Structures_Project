// Declaring the class Atom.

#pragma once

#include <opencv2/highgui/highgui.hpp>

using namespace cv;

class Atom
{
public:
    
    int index;
    Point3d frac_coords;
    Point3d cart_coords;
    
    Atom ( int i, Point3d f_c );
    
    Atom();
    ~Atom();
};
