#pragma once

#include "Atom.h"

#ifndef Tiny_Num
#define Tiny_Num 1
const double tiny_num = 1e-10;
#endif

bool Sort_Min_Distances ( pair<double, pair<Atom, Atom>>const& i, pair<double, pair<Atom, Atom>>const& j )
{
    return i.first > j.first;
}

bool Sort_Min_Distances_2 ( pair<double, pair<Atom, Atom>>const& i, pair<double, pair<Atom, Atom>>const& j )
{
    Point3d v1 = i.second.first.cart_coords - i.second.second.cart_coords;
    Point3d v2 = j.second.first.cart_coords - j.second.second.cart_coords;
    
    if (abs( v1.x ) < tiny_num)
    {
        if (abs( v1.y ) < tiny_num)
        {
            if (v1.z < -tiny_num) v1.z *= -1;
            
            while (v1.z > 100 - tiny_num) v1.z -= 100;
        }
        
        else
        {
            if (v1.y < -tiny_num)
            {
                v1.y *= -1;
                v1.z *= -1;
            }
            
            while (v1.y > 100 - tiny_num) v1.y -= 100;
            
            while (v1.z < -tiny_num) v1.z += 100;
            while (v1.z > 100 - tiny_num) v1.z -= 100;
        }
    }
    
    else
    {
        if (v1.x < -tiny_num)
        {
            v1.x *= -1;
            v1.y *= -1;
            v1.z *= -1;
        }
        
        while (v1.x > 100 - tiny_num) v1.x -= 100;
        
        while (v1.y < -tiny_num) v1.y += 100;
        while (v1.y > 100 - tiny_num) v1.y -= 100;
        while (v1.z < -tiny_num) v1.z += 100;
        while (v1.z > 100 - tiny_num) v1.z -= 100;
    }
    
    if (abs( v2.x ) < tiny_num)
    {
        if (abs( v2.y ) < tiny_num)
        {
            if (v2.z < -tiny_num) v2.z *= -1;
            
            while (v2.z > 100 - tiny_num) v2.z -= 100;
        }
        
        else
        {
            if (v2.y < -tiny_num)
            {
                v2.y *= -1;
                v2.z *= -1;
            }
            
            while (v2.y > 100 - tiny_num) v2.y -= 100;
            
            while (v2.z < -tiny_num) v2.z += 100;
            while (v2.z > 100 - tiny_num) v2.z -= 100;
        }
    }
    
    else
    {
        if (v2.x < -tiny_num)
        {
            v2.x *= -1;
            v2.y *= -1;
            v2.z *= -1;
        }
        
        while (v2.x > 100 - tiny_num) v2.x -= 100;
        
        while (v2.y < -tiny_num) v2.y += 100;
        while (v2.y > 100 - tiny_num) v2.y -= 100;
        while (v2.z < -tiny_num) v2.z += 100;
        while (v2.z > 100 - tiny_num) v2.z -= 100;
    }
    
    if (abs( v1.x - v2.x ) < tiny_num)
    {
        if (v1.y < v2.y - tiny_num) return true;
        if (v2.y < v1.y - tiny_num) return false;
        
        else if (v1.z < v2.z - tiny_num) return true;
        else return false;
    }
    
    else if (v1.x < v2.x - tiny_num) return true;
    else return false;
}
