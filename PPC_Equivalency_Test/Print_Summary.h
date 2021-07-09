#pragma once

#include <iostream>
#include <vector>

#include "Atom.h"

using namespace std;

void Print_Info ( int iteration, vector<pair<double, pair<Atom, Atom>>>const& min_dist )
{
    cout << "Crystal " << iteration << endl << endl;
    
    int counter = 1;
    
    for (auto iter = min_dist.rbegin(); iter != min_dist.rend(); ++iter, ++counter)
    {
        cout << "The ";
        cout << counter;
        if (counter % 10 == 1 && counter % 100 != 11) cout << "st";
        else if (counter % 10 == 2 && counter % 100 != 12) cout << "nd";
        else if (counter % 10 == 3 && counter % 100 != 13) cout << "rd";
        else cout << "th";
        cout << " minimum distance is " << iter->first << "." << endl;
        
        Point3d v = iter->second.first.cart_coords - iter->second.second.cart_coords;
        
        //if (v.x < -tiny_num) v *= -1;
        
        /*if (v.x < - tiny_num) v.x += 100;
        if (v.y < - tiny_num) v.y += 100;
        if (v.z < - tiny_num) v.z += 100;*/
        
        if (v.x < - tiny_num) v.x *= -1;
        if (v.y < - tiny_num) v.y *= -1;
        if (v.z < - tiny_num) v.z *= -1;
        
        cout << v << endl;
    }
    
    cout << endl;
}

void Print_Summary ( clock_t start_time )
{
    clock_t end_time = clock();
    
    double duration = (end_time - start_time) * 1000 / (double)CLOCKS_PER_SEC;
    
    cout << "Code duration = " << duration << "ms." << endl << endl;
}
