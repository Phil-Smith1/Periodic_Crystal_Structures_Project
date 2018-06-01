#pragma once

#include <iostream>

void Print_Summary ( clock_t start_time )
{
    clock_t end_time = clock();
    
    double duration = (end_time - start_time) * 1000 / (double)CLOCKS_PER_SEC;
    
    cout << "Code duration = " << duration << "ms." << endl << endl;
}
