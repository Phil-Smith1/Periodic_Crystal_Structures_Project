// Returns the index of side of the cell with the minimum cell length.

#pragma once

int Min_Cell_Length ( double length_1, double length_2, double length_3 )
{
    int min = (length_1 < length_2) ? 0 : 1;
    
    if (min == 0)
    {
        min = (length_1 < length_3) ? min : 2;
    }
    
    else
    {
        min = (length_2 < length_3) ? min : 2;
    }
    
    return min;
}
