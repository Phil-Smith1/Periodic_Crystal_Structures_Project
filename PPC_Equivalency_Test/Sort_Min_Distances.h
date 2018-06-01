#pragma once

#include "Atom.h"

bool Sort_Min_Distances ( pair<double, pair<Atom, Atom>>const& i, pair<double, pair<Atom, Atom>>const& j )
{
    return i.first > j.first;
}
