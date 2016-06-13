#include "UserDefinedFunction.hpp"

namespace qahwa
{

// this is for example. anchor.
double PerturbingPotential::operator()(const SnapShot& snapshot) const
{
    return this->coefficient_ *
           len_square(snapshot[this->index_] - this->anchored_position_);
}

// this is for example. distance from origin.
template<>
double ReactionCoordinate<double>::operator()(const SnapShot& snapshot) const
{
    return length(snapshot[1]);
}


}
