#include "UserDefinedFunctions.hpp"

namespace qahwa
{

// this is for example. anchor.
double PerturbingPotential::operator()(const SnapShot& snapshot) const
{
    return this->coefficient_ *
           len_square(snapshot.at(this->index_) - this->anchored_position_);
}

// this is for example. distance from origin.
double ReactionCoordinate::operator()(const SnapShot& snapshot) const
{
    return length(snapshot.at(0) - snapshot.at(1));
}


}
