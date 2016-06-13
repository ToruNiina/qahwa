#include "HistogramMaker.hpp"
#include <limits>

namespace qahwa
{

Histogram<std::size_t> HistogramMaker::make_histogram(
        const Trajectory& traj,
        const ReactionCoordinate<double>& reaction_coord) const
{
    histogram_type hist(this->bins_, this->begin_, this->end_);

    for(auto iter = traj.cbegin(); iter != traj.cend(); ++iter)
        ++(hist.at(reaction_coord(*iter)));

    return hist;
}

std::pair<double, double> HistogramMaker::find_range(
        const Trajectory& traj,
        const ReactionCoordinate<double>& reaction_coord) const
{
    double range_begin = std::numeric_limits<double>::max();
    double range_end   = std::numeric_limits<double>::min();

    for(auto iter = traj.cbegin(); iter != traj.cend(); ++iter)
    {
        const double crd = reaction_coord(*iter);
        if(crd < range_begin)
            range_begin = crd;
        if(range_end < crd)
            range_end = crd;
    }

    return std::make_pair(range_begin, range_end);
}

}
