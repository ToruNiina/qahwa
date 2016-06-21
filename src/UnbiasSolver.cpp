#include "UnbiasSolver.hpp"

namespace qahwa
{

ProbabilityDensityFunction UnbiasSolver::unbias(const Trajectory& traj,
                                      const PerturbingPotential perturb,
                                      const ReactionCoordinate& rctcrd) const 
{
    HistogramMaker histmaker(this->bins_);
    const auto range = histmaker.find_range(traj, rctcrd);
    histmaker.set_range(range);

    ProbabilityDensityFunction unbiased(this->bins_, range.first, range.second);
//         histmaker.make_prob_dens_func(traj, rctcrd);

    const double min_p = 1.0 / traj.size();
    for(auto iter = traj.cbegin(); iter != traj.cend(); ++iter)
    {
        const double x = rctcrd(*iter);
        const double W = perturb(*iter);
        unbiased.at(x) += min_p * std::exp(this->beta_ * W);
    }

    return unbiased;
}

ProbabilityDensityFunction
UnbiasSolver::reconstruct(const std::vector<window_type>& windows,
                          const ReactionCoordinate& rctcrd)
{
    double begin = std::numeric_limits<double>::max();
    double end   = std::numeric_limits<double>::min();
    for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
    {
        for(auto snap = iter->first.cbegin(); snap != iter->first.cend(); ++snap)
        {
            const double x = rctcrd(*snap);
            if(begin > x) begin = x;
            if(end   < x) end   = x;
        }
    }
    end += (end - begin) / this->bins_ * 0.01;

    ProbabilityDensityFunction reconst(this->bins_, begin, end);

    return reconst;
}


Histogram<double>
UnbiasSolver::make_pmf(const ProbabilityDensityFunction& pdf) const
{
    Histogram<double> pmf(pdf.bins(), pdf.range_begin(), pdf.range_end());
    auto potential = pmf.begin();
    for(auto iter = pdf.cbegin(); iter != pdf.cend(); ++iter)
    {
        *potential = -1.0 * std::log(*iter) / this->beta_;
        ++potential;
    }
    return pmf;
}




}
