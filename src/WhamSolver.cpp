#include "WhamSolver.hpp"

namespace qahwa
{

std::vector<double> WhamSolver::solve(
        const std::vector<window_type>& windows) const
{
    std::vector<double> expfs_prev(windows.size(), 0.5);

    while(true)
    {
        std::vector<double> expfs(windows.size(), 0.0);
        auto expf = expfs.begin();
        for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
        {
            for(auto win = windows.cbegin(); win != windows.cend(); ++win)
            {
                const PerturbingPotential& pot = win->second;
                for(auto step = win->first.cbegin(); step != win->first.cend(); ++step)
                {
                    const double numer = std::exp(-1e0 * beta_ * pot(*step));
                    const double denom = calc_denominator(windows, expfs_prev, *step);
                    *expf += numer / denom;
                }
            }
            ++expf;
        }
        if(is_converged(expfs_prev, expfs))
            return expfs;
        else
            expfs_prev.swap(expfs);
    }
    throw std::logic_error("never reach here");
}

double WhamSolver::calc_denominator(const std::vector<window_type>& windows,
                                    const std::vector<double>& expfs,
                                    const SnapShot& snapshot) const
{
    assert(expfs.size() == windows.size());

    const double minus_beta = -1.0 * beta_;
    double denom = 0.0;
    auto expf = expfs.cbegin();
    for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
    {
        denom += (iter->first).size() *
                 std::exp(minus_beta * (iter->second(snapshot))) * (*expf);
        ++expf;
    }
    return denom;
}


bool WhamSolver::is_converged(const std::vector<double>& f_prev,
                              const std::vector<double>& f) const
{
    assert(f_prev.size() == f.size());
    auto iter_prev = f_prev.cbegin();
    for(auto iter = f.cbegin(); iter != f.cend(); ++iter_prev, ++iter)
    {
        if(std::abs(*iter - *iter_prev) > tolerance) return false;
    }
    return true;
}

ProbabilityDensityFunction
WhamSolver::reconstruct(const std::vector<double>& free_energy_parameter,
                        const std::vector<window_type>& windows,
                        const ReactionCoordinate& rctcrd) const
{
    assert(windows.size() == free_energy_parameter.size());

    // find total range
    double begin = std::numeric_limits<double>::max();
    double end   = std::numeric_limits<double>::min();
    for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
    {
        for(auto step = iter->first.cbegin(); step != iter->first.cend(); ++step)
        {
            const double x = rctcrd(*step);
            if(x < begin) begin = x;
            if(end < x)   end   = x;
        }
    }
    end += (end - begin) * 1e-3;

    // reconst unbiased pdf(R) and project it into reaction coordinate
    ProbabilityDensityFunction unbiased(this->bins_, begin, end);
    for(auto win = windows.cbegin(); win != windows.cend(); ++win)
    {
        for(auto step = win->first.cbegin(); step != win->first.cend(); ++step)
        {
            double denom = 0.0;
            auto expf = free_energy_parameter.cbegin();
            for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
            {
                denom += iter->first.size() *
                         std::exp(-1.0 * this->beta_ * iter->second(*step)) *
                         (*expf);
                ++expf;
            }
            unbiased.at(rctcrd(*step)) += 1.0 / denom;
        }
    }

    return unbiased;
}

}
