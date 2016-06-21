#include "WhamSolver.hpp"
#include <iomanip>

namespace qahwa
{

/* solve this for F self-consistently.
 *
 * win  = windows.size();
 * step = windows.at(i).first.size();
 * W is perturbing potential
 * F is exp(beta f)
 *
 * for l = {0 to win}
 *    1                  /               /                 exp[-beta * W_l(R_ij)]             \ \
 *  ----- = sum_{i}^{win}| sum_{j}^{step}| -------------------------------------------------- | |
 *   F_l                 \               \  sum_{k}^{win}(n_k * exp[-beta * W_k(R_ij)] * F_k) / /
 *
 *  the value of exp[-beta * W_l(R_ij)] is constant through the iteration, 
 *  so cash_expW stores this value to solve fast.
 *  cash_expW[l][i][j] = exp[-beta * W_l(R_ij)];
 *  ==> cash_expW[win][win][step]
 */

std::vector<double> WhamSolver::solve(
        const std::vector<window_type>& windows) const
{
    const std::size_t win_size   = windows.size();
    const double      minus_beta = -1.0 * this->beta_;
    std::vector<double> expfs_prev(win_size, 1.0);

    // make cash
    std::vector<std::vector<std::vector<double>>> cash_expW(
            win_size, std::vector<std::vector<double>>(win_size));

    for(std::size_t l = 0; l < win_size; ++l)
    {
        const PerturbingPotential& W_l = windows.at(l).second;
        for(std::size_t i = 0; i < win_size; ++i)
        {
            (cash_expW.at(l).at(i)).resize(windows.at(i).first.size());
            for(std::size_t j = 0; j < windows.at(i).first.size(); ++j)
            {
                cash_expW.at(l).at(i).at(j) = 
                    std::exp(minus_beta * W_l(windows.at(i).first.at(j)));
            }
        }
    }
    // end

    std::size_t iteration_times = 0;
    while(true)
    {
        std::vector<double> expfs(windows.size(), 0.0);
        auto expf = expfs.begin();
        for(std::size_t l = 0; l < win_size; ++l)
        {
            for(std::size_t i = 0; i < win_size; ++i)
            {
                for(std::size_t j = 0; j < windows.at(i).first.size(); ++j)
                {
                    const double numer = cash_expW.at(l).at(i).at(j);
                    double denom = 0.0;
                    for(std::size_t k = 0; k < win_size; ++k)
                    {
                        denom += windows.at(k).first.size() * 
                                 cash_expW.at(k).at(i).at(j) *
                                 expfs_prev.at(k);
                    }
                    *expf += numer / denom;
                }
            }
            ++expf;
        }
        for(auto iter = expfs.begin(); iter != expfs.end(); ++iter)
            *iter = 1.0 / (*iter);

        ++iteration_times;
        std::cerr << "now " << iteration_times << "-th iteration" << std::endl;
        for(auto iter = expfs.cbegin(); iter != expfs.cend(); ++iter)
            std::cout << std::setprecision(6) << *iter << " ";
        std::cout << std::endl;

        if(is_converged(expfs_prev, expfs))
        {
            return expfs;
        }
        else
        {
            expfs_prev = expfs;
        }

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
        if(std::abs((*iter - *iter_prev) / (*iter)) > tolerance) return false;
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

    // equation 7
//     ProbabilityDensityFunction unbiased(this->bins_, begin, end);
//
//     HistogramMaker maker(this->bins_, begin, end);
//     std::vector<ProbabilityDensityFunction> biased;
//     biased.reserve(windows.size());
//     for(auto iter = windows.cbegin(); iter != windows.cend(); ++iter)
//     {
//         biased.push_back(maker.make_prob_dens_func(iter->first, rctcrd));
//     }
//
//     const double dx = (end - begin) / this->bins_;
//     for(std::size_t index = 0; index < this->bins_; ++index)
//     {
//         const double x = begin + dx * (index + 0.5);
//
//         double numer = 0.0;
//         for(std::size_t j = 0; j < windows.size(); ++j)
//         {
//             numer += windows.at(j).first.size() *
//                 std::exp(-1.0 * this->beta_ * windows.at(j).second(x)) *
//                 free_energy_parameter.at(j);
//         }
//
//         double denom = 0.0;
//         for(std::size_t i = 0; i < windows.size(); ++i)
//         {
//             denom += windows.at(i).first.size() * biased.at(i).at(x);
//         }
//
//         unbiased.at(x) = denom / numer;
//     }

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
