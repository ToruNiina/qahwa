#ifndef QAHWA_SAMPLE
#define QAHWA_SAMPLE
#include "DCDWriter.hpp"
#include <random>
#include <cstddef>

namespace qahwa
{

using DCDData = coffeemill::DCDData;
using DCDWriter = coffeemill::DCDWriter;
using Vector3d = ax::Vector3d;
using SnapShot = std::vector<Vector3d>;
using Trajectory = std::vector<SnapShot>;

namespace constants
{
static constexpr double      k_B     = 1.986231313e-3;// [kcal/mol/K]
}


class PerturbingPotential
{
  public:
    PerturbingPotential(const double coef, const double anch_dist)
        : coef_(coef), anch_(anch_dist)
    {}

    double energy(const double pos) const
    {
//         return 0;
        return coef_ * (pos - anch_) * (pos - anch_);
    }

  private:

    const double coef_;
    const double anch_;
};

// to find (1/x potential, 1/x^2 force)
class TruePotential
{
  public:

    TruePotential()
        : ke(10.0), q1(-1.0), q2(1.0)
    {}

    double energy(const double pos) const
    {
        return q1 * q2 * ke / pos;
    }

  private:

    const double ke;
    const double q1;
    const double q2;
};

class RandomNumberGenerator
{
  public:
    RandomNumberGenerator(unsigned int seed = 10)
        : seed_(seed), mt_(seed)
    {}

    const double uniform(const double min, const double max)
    {
        return (std::uniform_real_distribution<double>(min, max))(this->mt_);
    }

    const unsigned int get_seed() const {return seed_;}

  private:
    const unsigned int seed_;
    std::mt19937       mt_;
};

class Metropolis
{
  public:

    Metropolis(const double temperature)
    {
        this->beta_ = 1.0 / (temperature * constants::k_B);
    }

    bool accept(const double dE, RandomNumberGenerator& rng)
    {
        if(dE <= 0.0)
            return true;
        else
            return (rng.uniform(0.0, 1.0) < std::exp(-1.0 * this->beta_ * dE));
    }

    double  beta() const {return beta_;}
    double& beta()       {return beta_;}

  private:
    double beta_;

};

struct System
{
    double pos;
    double ene;
};

SnapShot make_snapshot(const System& sys)
{
    SnapShot ss;
    ss.push_back(Vector3d(0e0, 0e0, 0e0));
    ss.push_back(Vector3d(sys.pos, 0e0, 0e0));
    return ss;
}


}

#endif /* QAHWA_SAMPLE */
