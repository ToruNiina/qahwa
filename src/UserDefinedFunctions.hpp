#ifndef QAHWA_USER_DEFINITION_FUNCTION
#define QAHWA_USER_DEFINITION_FUNCTION
#include "Definitions.hpp"

namespace qahwa
{

/*!
 *  @brief functor that returns perturbing potential.
 *  This returns the perturbing potential (if you want to use qahwa to umbrella
 *  sampling, this will be umbrella potential).
 *  Normally this will be anchor potential. @n
 *  @n
 *  In the case that you add some abnormal potential to your simulation,
 *  modify this class to be able to calculate your abnormal potential. @n
 *  @n 
 *  It must be guaranteed that the definition of operator() is kept.
 *  operator() accepts SnapShot type and returns double type.
 *
 *  @sa SnapShot
 */
class PerturbingPotential
{
  public:
    PerturbingPotential()  = default;

    PerturbingPotential(const std::size_t index, const double coef, const Vector3d& pos)
        : index_(index), coefficient_(coef), anchored_position_(pos)
    {}

    ~PerturbingPotential() = default;

    double operator()(const SnapShot& snapshot) const;

    double  coef() const {return coefficient_;}
    double& coef()       {return coefficient_;}

    std::size_t  particle_index() const {return index_;}
    std::size_t& particle_index()       {return index_;}

  private:

    std::size_t  index_; //!< index of anchored particle :NOTE: 0-based index!
    double       coefficient_; //!< coefficient of the potential
    Vector3d     anchored_position_; //!< 3d vector the particle is anchored
};


/*!
 *  @brief functor that returns Reaction Coordinate.
 */
class ReactionCoordinate
{
  public:
    ReactionCoordinate() = default;
    ~ReactionCoordinate() = default;

    double operator()(const SnapShot& snapshot) const;
};


}


#endif /* QAHWA_USER_DEFINITION_FUNCTION */
