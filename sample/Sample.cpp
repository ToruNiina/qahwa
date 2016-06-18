#include "Sample.hpp"
using namespace qahwa;

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <anchor distance> <anchor coef>" << std::endl;
        return EXIT_FAILURE;
    }

    const double anchor_dist = std::stod(std::string(argv[1]));
    const double anchor_coef = std::stod(std::string(argv[2]));
    const std::string output_name("sample_dist_" + std::string(argv[1]) +
                                  "_coef_" + std::string(argv[2]) + ".dcd");
    std::cout << "distance = "    << anchor_dist << std::endl;
    std::cout << "coefficient = " << anchor_coef << std::endl;
    std::cout << "output = " << output_name << std::endl;

    const PerturbingPotential perturb(anchor_coef, anchor_dist);
    const TruePotential potential;

    std::random_device rnd;
    const unsigned int seed = rnd();
    std::cout << "seed = " << seed << std::endl;
    RandomNumberGenerator rng(seed);

    Metropolis metropolis(/*T = */ 300.0);
    const double top = 10.0;
    const double btm = 5.0;

    // the state of the system
    System system;
    system.pos = anchor_dist;
    system.ene = perturb.energy(system.pos) + potential.energy(system.pos);

    DCDData data;
    data.nset() = 100001;
    data.istart() = 0;
    data.nstep_save() = 1;
    data.nstep() = 100000;
    data.nunit() = 1;
    data.verCHARMM() = 24;
    data.nparticle() = 2;
    data.delta_t() = 0.1;
    data.signeture() = "CORD";
    data.push_header("========================== Qahwa sample dcd ");
    data.traj().push_back(make_snapshot(system));

    for(std::size_t tstep = 0; tstep < 100000; ++tstep)
    {
        while(true)
        {
            const double trial_pos = rng.uniform(btm, top);
            const double trial_ene = perturb.energy(trial_pos) + 
                                     potential.energy(trial_pos);
            const double dE = trial_ene - system.ene;
            if(metropolis.accept(dE, rng))
            {
                system.pos = trial_pos;
                system.ene = trial_ene;
                data.traj().push_back(make_snapshot(system));
                break;
            }
        }
    }
    DCDWriter writer(output_name);
    writer.data() = data;
    writer.write();

    return EXIT_SUCCESS;
}
