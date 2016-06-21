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

    // numerical integration
    const double T = 300.0;
    const double beta = 1.0 / (constants::k_B * T);

    const double top = 10.0;
    const double btm = 5.0;
    const std::size_t bins = 10000;
    const double dx = (top - btm) / bins;

    double integration = 0.0;
    for(std::size_t i=0; i < bins; ++i)
    {
        const double x  = btm + i * dx;
        const double f0 =
            std::exp(-1.0 * beta * (perturb.energy(x) + potential.energy(x)));
        const double f1 =
            std::exp(-1.0 * beta * (perturb.energy(x + dx) + potential.energy(x + dx)));

//         integration += dx * (f1 + f0) * 0.5;
        integration += (f1 + f0);
    }
    const double Z = dx * integration * 0.5;

    DCDData data;
//     data.nset() = 100001;
    data.istart() = 0;
    data.nstep_save() = 1;
    data.nstep() = 100000;
    data.nunit() = 1;
    data.verCHARMM() = 24;
    data.nparticle() = 2;
    data.delta_t() = 0.1;
    data.signeture() = "CORD";
    data.push_header("========================== Qahwa sample dcd ");

    const std::size_t total_frames = 100000;
    data.traj().reserve(total_frames);
    for(std::size_t i = 0; i < bins; ++i)
    {
        const double x = btm + i * dx;
        const double probability =
            std::exp(-1.0 * beta * (perturb.energy(x) + potential.energy(x))) / Z;

        const std::size_t frames = static_cast<std::size_t>(probability * dx * total_frames);
        for(std::size_t j = 0; j < frames; ++j)
            data.traj().push_back(SnapShot{Vector3d(0.0, 0.0, 0.0),
                                           Vector3d(rng.uniform(x, x+dx), 0.0, 0.0)});
    }
    data.nset() = data.traj().size();

    DCDWriter writer(output_name);
    writer.data() = data;
    writer.write();

    return EXIT_SUCCESS;
}
