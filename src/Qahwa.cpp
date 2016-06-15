#include "UserDefinedFunctions.hpp"
#include "HistogramMaker.hpp"
#include "UnbiasSolver.hpp"
#include "CmdParser.hpp"
#include "InputFileReader.hpp"
#include <fstream>

using namespace qahwa;

int main(int argc, char *argv[])
{
    try{

    CommandLineParser cmdline(argc, argv);
    cmdline.parse();

    switch(cmdline.mode())
    {
        case Mode::HISTOGRAM_DCD:
        {
            DCDReader dcdreader(cmdline.argv_at(2));
            dcdreader.read();
            auto trajectory = dcdreader.data().traj();
            HistogramMaker histmaker(100);
            ReactionCoordinate rctcrd;
            histmaker.set_range(histmaker.find_range(trajectory, rctcrd));
            auto histogram = histmaker.make_histogram(trajectory, rctcrd);
            std::ofstream ofs("qahwa_histogram.dat");
            if(!ofs.good())
                throw std::runtime_error("file open error : qahwa_histogram.dat");
            ofs << histogram;
            break;
        }
        case Mode::HISTOGRAM:
        {
            InputFileReader input(cmdline.argv_at(2));
            input.read();
            std::vector<std::string> fnames = 
                input.get_as_list<std::string>(input.at("histogram","dcdfiles"));

            std::size_t bins;
            try
            {
                bins =
                    input.get_as<std::size_t>(input.at("unbiased","bins"));
            }
            catch(std::exception& exp)
            {
                std::cerr << "the number of bins is set as default. 100"
                          << std::endl;
                bins = 100;
            }

            std::size_t index = 0;
            for(auto filename = fnames.cbegin(); filename != fnames.cend(); ++filename)
            {
                ++index;

                DCDReader dcdreader(*filename);
                dcdreader.read();
                auto trajectory = dcdreader.data().traj();
 
                HistogramMaker histmaker(bins);
                ReactionCoordinate rctcrd;
                histmaker.set_range(histmaker.find_range(trajectory, rctcrd));
                auto histogram = histmaker.make_histogram(trajectory, rctcrd);               

                std::ofstream ofs(input.get_as<std::string>(
                                      input.at("histogram","output")) +
                                  std::to_string(index) + ".dat");
                if(!ofs.good())
                    throw std::runtime_error("file open error : qahwa_histogram.dat");
                ofs << histogram;
            }
            break;
        }
        case Mode::UNBIASED:
        {
            InputFileReader input(cmdline.argv_at(2));
            input.read();

            std::vector<std::string> fnames = 
                input.get_as_list<std::string>(input.at("unbiased","dcdfiles"));
            std::vector<std::string> parameter_list =
                input.split_list(input.at("unbiased", "parameters"));

            double temperature;
            std::size_t bins;
            try
            {
                temperature =
                    input.get_as<double>(input.at("unbiased","temperature"));
            }
            catch(std::exception& exp)
            {
                std::cerr << "temperature is set as default. 300K" << std::endl;
                temperature = 300.0;
            }
            try
            {
                bins =
                    input.get_as<std::size_t>(input.at("unbiased","bins"));
            }
            catch(std::exception& exp)
            {
                std::cerr << "the number of bins is set as default. 100"
                          << std::endl;
                bins = 100;
            }

            if(fnames.size() != parameter_list.size())
                throw std::runtime_error("files and parameters are different size");

            UnbiasSolver solver(bins, temperature);

            for(std::size_t i=0; i<fnames.size(); ++i)
            {
                DCDReader dcdreader(fnames.at(i));
                dcdreader.read();
                auto trajectory = dcdreader.data().traj();

                ReactionCoordinate rctcrd;
                std::vector<double> parameters =
                    input.get_as_list<double>(parameter_list.at(i));
                const std::size_t anch_index = 1;
//                     static_cast<std::size_t>(parameters.at(0));
                const Vector3d anch = Vector3d(parameters.at(2),
                                               parameters.at(3),
                                               parameters.at(4));

                PerturbingPotential perturb(anch_index, parameters.at(1), anch);
 
                auto unbiased = solver.unbias(trajectory, perturb, rctcrd);

                std::ofstream ofs(input.get_as<std::string>(
                                      input.at("unbiased","output")) +
                                  std::to_string(i) + ".dat");
                if(!ofs.good())
                    throw std::runtime_error("file open error : qahwa_histogram.dat");
                ofs << unbiased;
            }
            break;
        }
        case Mode::WHAM:
        {
            throw std::runtime_error("not implemented yet");
        }  
        default:
            throw std::logic_error("invalid enum");
    }

    }
    catch(std::exception& excpt)
    {
        std::cerr << "Error: " << excpt.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
