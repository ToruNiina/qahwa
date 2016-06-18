#include "UserDefinedFunctions.hpp"
#include "HistogramMaker.hpp"
#include "UnbiasSolver.hpp"
#include "WhamSolver.hpp"
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
                auto pdf = make_pdf(histogram);

                std::ofstream ofs(input.get_as<std::string>(
                                      input.at("histogram","output")) +
                                  std::to_string(index) + ".dat");
                if(!ofs.good())
                    throw std::runtime_error("file open error : qahwa_histogram.dat");
                ofs << pdf;
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
            InputFileReader input(cmdline.argv_at(2));
            input.read();

            std::vector<std::string> fnames = 
                input.get_as_list<std::string>(input.at("wham","dcdfiles"));
            std::vector<std::string> parameter_list =
                input.split_list(input.at("wham", "parameters"));

            double temperature;
            std::size_t bins;
            try
            {
                temperature =
                    input.get_as<double>(input.at("wham","temperature"));
                std::cerr << "temperature = " << temperature << std::endl;
            }
            catch(std::exception& exp)
            {
                std::cerr << "temperature is set as default. 300K" << std::endl;
                temperature = 300.0;
            }
            try
            {
                bins =
                    input.get_as<std::size_t>(input.at("wham","bins"));
                std::cerr << "bins  = " << bins << std::endl;
            }
            catch(std::exception& exp)
            {
                std::cerr << "the number of bins is set as default. 100"
                          << std::endl;
                bins = 100;
            }

            if(fnames.size() != parameter_list.size())
                throw std::runtime_error("files and parameters are different size");

            WhamSolver solver(bins, temperature);

            ReactionCoordinate rctcrd;
            std::vector<std::pair<Trajectory, PerturbingPotential>> windows;
            windows.reserve(fnames.size());
            for(std::size_t i=0; i<fnames.size(); ++i)
            {
                DCDReader dcdreader(fnames.at(i));
                dcdreader.read();
                auto trajectory = dcdreader.data().traj();

                std::vector<double> parameters =
                    input.get_as_list<double>(parameter_list.at(i));
                const std::size_t anch_index = 
                    static_cast<std::size_t>(parameters.at(0));
                const Vector3d anch = Vector3d(parameters.at(2),
                                               parameters.at(3),
                                               parameters.at(4));

                PerturbingPotential perturb(anch_index, parameters.at(1), anch);
                windows.push_back(std::make_pair(std::move(trajectory), perturb));
            }
 
            std::cerr << "trajectory reading end. start solving" << std::endl;
            auto parameter = solver.solve(windows);

            std::cout << "parameter solved" << std::endl;
            for(auto iter = parameter.cbegin(); iter != parameter.cend(); ++iter)
            {
                std::cout << *iter << std::endl;
            }

            auto unbiased = solver.reconstruct(parameter, windows, rctcrd);

            std::ofstream pdf(input.get_as<std::string>(
                                  input.at("wham","output")) +
                              "_pdf.dat");
            if(!pdf.good())
                throw std::runtime_error("file open error : " + 
                        input.get_as<std::string>(input.at("wham","output")) +
                                         "_pdf.dat");
            pdf << unbiased;
            pdf.close();

            auto potential = solver.make_pmf(unbiased);
            std::ofstream pmf(input.get_as<std::string>(
                                  input.at("wham","output")) +
                              "_pmf.dat");
            if(!pmf.good())
                throw std::runtime_error("file open error : " +
                        input.get_as<std::string>(input.at("wham","output")) +
                                         "_pmf.dat");
            pmf << potential;
            pmf.close();

            break;
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
