#ifndef QAHWA_COMMAND_LINE_PARSER
#define QAHWA_COMMAND_LINE_PARSER
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace qahwa
{

enum class Mode
{
    HISTOGRAM_DCD,
    HISTOGRAM,
    UNBIASED,
    WHAM,
};

class CommandLineParser
{
  public:
    CommandLineParser(int argc, char **argv)
        : argc_(argc)
    {
        if(argc < 3)
        {
            std::cerr << this->help() << std::endl;
            throw std::invalid_argument("invalid command");
        }

        for(int i=0; i<argc; ++i)
        {
            argv_.push_back(std::string(argv[i]));
        }
    }
    ~CommandLineParser() = default;

    Mode mode() const {return mode_;}
          std::string& argv_at(std::size_t i)       {return argv_.at(i);}
    const std::string& argv_at(std::size_t i) const {return argv_.at(i);}

    void parse()
    {
        if(argv_.at(1) == "--make-histogram")
        {
            if(argv_.at(2).substr(argv_.at(2).size() - 3, 3) == "dcd")
                mode_ = Mode::HISTOGRAM_DCD;
            else
                mode_ = Mode::HISTOGRAM;
        }
        else if(argv_.at(1) == "--make-unbiased")
        {
            mode_ = Mode::UNBIASED;
        }
        else if(argv_.at(1) == "--wham")
        {
            mode_ = Mode::WHAM;
        }
        else
            throw std::invalid_argument("unknown mode: " + argv_.at(1));
    }

    std::string help() const
    {
        std::ostringstream oss;
        oss << "Usage: " << "qahwa --mode input.file" << std::endl;
        oss << "for example, " << std::endl;
        oss << "--make-histogram traj.dcd" << std::endl;
        oss << "--make-histogram input.toml" << std::endl;
        oss << "--make-unbiased  input.toml" << std::endl;
        oss << "--wham           input.toml" << std::endl;
        oss << std::endl << "see README." << std::endl;
        return oss.str();
    }

  private:
    Mode mode_;
    int argc_;
    std::vector<std::string> argv_;
};



}

#endif /* COMMAND_LINE_PARSER */
