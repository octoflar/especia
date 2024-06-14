/// @file runner.h
/// The model runner.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include "runner.h"

especia::Runner::Runner (int argc, char *argv[])
{
  using std::string;

  for (int i = 0; i < argc; ++i)
    {
      args.emplace_back (argv[i]);
    }
}

especia::Runner::~Runner () = default;

void
especia::Runner::write_command_line (std::ostream &os) const
{
  using std::endl;

  os << "<!DOCTYPE html>" << endl;
  os << "<html>" << endl;
  os << "<!--" << endl;
  os << "<command>" << endl;

  for (const auto &arg : args)
    {
      os << " " << arg;
    }

  os << std::endl;
  os << "</command>" << endl;
  os << "-->" << endl;
  os << "</html>" << endl;
}

void
especia::Runner::write_result_messages (std::ostream &os,
                                        const Optimizer::Result &result) const
{
  using std::cout;
  using std::endl;

  os << "<!--" << endl;
  os << "<message>" << endl;

  if (result.is_optimized ())
    {
      os << "especia::Runner::run() Message: optimization completed "
            "successfully"
         << endl;
    }
  else
    {
      os << "especia::Runner::run() Warning: optimization stopped at "
            "generation "
         << result.get_generation_number () << endl;
    }
  if (result.is_underflow ())
    {
      os << "especia::Runner::run() Warning: optimization stopped because of "
            "an "
         << "underflow of the mutation variance" << endl;
    }

  os << "</message>" << endl;
  os << "-->" << endl;
}

void
especia::Runner::write_usage_message (std::ostream &os) const
{
  using std::endl;

  os << project_long_name << " " << project_doi << endl;
  os << "usage: " << get_program_name () << ": "
     << "{seed} {parents} {population} {step} {accuracy} {stop} {trace} < "
        "{model file} [> {result file}]"
     << endl;
}
