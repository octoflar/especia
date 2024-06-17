/// @file exitcodes.h
/// Application error exit codes.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#ifndef ESPECIA_EXITCODES_H
#define ESPECIA_EXITCODES_H

namespace especia
{

/// Application error exit codes.
class Exit_Codes
{
public:
  /// The optimization stopped due to an underflow of the mutation variance
  /// (exit code = 1).
  static const int optimization_underflow = 00001;

  /// The optimization failed to reach the required accuracy goal within the
  /// prescribed number of generations (exit code = 2).
  static const int optimization_stopped = 00002;

  /// A logic error occurred (exit code = 8).
  static const int logic_error = 00010;

  /// A runtime error occurred (exit code = 16).
  static const int runtime_error = 00020;

  /// An unspecific exception occurred (exit code = 64).
  static const int unspecific_exception = 00100;

private:
  Exit_Codes () = default; // private constructor prevents instantiation
};

} // namespace especia

#endif // ESPECIA_EXITCODES_H
