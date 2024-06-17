/// @file optimizer.h
/// CMA-ES classes for nonlinear function optimization.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#ifndef ESPECIA_OPTIMIZER_H
#define ESPECIA_OPTIMIZER_H

#include <cstddef>
#include <functional>

#include "decompose.h"
#include "deviates.h"
#include "optimize.h"
#include "random.h"

using std::valarray;

namespace especia
{

/// No constraint.
///
/// @tparam T The number type.
template <class T = real> class No_Constraint
{
public:
  /// The constructor.
  No_Constraint () = default;

  /// The destructor.
  ~No_Constraint () = default;

  /// Tests if a given parameter vector violates the constraint.
  ///
  /// @param[in] x The parameter vector.
  /// @param[in] n The number of parameters to test.
  /// @return always @c false.
  bool
  is_violated (const T x[], natural n) const
  {
    return false;
  }

  /// Computes the cost associated with the constraint.
  ///
  /// @param[in] x The parameter vector.
  /// @param[in] n The number of parameters to take account of.
  /// @return always zero.
  T
  cost (const T x[], natural n) const
  {
    return T (0);
  }
};

/// No tracing.
///
/// @tparam T The number type.
template <class T = real> class No_Tracing
{
public:
  /// The constructor.
  No_Tracing () = default;

  /// The destructor.
  ~No_Tracing () = default;

  /// Tests if tracing is enabled.
  ///
  /// @param[in] g The generation number.
  /// @return always @c false.
  bool
  is_tracing (natural g) const
  {
    return false;
  }

  /// Traces state information.
  ///
  /// @param[in] g The generation number.
  /// @param[in] y The value of the objective function.
  /// @param[in] min_step The minimum mutation step size.
  /// @param[in] max_step The maximum mutation step size.
  void
  trace (natural g, T y, T min_step, T max_step) const
  {
  }
};

/// An optimizer based on the CMA-ES developed by Hansen and Ostermeier (2001).
///
/// Further reading:
///
/// N. Hansen, S. D. Müller, P. Koumoutsakos (2003).
///  *Reducing the Increasing the Time Complexity of the Derandomized Evolution
///     Strategy with Covariance Matrix Adaption (CMA-ES).*
///  Evolutionary Computation, 11, 1, ISSN 1063-6560.
///
/// N. Hansen, A. Ostermeier (2001).
///   *Completely Derandomized Self-Adaption in Evolution Strategies.*
///   Evolutionary Computation, 9, 159, ISSN 1063-6560.
class Optimizer
{
public:
  /// Builds a new optimizer.
  class Builder
  {
  public:
    /// Default constructor.
    Builder ();

    /// The destructor.
    ~Builder ();

    /// Builds a new optimizer.
    ///
    /// @return the optimizer.
    Optimizer build ();

    /// Returns the problem dimension.
    ///
    /// @return the problem dimension.
    natural
    get_problem_dimension () const
    {
      return n;
    }

    /// Returns the parent number.
    ///
    /// @return the parent number.
    natural
    get_parent_number () const
    {
      return parent_number;
    }

    /// Returns the population size.
    ///
    /// @return the population size.
    natural
    get_population_size () const
    {
      return population_size;
    }

    /// Returns the covariance matrix update modulus.
    ///
    /// @return the covariance matrix update modulus.
    natural
    get_covariance_update_modulus () const
    {
      return update_modulus;
    }

    /// Returns the accuracy goal.
    ///
    /// @return the accuracy goal.
    real
    get_accuracy_goal () const
    {
      return accuracy_goal;
    }

    /// Returns the random seed.
    ///
    /// @return the random seed.
    word64
    get_random_seed () const
    {
      return random_seed;
    }

    /// Returns the stop generation.
    ///
    /// @return the stop generation.
    natural
    get_stop_generation () const
    {
      return stop_generation;
    }

    /// Returns the recombination weights.
    ///
    /// @return the recombination weights.
    const std::valarray<real> &
    get_weights () const
    {
      return weights;
    }

    /// Returns the step size cumulation rate.
    ///
    /// @return the step size cumulation rate.
    real
    get_step_size_cumulation_rate () const
    {
      return cs;
    }

    /// Returns the distribution cumulation rate.
    ///
    /// @return the distribution cumulation rate.
    real
    get_distribution_cumulation_rate () const
    {
      return cc;
    }

    /// Returns the rank-1 covariance matrix adaption rate.
    ///
    /// @return the rank-1 covariance matrix adaption rate.
    real
    get_rank_1_covariance_matrix_adaption_rate () const
    {
      return ccov;
    }

    /// Returns the rank-µ covariance matrix adaption rate.
    ///
    /// @return the rank-µ covariance matrix adaption rate.
    real
    get_rank_m_covariance_matrix_adaption_rate () const
    {
      return acov;
    }

    /// Returns the step size damping.
    ///
    /// @return the step size damping.
    real
    get_step_size_damping () const
    {
      return step_size_damping;
    }

    /// Configures default settings.
    ///
    /// @return the problem dimension.
    Builder &with_defaults ();

    /// Configures the problem dimension.
    ///
    /// @param[in] n The problem dimension.
    /// @return this builder.
    Builder &with_problem_dimension (natural n = 10);

    /// Configures the parent number and sets the population size to twice the
    /// parent number.
    ///
    /// @param[in] parent_number The parent number.
    /// @return this builder.
    Builder &with_parent_number (natural parent_number = 20);

    /// Configures the population size.
    ///
    /// @attention The population size must be greater than or equal to twice
    /// the parent number.
    ///
    /// @param[in] population_size The population size.
    /// @return this builder.
    Builder &with_population_size (natural population_size);

    /// Configures the covariance matrix update modulus.
    ///
    /// @param[in] update_modulus The update modulus.
    /// @return this builder.
    Builder &with_covariance_update_modulus (natural update_modulus = 1);

    /// Configures the accuracy goal.
    ///
    /// @param[in] accuracy_goal The accuracy goal.
    /// @return this builder.
    Builder &with_accuracy_goal (real accuracy_goal = 1.0E-06);

    /// Configures the random seed.
    ///
    /// @param[in] seed The random seed.
    /// @return this builder.
    Builder &with_random_seed (word64 seed = 9600629759793949339ULL);

    /// Configures the stop generation.
    ///
    /// @param[in] stop_generation The stop generation.
    /// @return this builder.
    Builder &with_stop_generation (natural stop_generation = 1000);

  private:
    /// Returns a pointer to the recombination weights.
    ///
    /// @return a pointer to the recombination weights.
    const real *
    get_weights_pointer () const
    {
      return &weights[0];
    }

    /// Configures strategy parameters like recombination weights, cumulation
    /// and adaption rates according to Hansen (2014,
    /// http://cma.gforge.inria.fr/purecmaes.m).
    void with_strategy_parameters ();

    /// The problem dimension.
    natural n = 10;

    /// The parent number.
    natural parent_number = 20;

    /// The population size.
    natural population_size = 40;

    /// The covariance matrix update modulus.
    natural update_modulus = 1;

    /// The accuracy goal.
    real accuracy_goal = 1.0E-6;

    /// The random seed.
    word64 random_seed = 271828;

    /// The stop generation.
    natural stop_generation = 1000;

    /// The recombination weights.
    std::valarray<real> weights;

    /// The variance of the recombination weights.
    real wv;

    /// The step size cumulation rate.
    real cs;

    /// The distribution cumulation rate.
    real cc;

    /// The rank-µ covariance matrix adaption rate.
    real acov;

    /// The rank-1 covariance matrix adaption rate.
    real ccov;

    /// The step size damping.
    real step_size_damping;

    friend class Optimizer;
  };

  /// The optimization result.
  class Result
  {
  public:
    /// The destructor.
    ~Result ();

    /// Returns the covariance matrix (upper triangular part only, in
    /// column-major layout).
    ///
    /// @return the covariance matrix.
    const std::valarray<real> &
    get_covariance_matrix () const
    {
      return C;
    }

    /// Returns the distribution cumulation path.
    ///
    /// @return the distribution cumulation path.
    const std::valarray<real> &
    get_distribution_cumulation_path () const
    {
      return pc;
    }

    /// Returns the optimized fitness.
    ///
    /// @return the optimized fitness.
    real
    get_fitness () const
    {
      return y;
    }

    /// Returns the final generation number.
    ///
    /// @return the final generation number.
    natural
    get_generation_number () const
    {
      return g;
    }

    /// Returns the final global step size.
    ///
    /// @return the final global step size.
    real
    get_global_step_size () const
    {
      return s;
    }

    /// Returns the final local step sizes.
    ///
    /// @return the final local step sizes.
    const std::valarray<real> &
    get_local_step_sizes () const
    {
      return d;
    }

    /// Returns the optimized parameter values.
    ///
    /// @return the optimized parameter values.
    const std::valarray<real> &
    get_parameter_values () const
    {
      return x;
    }

    /// Returns the parameter uncertainties.
    ///
    /// @return the parameter uncertainties.
    const std::valarray<real> &
    get_parameter_uncertainties () const
    {
      return z;
    }

    /// Returns the final rotation matrix (in column-major layout).
    ///
    /// @return the final rotation matrix.
    const std::valarray<real> &
    get_rotation_matrix () const
    {
      return B;
    }

    /// Returns the step size cumulation path.
    ///
    /// @return the step size cumulation path.
    const std::valarray<real> &
    get_step_size_cumulation_path () const
    {
      return ps;
    }

    /// Returns the optimization status flag.
    ///
    /// @return the optimization status flag.
    bool
    is_optimized () const
    {
      return optimized;
    }

    /// Returns the mutation variance underflow status flag.
    ///
    /// @return the mutation variance underflow status flag.
    bool
    is_underflow () const
    {
      return underflow;
    }

  private:
    /// The constructor.
    ///
    /// @param[in] n The problem dimension.
    /// @param[in] x The initial parameter values.
    /// @param[in] d The initial local step sizes.
    /// @param[in] s The initial global step size.
    Result (natural n, const std::valarray<real> &x,
            const std::valarray<real> &d, real s);

    /// Returns a pointer to the covariance matrix.
    ///
    /// @return a pointer to the covariance matrix.
    real *
    get_covariance_matrix_pointer ()
    {
      return &C[0];
    }

    /// Returns a pointer to the distribution cumulation path.
    ///
    /// @return a pointer to the distribution cumulation path.
    real *
    get_distribution_cumulation_path_pointer ()
    {
      return &pc[0];
    }

    /// Returns a reference to the fitness.
    ///
    /// @return a reference to the fitness.
    real &
    __fitness ()
    {
      return y;
    }

    /// Returns a reference to the generation number.
    ///
    /// @return a reference to the generation number.
    natural &
    __generation_number ()
    {
      return g;
    }

    /// Returns a reference to the global step size.
    ///
    /// @return a reference to the global step size.
    real &
    __global_step_size ()
    {
      return s;
    }

    /// Returns a pointer to the local step sizes.
    ///
    /// @return a pointer to the local step sizes.
    real *
    get_local_step_sizes_pointer ()
    {
      return &d[0];
    }

    /// Returns a pointer to the parameter values.
    ///
    /// @return a pointer to the parameter values.
    real *
    get_parameter_values_pointer ()
    {
      return &x[0];
    }

    /// Returns a pointer to the parameter uncertainties.
    ///
    /// @return a pointer to the parameter uncertainties.
    real *
    get_parameter_uncertainties_pointer ()
    {
      return &z[0];
    }

    /// Returns a pointer to the rotation matrix.
    ///
    /// @return a pointer to the rotation matrix.
    real *
    get_rotation_matrix_pointer ()
    {
      return &B[0];
    }

    /// Returns a pointer to the step size cumulation path.
    ///
    /// @return a pointer to the step size cumulation path.
    real *
    get_step_size_cumulation_path_pointer ()
    {
      return &ps[0];
    }

    /// Returns a reference to the optimization status flag.
    ///
    /// @return a reference to the optimization status flag.
    bool &
    __optimized ()
    {
      return optimized;
    }

    /// Returns a reference to the mutation variance underflow status flag.
    ///
    /// @return a reference to the mutation variance underflow status flag.
    bool &
    __underflow ()
    {
      return underflow;
    }

    /// The optimized parameter values.
    std::valarray<real> x;

    /// The final local step sizes.
    std::valarray<real> d;

    /// The final global step size.
    real s;

    /// The parameter uncertainties.
    std::valarray<real> z;

    /// The optimized fitness.
    real y;

    /// The final rotation matrix (in column-major layout).
    std::valarray<real> B;

    /// The final covariance matrix (upper triangular part only, in
    /// column-major layout).
    std::valarray<real> C;

    /// The distribution cumulation path.
    std::valarray<real> pc;

    /// The step size cumulation path.
    std::valarray<real> ps;

    /// The optimization status flag.
    bool optimized;

    /// The mutation variance underflow status flag.
    bool underflow;

    /// The final generation number.
    natural g;

    friend class Optimizer;
  };

  /// The destructor.
  ~Optimizer ();

  /// Maximizes an objective function.
  ///
  /// @tparam F The function type.
  /// @tparam Constraint The constraint type.
  /// @tparam Tracing The tracer type.
  ///
  /// @param[in] f The objective function.
  /// @param[in] x The initial parameter values.
  /// @param[in] d The initial local step sizes.
  /// @param[in] s The initial global step size.
  /// @param[in] constraint The constraint.
  /// @param[in] tracer The tracer.
  ///
  /// @return the maximization result.
  template <class F, class Constraint, class Tracing>
  Result
  maximize (const F &f, const std::valarray<real> &x,
            const std::valarray<real> &d, const real &s,
            const Constraint &constraint, const Tracing &tracer) const
  {
    return optimize (f, x, d, s, constraint, tracer, std::greater<real> ());
  }

  /// @overload
  template <class F>
  Result
  maximize (const F &f, const std::valarray<real> &x,
            const std::valarray<real> &d, const real &s) const
  {
    return optimize (f, x, d, s, No_Constraint<real> (), No_Tracing<real> (),
                     std::greater<real> ());
  }

  /// Minimizes an objective function.
  ///
  /// @tparam F The function type.
  /// @tparam Constraint The constraint type.
  /// @tparam Tracing The tracer type.
  ///
  /// @param[in] f The objective function.
  /// @param[in] x The initial parameter values.
  /// @param[in] d The initial local step sizes.
  /// @param[in] s The initial global step size.
  /// @param[in] constraint The constraint.
  /// @param[in] tracer The tracer.
  ///
  /// @return the minimization result.
  template <class F, class Constraint, class Tracing>
  Result
  minimize (const F &f, const std::valarray<real> &x,
            const std::valarray<real> &d, const real &s,
            const Constraint &constraint, const Tracing &tracer) const
  {
    return optimize (f, x, d, s, constraint, tracer, std::less<real> ());
  }

  /// @overload
  template <class F>
  Result
  minimize (const F &f, const std::valarray<real> &x,
            const std::valarray<real> &d, const real &s) const
  {
    return optimize (f, x, d, s, No_Constraint<real> (), No_Tracing<real> (),
                     std::less<real> ());
  }

private:
  /// Creates a new instance of this class with the build configuration
  /// supplied as argument.
  ///
  /// @param[in] builder The build configuration.
  explicit Optimizer (const Builder &builder);

  /// Optimizes an objective function.
  ///
  /// @tparam F The function type.
  /// @tparam Constraint The constraint type.
  /// @tparam Tracing The tracer type.
  /// @tparam Compare The fitness comparator type.
  ///
  /// @param[in] f The objective function.
  /// @param[in] x The initial parameter values.
  /// @param[in] d The initial local step sizes.
  /// @param[in] s The initial global step size.
  /// @param[in] constraint The constraint.
  /// @param[in] tracer The tracer.
  /// @param[in] compare The fitness comparator.
  ///
  /// @return the optimization result.
  template <class F, class Constraint, class Tracing, class Compare>
  Result
  optimize (const F &f, const std::valarray<real> &x,
            const std::valarray<real> &d, const real &s,
            const Constraint &constraint, const Tracing &tracer,
            const Compare &compare) const
  {
    using especia::optimize;
    using especia::postopti;

    const natural n = config.get_problem_dimension ();

    Result result (n, x, d, s);

    optimize (
        f, constraint, n, config.get_parent_number (),
        config.get_population_size (), config.get_weights_pointer (),
        config.get_step_size_damping (),
        config.get_step_size_cumulation_rate (),
        config.get_distribution_cumulation_rate (),
        config.get_rank_1_covariance_matrix_adaption_rate (),
        config.get_rank_m_covariance_matrix_adaption_rate (),
        config.get_covariance_update_modulus (), config.get_accuracy_goal (),
        config.get_stop_generation (), result.__generation_number (),
        result.get_parameter_values_pointer (), result.__global_step_size (),
        result.get_local_step_sizes_pointer (),
        result.get_rotation_matrix_pointer (),
        result.get_covariance_matrix_pointer (),
        result.get_step_size_cumulation_path_pointer (),
        result.get_distribution_cumulation_path_pointer (),
        result.__fitness (), result.__optimized (), result.__underflow (),
        deviate, decompose, compare, tracer);

    if (result.__optimized ())
      {
        postopti (f, constraint, n, result.get_parameter_values_pointer (),
                  result.get_local_step_sizes_pointer (),
                  result.get_rotation_matrix_pointer (),
                  result.get_covariance_matrix_pointer (),
                  result.get_global_step_size (),
                  result.get_parameter_uncertainties_pointer ());
      }

    return result;
  }

  /// The build configuration.
  const Builder config;

  /// The eigenvalue decomposition strategy.
  const Decompose decompose;

  /// The random number generator.
  const Normal_Deviate<Melg19937_64> deviate;
};

} // namespace especia

#endif // ESPECIA_OPTIMIZER_H
