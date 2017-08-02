#include "MPC.h"
#include "Constants.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;

namespace
{
  class FG_eval final {
  public:
    // Fitted polynomial coefficients
    Eigen::Ref<const Eigen::VectorXd> coeffs;
    FG_eval(const Eigen::VectorXd& coeffs) : coeffs(coeffs) {}
    
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& vars) {
      fg[0] = 0;

      // The part of the cost based on the reference state.
      for (size_t t = 0; t < Constants::N; ++t) {
        fg[0] += Constants::wCte * CppAD::pow(vars[Constants::idxCte + t], 2);
        fg[0] += Constants::wEpsi * CppAD::pow(vars[Constants::idxEpsi + t], 2);
        fg[0] += Constants::wV * CppAD::pow(vars[Constants::idxV + t] - Constants::refV, 2);
      }


      // Minimize the value gap between sequential errors (cte and epsi).
      for (size_t t = 0; t < Constants::N - 1; ++t) {
        fg[0] += Constants::wDCte * CppAD::pow(vars[Constants::idxCte + t + 1] - vars[Constants::idxCte + t], 2);
        fg[0] += Constants::wDEpsi * CppAD::pow(vars[Constants::idxEpsi + t + 1] - vars[Constants::idxEpsi + t], 2);
      }

      // Minimize the use of actuators.
      for (size_t t = 0; t < Constants::N - 1; ++t) {
        fg[0] += Constants::wSteering * CppAD::pow(vars[Constants::idxSteering + t], 2);
        fg[0] += Constants::wAcc * CppAD::pow(vars[Constants::idxAcc + t], 2);
      }

      // Minimize the value gap between sequential actuations.
      for (size_t t = 0; t < Constants::N - 2; ++t) {
        fg[0] += Constants::wDsteering * CppAD::pow(vars[Constants::idxSteering + t + 1] - vars[Constants::idxSteering + t], 2);
        fg[0] += Constants::wDacc * CppAD::pow(vars[Constants::idxAcc + t + 1] - vars[Constants::idxAcc + t], 2);
      }

      //
      // Setup Constraints
      //
      // NOTE: In this section you'll setup the model constraints.

      // Initial constraints
      //
      // We add 1 to each of the starting indices due to cost being located at
      // index 0 of `fg`.
      // This bumps up the position of all the other values.
      fg[1 + Constants::idxX] = vars[Constants::idxX];
      fg[1 + Constants::idxY] = vars[Constants::idxY];
      fg[1 + Constants::idxPsi] = vars[Constants::idxPsi];
      fg[1 + Constants::idxV] = vars[Constants::idxV];
      fg[1 + Constants::idxCte] = vars[Constants::idxCte];
      fg[1 + Constants::idxEpsi] = vars[Constants::idxEpsi];
      
      // The rest of the constraints
      for (size_t t = 1; t < Constants::N; ++t) {
        // The state at time t+1 .
        const AD<double> x1 = vars[Constants::idxX + t];
        const AD<double> y1 = vars[Constants::idxY + t];
        const AD<double> psi1 = vars[Constants::idxPsi + t];
        const AD<double> v1 = vars[Constants::idxV + t];
        const AD<double> cte1 = vars[Constants::idxCte + t];
        const AD<double> epsi1 = vars[Constants::idxEpsi + t];

        // The state at time t.
        const AD<double> x0 = vars[Constants::idxX + t - 1];
        const AD<double> y0 = vars[Constants::idxY + t - 1];
        const AD<double> psi0 = vars[Constants::idxPsi + t - 1];
        const AD<double> v0 = vars[Constants::idxV + t - 1];
        const AD<double> cte0 = vars[Constants::idxCte + t - 1];
        const AD<double> epsi0 = vars[Constants::idxEpsi + t - 1];

        // Only consider the actuation at time t.
        const AD<double> delta0 = -vars[Constants::idxSteering + t - 1];
        const AD<double> a0 = vars[Constants::idxAcc + t - 1];

        AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0*x0 + coeffs[3] * x0*x0*x0;
        AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0*x0);

        // Here's `x` to get you started.
        // The idea here is to constraint this value to be 0.
        //
        // Recall the equations for the model:
        // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
        // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
        // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
        // v_[t+1] = v[t] + a[t] * dt
        // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
        // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
        fg[1 + Constants::idxX + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * Constants::dt);
        fg[1 + Constants::idxY + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * Constants::dt);
        fg[1 + Constants::idxPsi + t] = psi1 - (psi0 + v0 * delta0 / Constants::Lf * Constants::dt);
        fg[1 + Constants::idxV + t] = v1 - (v0 + a0 * Constants::dt);
        fg[1 + Constants::idxCte + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * Constants::dt));
        fg[1 + Constants::idxEpsi + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Constants::Lf * Constants::dt);
      }
    }

    FG_eval(const FG_eval&) = delete;
    FG_eval(FG_eval&&) = delete;
    FG_eval& operator=(const FG_eval&) = delete;
    FG_eval& operator=(FG_eval&&) = delete;
  };
}

std::vector<double> MPC::solve(const Eigen::VectorXd& state, const Eigen::VectorXd& coeffs,
  std::vector<double>& mpcX, std::vector<double>& mpcY) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  const double x = state[0];
  const double y = state[1];
  const double psi = state[2];
  const double v = state[3];
  const double cte = state[4];
  const double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(Constants::numVars);
  std::fill(vars.data(), vars.data() + Constants::numVars, 0);

  // Set the initial variable values
  vars[Constants::idxX] = x;
  vars[Constants::idxY] = y;
  vars[Constants::idxPsi] = psi;
  vars[Constants::idxV] = v;
  vars[Constants::idxCte] = cte;
  vars[Constants::idxEpsi] = epsi;

  Dvector vars_lowerbound(Constants::numVars);
  Dvector vars_upperbound(Constants::numVars);

  // TODO: Set lower and upper limits for variables.
  std::fill(vars_lowerbound.data(), vars_lowerbound.data() + Constants::idxSteering, -1e7);
  std::fill(vars_upperbound.data(), vars_upperbound.data() + Constants::idxSteering, 1e7);

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  std::fill(vars_lowerbound.data() + Constants::idxSteering, vars_lowerbound.data() + Constants::idxAcc, -0.436332);
  std::fill(vars_upperbound.data() + Constants::idxSteering, vars_upperbound.data() + Constants::idxAcc, 0.436332);

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  std::fill(vars_lowerbound.data() + Constants::idxAcc, vars_lowerbound.data() + Constants::numVars, -1.0);
  std::fill(vars_upperbound.data() + Constants::idxAcc, vars_upperbound.data() + Constants::numVars, 1.0);

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(Constants::numConstraints);
  Dvector constraints_upperbound(Constants::numConstraints);
  std::fill(constraints_lowerbound.data(), constraints_lowerbound.data() + Constants::numConstraints, 0);
  std::fill(constraints_upperbound.data(), constraints_upperbound.data() + Constants::numConstraints, 0);

  constraints_lowerbound[Constants::idxX] = x;
  constraints_lowerbound[Constants::idxY] = y;
  constraints_lowerbound[Constants::idxPsi] = psi;
  constraints_lowerbound[Constants::idxV] = v;
  constraints_lowerbound[Constants::idxCte] = cte;
  constraints_lowerbound[Constants::idxEpsi] = epsi;

  constraints_upperbound[Constants::idxX] = x;
  constraints_upperbound[Constants::idxY] = y;
  constraints_upperbound[Constants::idxPsi] = psi;
  constraints_upperbound[Constants::idxV] = v;
  constraints_upperbound[Constants::idxCte] = cte;
  constraints_upperbound[Constants::idxEpsi] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
    options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
    constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  mpcX.clear();
  mpcY.clear();

  for (size_t i = 0; Constants::N > i; ++i) {
    mpcX.push_back(solution.x[Constants::idxX + i]);
    mpcY.push_back(solution.x[Constants::idxY + i]);
  }

  return { solution.x[Constants::idxSteering], solution.x[Constants::idxAcc] };
}
