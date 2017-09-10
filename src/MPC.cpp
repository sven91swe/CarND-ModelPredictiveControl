#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "tools.h"

using CppAD::AD;

// Evaluate a polynomial. -- copy from main, should be in a seperate file.
AD<double> polyeval(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}

AD<double> polyeval_derivative(Eigen::VectorXd coeffs, AD<double> x) {
  AD<double> result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * CppAD::pow(x, i-1);
  }
  return result;
}

double polyeval_double(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

double polyeval_derivative_double(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 1; i < coeffs.size(); i++) {
    result += i * coeffs[i] * pow(x, i-1);
  }
  return result;
}

size_t N = 20;
double dt = 0.10;

int numberOfStateVariables = 6; //Position x, position y, speed, heading, cte, e_yaw
int numberOfActuatorVariables = 2; //Acceleration and steering

int x_start = 0;
int y_start = x_start + N;
int v_start = y_start + N;
int yaw_start = v_start + N;
int cte_start = yaw_start + N;
int eyaw_start = cte_start + N;


int acc_start = eyaw_start + N;
int steer_start = acc_start + N - 1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0.0;



    //Cost function - States
    for(int t=0; t<N; t++){
      fg[0] += 10*CppAD::pow(vars[v_start + t] - 50 * MPC::mph_to_meters_per_second, 2); //Speed
      fg[0] += 2000*CppAD::pow(vars[cte_start + t], 2); //Distance to center
      fg[0] += 500*CppAD::pow(vars[eyaw_start + t], 2); //Angle
    }

    //Cost function - actuator values
    for(int t=0; t<(N - 1); t++){
      fg[0] += 5 * CppAD::pow(vars[acc_start + t] - 1, 2); //Low acceleration
      fg[0] += 50 * CppAD::pow(vars[steer_start + t], 2); //Low steering values
    }

    //Cost function - actuator differences
    for(int t=1; t<(N - 1); t++){
      fg[0] += 1 * CppAD::pow(vars[acc_start + t] - vars[acc_start + t - 1], 2); //Smoother acceleration
      fg[0] += 10 * CppAD::pow(vars[steer_start + t] - vars[steer_start + t - 1], 2); //Smoother steering
    }

    //Constraints
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[v_start + 1] = vars[v_start];
    fg[yaw_start + 1] = vars[yaw_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[eyaw_start + 1] = vars[eyaw_start];

    for(int t=1; t<N; t++){
      //State at time t + 1
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> yaw1 = vars[yaw_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> eyaw1 = vars[eyaw_start + t];

      //State at time t
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> yaw0 = vars[yaw_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> eyaw0 = vars[eyaw_start + t - 1];

      //Acctuators
      AD<double> acc0 = vars[acc_start + t - 1];
      AD<double> steering0  = vars[steer_start + t - 1];

      //Polynomial evaluations
      AD<double> f0 = polyeval(coeffs, x0);
      AD<double> f0_derivative = polyeval_derivative(coeffs, x0);
      AD<double> yaw_line = CppAD::atan(coeffs[1]);

      //Constraints between states
      fg[x_start + t + 1] = x1 - (x0 + v0 * CppAD::cos(yaw0) * dt);
      fg[y_start + t + 1] = y1 - (y0 + v0 * CppAD::sin(yaw0) * dt);
      fg[v_start + t + 1] = v1 - (v0 + acc0 * dt);
      fg[yaw_start + t + 1] = yaw1 - (yaw0 + v0 / Lf * steering0 * dt);
      fg[cte_start + t + 1] = cte1 - (f0 - y0  + v0 * CppAD::sin(eyaw0) * dt);
      fg[eyaw_start + t +1] = eyaw1 - (yaw0 - yaw_line + (v0 / Lf * steering0 * dt));
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x0 = state[0];
  double y0 = state[1];
  double v0 = state[2];
  double yaw0 = state[3];

  //Single update before optimization due to latency
  double x = x0 + v0 * cos(yaw0) * 0.1;
  double y = y0 + v0 * sin(yaw0) * 0.1;
  double v = v0 + lastAccelerator * 0.1;
  double yaw = yaw0 + v0 / Lf * lastSteering * 0.1;
  double cte = polyeval_double(coeffs, x) - y;
  double eyaw = atan(polyeval_derivative_double(coeffs, x)) - yaw;

  size_t n_vars = N*(numberOfStateVariables + numberOfActuatorVariables) - numberOfActuatorVariables;

  size_t n_constraints = N*numberOfStateVariables;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++){
    vars[i] = 0;
  }

  vars[x_start] = x;
  vars[y_start] = y;
  vars[v_start] = v;
  vars[yaw_start] = yaw;
  vars[cte_start] = cte;
  vars[eyaw_start] = eyaw;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for(int i = x_start; i < (x_start + N); i++){
    vars_lowerbound[i] = -100;
    vars_upperbound[i] = 100;
  }
  for(int i = y_start; i < (y_start + N); i++){
    vars_lowerbound[i] = -100;
    vars_upperbound[i] = 100;
  }
  for(int i = v_start; i < (v_start + N); i++){
    vars_lowerbound[i] = -10;
    vars_upperbound[i] = 100;
  }
  for(int i = yaw_start; i < (yaw_start + N); i++){
    vars_lowerbound[i] = -2.0*M_PI;
    vars_upperbound[i] = 2.0*M_PI;
  }
  for(int i = cte_start; i < (cte_start + N); i++){
    vars_lowerbound[i] = -100;
    vars_upperbound[i] = 100;
  }
  for(int i = eyaw_start; i < (eyaw_start + N); i++){
    vars_lowerbound[i] = -2.0*M_PI;
    vars_upperbound[i] = 2.0*M_PI;
  }

  for(int i = acc_start; i < (acc_start + N - 1); i++){
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  for(int i = steer_start; i < (steer_start + N - 1); i++){
    vars_lowerbound[i] = -25.0/180 * M_PI;
    vars_upperbound[i] = 25.0/180 * M_PI;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_upperbound[x_start] = x;

  constraints_lowerbound[y_start] = y;
  constraints_upperbound[y_start] = y;

  constraints_lowerbound[v_start] = v;
  constraints_upperbound[v_start] = v;

  constraints_lowerbound[yaw_start] = yaw;
  constraints_upperbound[yaw_start] = yaw;

  constraints_lowerbound[cte_start] = cte;
  constraints_upperbound[cte_start] = cte;

  constraints_lowerbound[eyaw_start] = eyaw;
  constraints_upperbound[eyaw_start] = eyaw;


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
  options += "Numeric max_cpu_time          5.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;


  mpc_x = vector<double>(N);
  mpc_y = vector<double>(N);

  for(int t=0; t<N; t++){
    mpc_x.at(t) = solution.x[x_start + t];
    mpc_y.at(t) = solution.x[y_start + t];
  }

  double steeringToReturn = solution.x[steer_start + 0] + solution.x[steer_start + 1] + solution.x[steer_start + 2];
  steeringToReturn /= 3.0;


  double accelerationToReturn = solution.x[acc_start + 0] + solution.x[acc_start + 1] + solution.x[acc_start + 2];
  accelerationToReturn /= 3.0;

  lastSteering = steeringToReturn;
  lastAccelerator = accelerationToReturn;

  return {steeringToReturn, accelerationToReturn};
}
