#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
public:
  vector<double> mpc_x;
  vector<double> mpc_y;
  double lastSteering = 0;
  double lastAccelerator = 0;

  static constexpr const double mph_to_meters_per_second = 0.44704;

  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
