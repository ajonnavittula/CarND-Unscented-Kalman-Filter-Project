#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
  	std::cout << "Invalid estimations size" << std::endl;
  	return rmse;
  }

  for (int i=0; i<estimations.size(); i++)
  {
  	VectorXd residuals = estimations[i] - ground_truth[i];

  	residuals = residuals.array() * residuals.array();
  	rmse += residuals;
  }

  rmse = rmse/estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}