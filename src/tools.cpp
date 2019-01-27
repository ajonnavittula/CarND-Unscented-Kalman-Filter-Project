#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
  	cout << "Invalid estimations size" << endl;
  	return rmse;
  }

  for (unsigned int i = 0; i < estimations.size(); i++)
  {
  	VectorXd residuals = estimations[i] - ground_truth[i];

  	residuals = residuals.array() * residuals.array();
  	rmse += residuals;
  }

  rmse = rmse/estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}