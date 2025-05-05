#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// Generic version of the vector field
void vectorField(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];
  Node d = params[3];

  Node Y_r_1 = in[0];
  Node Y_r_2 = in[1];
  Node Y_i_1 = in[2];
  Node Y_i_2 = in[3];
  Node Z_r_1 = in[4];
  Node Z_r_2 = in[5];
  Node Z_i_1 = in[6];
  Node Z_i_2 = in[7];
  Node lambda_real = in[8];
  Node lambda_imag = in[9];
  Node kappa = in[10];
  Node epsilon = in[11];

  // FIXME: Implement these
  Node N1_a_real = 0 * delta;
  Node N1_a_imag = 0 * delta;
  Node N1_b_real = 0 * delta;
  Node N1_b_imag = 0 * delta;
  Node N2_a_real = 0 * delta;
  Node N2_a_imag = 0 * delta;
  Node N2_b_real = 0 * delta;
  Node N2_b_imag = 0 * delta;

  // The indexing notation _ij corresponds to the Julia indexing convention.

  // M1
  Node M1_11_real = (-epsilon * (kappa / sigma + N1_a_real - lambda_real) - N2_a_real - omega) / (1 + (epsilon^2));
  Node M1_11_imag = (-epsilon * (N1_a_imag - lambda_imag) - N2_a_imag) / (1 + (epsilon^2));

  Node M1_12_real = (-kappa / sigma - epsilon * (N1_b_real - omega) - N2_b_real + lambda_real) / (1 + (epsilon^2));
  Node M1_12_imag = (-epsilon * N1_b_imag - N2_b_imag + lambda_imag) / (1 + (epsilon^2));

  Node M1_21_real = (kappa / sigma - epsilon * (N2_a_real + omega) + N1_a_real - lambda_real) / (1 + (epsilon^2));
  Node M1_21_imag = (epsilon * N2_a_imag + N1_a_imag - lambda_imag) / (1 + (epsilon^2));

  Node M1_22_real = (-epsilon * (kappa / sigma + N2_b_real - lambda_real) + N1_b_real - omega) / (1 + (epsilon^2));
  Node M1_22_imag = (-epsilon * (N2_b_imag - lambda_imag) - N1_b_imag) / (1 + (epsilon^2));

  // M2
  Node M2_11 = kappa / (1 + (epsilon^2)) * (-epsilon) * xi - (d - 1) / xi;
  Node M2_12 = kappa / (1 + (epsilon^2)) * (-1) * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[0] = Z_r_1;
  out[1] = Z_r_2;
  out[2] = Z_i_1;
  out[3] = Z_i_2;

  out[4] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[5] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[6] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[7] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[8] = 0 * lambda_real;
  out[9] = 0 * lambda_imag;
  out[10] = 0 * kappa;
  out[11] = 0 * epsilon;
}

// Vector field handling the case d == 1
void vectorField_d1(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];

  Node Y_r_1 = in[0];
  Node Y_r_2 = in[1];
  Node Y_i_1 = in[2];
  Node Y_i_2 = in[3];
  Node Z_r_1 = in[4];
  Node Z_r_2 = in[5];
  Node Z_i_1 = in[6];
  Node Z_i_2 = in[7];
  Node lambda_real = in[8];
  Node lambda_imag = in[9];
  Node kappa = in[10];
  Node epsilon = in[11];

  // FIXME: Implement these
  Node N1_a_real = 0 * delta;
  Node N1_a_imag = 0 * delta;
  Node N1_b_real = 0 * delta;
  Node N1_b_imag = 0 * delta;
  Node N2_a_real = 0 * delta;
  Node N2_a_imag = 0 * delta;
  Node N2_b_real = 0 * delta;
  Node N2_b_imag = 0 * delta;

  // M1
  // The indices 11, 12, 21 and 22 correspond to Julia indices
  Node M1_11_real = (-epsilon * (kappa / sigma + N1_a_real - lambda_real) - N2_a_real - omega) / (1 + (epsilon^2));
  Node M1_11_imag = (-epsilon * (N1_a_imag - lambda_imag) - N2_a_imag) / (1 + (epsilon^2));

  Node M1_12_real = (-kappa / sigma - epsilon * (N1_b_real - omega) - N2_b_real + lambda_real) / (1 + (epsilon^2));
  Node M1_12_imag = (-epsilon * N1_b_imag - N2_b_imag + lambda_imag) / (1 + (epsilon^2));

  Node M1_21_real = (kappa / sigma - epsilon * (N2_a_real + omega) + N1_a_real - lambda_real) / (1 + (epsilon^2));
  Node M1_21_imag = (epsilon * N2_a_imag + N1_a_imag - lambda_imag) / (1 + (epsilon^2));

  Node M1_22_real = (-epsilon * (kappa / sigma + N2_b_real - lambda_real) + N1_b_real - omega) / (1 + (epsilon^2));
  Node M1_22_imag = (-epsilon * (N2_b_imag - lambda_imag) - N1_b_imag) / (1 + (epsilon^2));

  // M2
  Node M2_11 = kappa / (1 + (epsilon^2)) * (-epsilon) * xi;
  Node M2_12 = kappa / (1 + (epsilon^2)) * (-1) * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[0] = Z_r_1;
  out[1] = Z_r_2;
  out[2] = Z_i_1;
  out[3] = Z_i_2;

  out[4] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[5] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[6] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[7] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[8] = 0 * lambda_real;
  out[9] = 0 * lambda_imag;
  out[10] = 0 * kappa;
  out[11] = 0 * epsilon;
}

int main()
{
  cout.precision(17); // Enough to exactly recover Float64 values
  cerr.precision(17); // Enough to exactly recover Float64 values

  // Read initial value
  // 8 for initial data, 2 for lambda, 2 for kappa and epsilon
  IVector u0(12);
  cin >> u0[0] >> u0[1] >> u0[2] >> u0[3] >> u0[4] >> u0[5] >> u0[6] >> u0[7];

  // Read parameter values
  int d;
  interval lambda_real, lambda_imag, nu_real, nu_imag, kappa, epsilon, omega, sigma, delta;
  cin >> d;
  cin >> lambda_real;
  cin >> lambda_imag;
  cin >> nu_real;
  cin >> nu_imag;
  cin >> kappa;
  cin >> epsilon;
  cin >> omega;
  cin >> sigma;
  cin >> delta;

  u0[8] = lambda_real;
  u0[9] = lambda_imag;
  u0[10] = kappa;
  u0[11] = epsilon;

  // Read time span
  interval T0, T1;
  cin >> T0 >> T1;

  // Read flag for if to output Jacobian
  int output_jacobian;
  cin >> output_jacobian;

  // Read tolerance to use
  double tol;
  cin >> tol;

  // Create the vector field and the parameters
  int dim = 12;
  IMap vf;

  // TODO: Implement optimized versions
   if (d == 1) {
    vf = IMap(vectorField_d1, dim, dim, 3);
    vf.setParameter(0, omega);
    vf.setParameter(1, sigma);
    vf.setParameter(2, delta);
  } else {
    vf = IMap(vectorField, dim, dim, 4);
    vf.setParameter(0, omega);
    vf.setParameter(1, sigma);
    vf.setParameter(2, delta);
    vf.setParameter(3, interval(d));
  }

  // Create the solver and the time map
  IOdeSolver solver(vf, 20);

  solver.setAbsoluteTolerance(tol);
  solver.setRelativeTolerance(tol);

  ITimeMap timeMap(solver);

  try {
    if (output_jacobian) {
      // Define a representation of the initial value
      C1HORect2Set s(u0, T0);

      // Solve the system
      IVector result = timeMap(T1, s);

      IMatrix m = (IMatrix)(s);

      for (int i = 0; i < 8; i++)
          // Only print derivatives of u[0], ..., u[7]
          for (int j = 0; j < 8; j++)
              cout << m[j][i] << endl;

      // Derivative w.r.t. lambda_real
      for (int j = 0; j < 8; j++)
	  cout << m[j][8] << endl;

      // Derivative w.r.t. lambda_imag
      for (int j = 0; j < 8; j++)
	  cout << m[j][9] << endl;
    } else {
      // Define a doubleton representation of the initial value
      C0HORect2Set s(u0, T0);

      // Solve the system
      IVector result = timeMap(T1, s);

      for (int i = 0; i < 8; i++)
        cout << result[i] << endl;
    }
  } catch(exception& e) {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
