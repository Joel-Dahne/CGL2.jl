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

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];

  Node Y_r_1 = in[4];
  Node Y_r_2 = in[5];
  Node Y_i_1 = in[6];
  Node Y_i_2 = in[7];
  Node Z_r_1 = in[8];
  Node Z_r_2 = in[9];
  Node Z_i_1 = in[10];
  Node Z_i_2 = in[11];

  Node lambda_real = in[12];
  Node lambda_imag = in[13];
  Node kappa = in[14];
  Node epsilon = in[15];

  // Compute forward ODE
  Node mkappa = -kappa;
  Node momega = -omega;

  Node a2b2_sigma = exp(sigma * log((a^2) + (b^2)));

  Node F1 = -(d - 1) / xi * (alpha + epsilon * beta) +
      mkappa * xi * beta +
      mkappa / sigma * b +
      momega * a -
      a2b2_sigma * a +
      delta * a2b2_sigma * b;

  Node F2 = -(d - 1) / xi * (beta - epsilon * alpha) -
      mkappa * xi * alpha -
      mkappa / sigma * a +
      momega * b -
      a2b2_sigma * b -
      delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;

  // Compute linearized ODE

  // Compute J_N
  Node a2b2_sigmam1 = exp((sigma - 1) * log((a^2) + (b^2)));
  Node N1_a = -a2b2_sigmam1 * (delta * (1 + 2 * sigma) * (a^2) + 2 * sigma * a * b + delta * (b^2));
  Node N1_b = -a2b2_sigmam1 * ((a^2) + 2 * delta * sigma * a * b + (1 + 2 * sigma) * (b^2));
  Node N2_a = a2b2_sigmam1 * ((1 + 2 * sigma) * (a^2) - 2 * delta * sigma * a * b + (b^2));
  Node N2_b = -a2b2_sigmam1 * (delta * (a^2) - 2 * sigma * a * b + delta * (1 + 2 * sigma) * (b^2));

  // The indexing notation _ij corresponds to the Julia indexing convention.

  // M1
  Node M1_11_real = (-epsilon * (kappa / sigma + N1_a - lambda_real) - N2_a - omega) / (1 + (epsilon^2));
  Node M1_11_imag = epsilon * lambda_imag / (1 + (epsilon^2));

  Node M1_12_real = (-kappa / sigma - epsilon * (N1_b - omega) - N2_b + lambda_real) / (1 + (epsilon^2));
  Node M1_12_imag = lambda_imag / (1 + (epsilon^2));

  Node M1_21_real = (kappa / sigma - epsilon * (N2_a + omega) + N1_a - lambda_real) / (1 + (epsilon^2));
  Node M1_21_imag = -lambda_imag / (1 + (epsilon^2));

  Node M1_22_real = (-epsilon * (kappa / sigma + N2_b - lambda_real) + N1_b - omega) / (1 + (epsilon^2));
  Node M1_22_imag = epsilon * lambda_imag / (1 + (epsilon^2));

  // M2
  Node M2_11 = kappa / (1 + (epsilon^2)) * (-epsilon) * xi - (d - 1) / xi;
  Node M2_12 = kappa / (1 + (epsilon^2)) * (-1) * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[4] = Z_r_1;
  out[5] = Z_r_2;
  out[6] = Z_i_1;
  out[7] = Z_i_2;

  out[8] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[9] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[10] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[11] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[12] = 0 * lambda_real;
  out[13] = 0 * lambda_imag;
  out[14] = 0 * kappa;
  out[15] = 0 * epsilon;
}

// Vector field handling the case d == 1
void vectorField_d1(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];

  Node Y_r_1 = in[4];
  Node Y_r_2 = in[5];
  Node Y_i_1 = in[6];
  Node Y_i_2 = in[7];
  Node Z_r_1 = in[8];
  Node Z_r_2 = in[9];
  Node Z_i_1 = in[10];
  Node Z_i_2 = in[11];
  Node lambda_real = in[12];
  Node lambda_imag = in[13];
  Node kappa = in[14];
  Node epsilon = in[15];

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

  // TODO: Implement forward ODE
  out[0] = 0 * in[0];
  out[1] = 0 * in[1];
  out[2] = 0 * in[2];
  out[3] = 0 * in[3];

  out[4] = Z_r_1;
  out[5] = Z_r_2;
  out[6] = Z_i_1;
  out[7] = Z_i_2;

  out[8] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[9] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[10] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[11] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[12] = 0 * lambda_real;
  out[13] = 0 * lambda_imag;
  out[14] = 0 * kappa;
  out[15] = 0 * epsilon;
}

// Vector field optimized for sigma = 1
void vectorField_sigma1(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];
  Node d = params[3];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];

  Node Y_r_1 = in[4];
  Node Y_r_2 = in[5];
  Node Y_i_1 = in[6];
  Node Y_i_2 = in[7];
  Node Z_r_1 = in[8];
  Node Z_r_2 = in[9];
  Node Z_i_1 = in[10];
  Node Z_i_2 = in[11];

  Node lambda_real = in[12];
  Node lambda_imag = in[13];
  Node kappa = in[14];
  Node epsilon = in[15];

  // Compute forward ODE
  Node mkappa = -kappa;
  Node momega = -omega;

  Node a2b2_sigma = (a^2) + (b^2);

  Node F1 = -(d - 1) / xi * (alpha + epsilon * beta) +
      mkappa * xi * beta +
      mkappa / sigma * b +
      momega * a -
      a2b2_sigma * a +
      delta * a2b2_sigma * b;

  Node F2 = -(d - 1) / xi * (beta - epsilon * alpha) -
      mkappa * xi * alpha -
      mkappa / sigma * a +
      momega * b -
      a2b2_sigma * b -
      delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;

  // Compute linearized ODE

  // Compute J_N
  Node a2b2_sigmam1 = 1 + 0 * a2b2_sigma;
  Node N1_a = -a2b2_sigmam1 * (delta * (1 + 2 * sigma) * (a^2) + 2 * sigma * a * b + delta * (b^2));
  Node N1_b = -a2b2_sigmam1 * ((a^2) + 2 * delta * sigma * a * b + (1 + 2 * sigma) * (b^2));
  Node N2_a = a2b2_sigmam1 * ((1 + 2 * sigma) * (a^2) - 2 * delta * sigma * a * b + (b^2));
  Node N2_b = -a2b2_sigmam1 * (delta * (a^2) - 2 * sigma * a * b + delta * (1 + 2 * sigma) * (b^2));

  // The indexing notation _ij corresponds to the Julia indexing convention.

  // M1
  Node M1_11_real = (-epsilon * (kappa / sigma + N1_a - lambda_real) - N2_a - omega) / (1 + (epsilon^2));
  Node M1_11_imag = (-epsilon * (-lambda_imag)) / (1 + (epsilon^2));

  Node M1_12_real = (-kappa / sigma - epsilon * (N1_b - omega) - N2_b + lambda_real) / (1 + (epsilon^2));
  Node M1_12_imag = (lambda_imag) / (1 + (epsilon^2));

  Node M1_21_real = (kappa / sigma - epsilon * (N2_a + omega) + N1_a - lambda_real) / (1 + (epsilon^2));
  Node M1_21_imag = (-lambda_imag) / (1 + (epsilon^2));

  Node M1_22_real = (-epsilon * (kappa / sigma + N2_b - lambda_real) + N1_b - omega) / (1 + (epsilon^2));
  Node M1_22_imag = (-epsilon * (-lambda_imag)) / (1 + (epsilon^2));

  // M2
  Node M2_11 = kappa / (1 + (epsilon^2)) * (-epsilon) * xi - (d - 1) / xi;
  Node M2_12 = kappa / (1 + (epsilon^2)) * (-1) * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[4] = Z_r_1;
  out[5] = Z_r_2;
  out[6] = Z_i_1;
  out[7] = Z_i_2;

  out[8] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[9] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[10] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[11] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[12] = 0 * lambda_real;
  out[13] = 0 * lambda_imag;
  out[14] = 0 * kappa;
  out[15] = 0 * epsilon;
}

// Vector field specialized for d = 3, sigma = 1, delta = 0, omega = 1
void vectorField_optimized_d3(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params[]*/, int /*noParams*/)
{
  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];

  Node Y_r_1 = in[4];
  Node Y_r_2 = in[5];
  Node Y_i_1 = in[6];
  Node Y_i_2 = in[7];
  Node Z_r_1 = in[8];
  Node Z_r_2 = in[9];
  Node Z_i_1 = in[10];
  Node Z_i_2 = in[11];

  Node lambda_real = in[12];
  Node lambda_imag = in[13];
  Node kappa = in[14];
  Node epsilon = in[15];

  Node a2 = (a^2);
  Node b2 = (b^2);
  Node a2b2 = a2 + b2;

  // Compute forward ODE
  Node F1 = -2 / xi * (alpha + epsilon * beta) -
      kappa * xi * beta -
      kappa * b -
      a -
      a2b2 * a;

  Node F2 = -2 / xi * (beta - epsilon * alpha) +
      kappa * xi * alpha +
      kappa * a -
      b -
      a2b2 * b;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;

  // Compute linearized ODE

  // Compute J_N
  Node N1_a = -2 * a * b;
  Node N1_b = -(a2 + 3 * b2);
  Node N2_a = 3 * a2 + b2;
  Node N2_b = -N1_a;

  // The indexing notation _ij corresponds to the Julia indexing convention.

  // M1
  Node M1_11_real = (-epsilon * (kappa + N1_a - lambda_real) - N2_a - 1) / one_p_epsilon2;
  Node M1_11_imag = (epsilon * lambda_imag) / one_p_epsilon2;

  Node M1_12_real = (-kappa - epsilon * (N1_b - 1) - N2_b + lambda_real) / one_p_epsilon2;
  Node M1_12_imag = lambda_imag / one_p_epsilon2;

  Node M1_21_real = (kappa - epsilon * (N2_a + 1) + N1_a - lambda_real) / one_p_epsilon2;
  Node M1_21_imag = -M1_12_imag;

  Node M1_22_real = (-epsilon * (kappa + N2_b - lambda_real) + N1_b - 1) / one_p_epsilon2;
  Node M1_22_imag = M1_11_imag;

  // M2
  Node M2_11 = -kappa / one_p_epsilon2 * epsilon * xi - 2 / xi;
  Node M2_12 = -kappa / one_p_epsilon2 * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[4] = Z_r_1;
  out[5] = Z_r_2;
  out[6] = Z_i_1;
  out[7] = Z_i_2;

  out[8] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[9] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[10] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[11] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[12] = 0 * lambda_real;
  out[13] = 0 * lambda_imag;
  out[14] = 0 * kappa;
  out[15] = 0 * epsilon;
}

// Vector field specialized for d = 3, sigma = 1, delta = 0, omega = 1
void vectorField_optimized_d3_fixed_epsilon(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node epsilon = params[0];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];

  Node Y_r_1 = in[4];
  Node Y_r_2 = in[5];
  Node Y_i_1 = in[6];
  Node Y_i_2 = in[7];
  Node Z_r_1 = in[8];
  Node Z_r_2 = in[9];
  Node Z_i_1 = in[10];
  Node Z_i_2 = in[11];

  Node lambda_real = in[12];
  Node lambda_imag = in[13];
  Node kappa = in[14];

  Node a2 = (a^2);
  Node b2 = (b^2);
  Node a2b2 = a2 + b2;

  // Compute forward ODE
  Node F1 = -2 / xi * (alpha + epsilon * beta) -
      kappa * xi * beta -
      kappa * b -
      a -
      a2b2 * a;

  Node F2 = -2 / xi * (beta - epsilon * alpha) +
      kappa * xi * alpha +
      kappa * a -
      b -
      a2b2 * b;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;

  // Compute linearized ODE

  // Compute J_N
  Node N1_a = -2 * a * b;
  Node N1_b = -(a2 + 3 * b2);
  Node N2_a = 3 * a2 + b2;
  Node N2_b = -N1_a;

  // The indexing notation _ij corresponds to the Julia indexing convention.

  // M1
  Node M1_11_real = (-epsilon * (kappa + N1_a - lambda_real) - N2_a - 1) / one_p_epsilon2;
  Node M1_11_imag = (epsilon * lambda_imag) / one_p_epsilon2;

  Node M1_12_real = (-kappa - epsilon * (N1_b - 1) - N2_b + lambda_real) / one_p_epsilon2;
  Node M1_12_imag = lambda_imag / one_p_epsilon2;

  Node M1_21_real = (kappa - epsilon * (N2_a + 1) + N1_a - lambda_real) / one_p_epsilon2;
  Node M1_21_imag = -M1_12_imag;

  Node M1_22_real = (-epsilon * (kappa + N2_b - lambda_real) + N1_b - 1) / one_p_epsilon2;
  Node M1_22_imag = M1_11_imag;

  // M2
  Node M2_11 = -kappa / one_p_epsilon2 * epsilon * xi - 2 / xi;
  Node M2_12 = -kappa / one_p_epsilon2 * xi;
  Node M2_21 = -M2_12;
  Node M2_22 = M2_11;

  out[4] = Z_r_1;
  out[5] = Z_r_2;
  out[6] = Z_i_1;
  out[7] = Z_i_2;

  out[8] = M1_11_real * Y_r_1 + M1_12_real * Y_r_2 - (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2) + M2_11 * Z_r_1 + M2_12 * Z_r_2;
  out[9] = M1_21_real * Y_r_1 + M1_22_real * Y_r_2 - (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2) + M2_21 * Z_r_1 + M2_22 * Z_r_2;
  out[10] = M1_11_imag * Y_r_1 + M1_12_imag * Y_r_2 + M1_11_real * Y_i_1 + M1_12_real * Y_i_2 + M2_11 * Z_i_1 + M2_12 * Z_i_2;
  out[11] = M1_21_imag * Y_r_1 + M1_22_imag * Y_r_2 + M1_21_real * Y_i_1 + M1_22_real * Y_i_2 + M2_21 * Z_i_1 + M2_22 * Z_i_2;

  out[12] = 0 * lambda_real;
  out[13] = 0 * lambda_imag;
  out[14] = 0 * kappa;
}

void Y_optimized_d3_fixed_epsilon(
    IVector Q_hat_0,
    IVector Y_0,
    interval lambda_real,
    interval lambda_imag,
    interval kappa,
    interval epsilon,
    interval T0,
    interval T1,
    double tol
) {
    IVector u0(15);

    u0[0] = Q_hat_0[0];
    u0[1] = Q_hat_0[1];
    u0[2] = Q_hat_0[2];
    u0[3] = Q_hat_0[3];

    u0[4] = Y_0[0];
    u0[5] = Y_0[1];
    u0[6] = Y_0[2];
    u0[7] = Y_0[3];
    u0[8] = Y_0[4];
    u0[9] = Y_0[5];
    u0[10] = Y_0[6];
    u0[11] = Y_0[7];

    u0[12] = lambda_real;
    u0[13] = lambda_imag;
    u0[14] = kappa;

    int dim = 15;
    IMap vf = IMap(vectorField_optimized_d3_fixed_epsilon, dim, dim, 1);
    vf.setParameter(0, epsilon);

    // Create the solver and the time map
    IOdeSolver solver(vf, 20);

    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    ITimeMap timeMap(solver);

  try {
      // Define a doubleton representation of the initial value
      C0HORect2Set s(u0, T0);

      // Solve the system
      IVector result = timeMap(T1, s);

      for (int i = 4; i < 12; i++)
          cout << result[i] << endl;
  } catch(exception& e) {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
}

int main()
{
  // Enough to exactly recover Float64 values
  cout.precision(17);
  cerr.precision(17);

  // Read initial data
  IVector Q_hat_0(4), Y_0(8);
  // Read Q_hat initial data
  cin >> Q_hat_0[0] >> Q_hat_0[1] >> Q_hat_0[2] >> Q_hat_0[3];
  // Read Y initial data
  cin >> Y_0[0] >> Y_0[1] >> Y_0[2] >> Y_0[3] >> Y_0[4] >> Y_0[5] >> Y_0[6] >> Y_0[7];

  // Read parameter values
  int d;
  interval lambda_real, lambda_imag, kappa, epsilon, omega, sigma, delta;
  cin >> d;
  cin >> lambda_real;
  cin >> lambda_imag;
  cin >> kappa;
  cin >> epsilon;
  cin >> omega;
  cin >> sigma;
  cin >> delta;

  // Read time span
  interval T0, T1;
  cin >> T0 >> T1;

  // Read tolerance to use
  double tol;
  cin >> tol;

  if (d == 3 && sigma == 1 && delta == 0 && omega == 1) {
      Y_optimized_d3_fixed_epsilon(Q_hat_0, Y_0, lambda_real, lambda_imag, kappa, epsilon, T0, T1, tol);
  } else {
      // TODO: Do we need to implement any other versions?
      cout << "Exception: Got unsupported parameters" << endl;
  }
}
