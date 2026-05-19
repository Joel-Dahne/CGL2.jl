#include "capd/capdlib.h"
#include <vector>

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// ===========================================================================
// VECTOR FIELD DEFINITIONS
// ===========================================================================

// Generic version of the vector field as a first-order real system.
// TODO: More description?
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

// Specialization for d == 1: the -(d-1)/xi * (...) singular term vanishes.
// TODO: More description?
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

// specialization for d == 3, omega == 1, sigma == 1, delta == 0.
// TODO: More description?
void vectorField_d3_optimized(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node* /*params*/, int /*noParams*/)
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

// ===========================================================================
// SHARED LOGIC
// ===========================================================================

// Construct the most specialized IMap for the given parameters.
IMap build_vector_field(int d, interval omega, interval sigma, interval delta, interval epsilon) {
    int dim = 16;
    IMap vf;

    bool d3_opt = (d == 3) && (omega == 1) && (sigma == 1) && (delta == 0);

    if (d == 1) {
	// Generic version for d = 1
        vf = IMap(vectorField_d1, dim, dim, 3);
        vf.setParameter(0, omega);
	vf.setParameter(1, sigma);
	vf.setParameter(2, delta);
    } else if (d3_opt) {
	// Specialized for the case considered in the paper
        vf = IMap(vectorField_d3_optimized, dim, dim, 0);
    } else {
	// Generic version
        vf = IMap(vectorField, dim, dim, 4);
        vf.setParameter(0, omega);
	vf.setParameter(1, sigma);
	vf.setParameter(2, delta);
	vf.setParameter(3, interval(d));
    }

    return vf;
}

// Integrate the CGL ODE from xi_0 to xi_1 and print the result to stdout.
// Outputs [a, b, alpha, beta] at xi_1 as intervals, one per line.
// Returns 0 on success, 1 on solver error (outputs NaN intervals).
int Y_zero(
    /** Initial value Q_hat(xi_0) **/
    IVector Q_hat_xi_0,
    /** Initial value Y(xi_0) **/
    IVector Y_xi_0,
    /** Parameters **/
    interval lambda_real,
    interval lambda_imag,
    int d,
    interval kappa,
    interval epsilon,
    interval omega,
    interval sigma,
    interval delta,
    /** Integration interval **/
    interval xi_0,
    interval xi_1,
    /** Settings **/
    double tol
) {
    // To get better enclosures for wide kappa and epsilon we treat
    // them as variables in the ODE.
    IVector Q_hat_Y_xi_0_with_parameters(16);
    for (int i = 0; i < 4; i++)
	Q_hat_Y_xi_0_with_parameters[i] = Q_hat_xi_0[i];
    for (int i = 0; i < 8; i++)
	Q_hat_Y_xi_0_with_parameters[4 + i] = Y_xi_0[i];
    Q_hat_Y_xi_0_with_parameters[12] = lambda_real;
    Q_hat_Y_xi_0_with_parameters[13] = lambda_imag;
    Q_hat_Y_xi_0_with_parameters[14] = kappa;
    Q_hat_Y_xi_0_with_parameters[15] = epsilon;

    // Choose optimized vector field based on parameters
    IMap vf = build_vector_field(d, omega, sigma, delta, epsilon);

    // Create the solver and the time map
    IOdeSolver solver(vf, 20);

    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    ITimeMap timeMap(solver);

    try {
	// Initial value for solver
	C0HORect2Set s(Q_hat_Y_xi_0_with_parameters, xi_0);

	// Solve the system
	IVector result = timeMap(xi_1, s);

	for (int i = 4; i < 12; i++) {
	    cout << result[i] << endl;
	}

	return 0; // Success
    } catch(...) {
	for (int i = 4; i < 12; i++) {
	    cout << interval(NAN) << endl;
	}

	return 1; // Error
    }
}

int main()
{
  cout.precision(17); // Enough to exactly recover Float64 values
  cerr.precision(17); // Enough to exactly recover Float64 values

  // Read input from stdin

  // Read initial data
  IVector Q_hat_xi_0(4), Y_xi_0(8);
  // Read Q_hat initial data
  cin >> Q_hat_xi_0[0] >> Q_hat_xi_0[1] >> Q_hat_xi_0[2] >> Q_hat_xi_0[3];
  // Read Y initial data
  cin >> Y_xi_0[0] >> Y_xi_0[1] >> Y_xi_0[2] >> Y_xi_0[3] >> Y_xi_0[4] >> Y_xi_0[5] >> Y_xi_0[6] >> Y_xi_0[7];

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
  interval xi_0, xi_1;
  cin >> xi_0 >> xi_1;

  // Read tolerance to use
  double tol;
  cin >> tol;

  return Y_zero(Q_hat_xi_0, Y_xi_0, lambda_real, lambda_imag, d, kappa, epsilon, omega, sigma, delta, xi_0, xi_1, tol);
}
