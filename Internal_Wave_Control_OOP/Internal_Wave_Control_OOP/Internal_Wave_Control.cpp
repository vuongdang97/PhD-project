#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <matplotlibcpp.h> // Optional for visualization, requires matplotlib-cpp

namespace plt = matplotlibcpp;
using namespace Eigen;

// Class for solving the coupled wave problem
class CoupledWaveFVM {
public:
    CoupledWaveFVM(double a, double b, double T0, double T, int n, int m, double alpha, double gamma, double mu,
        const std::vector<double>& bl, const std::vector<double>& br);
    std::pair<MatrixXd, MatrixXd> solve(const std::vector<double>& y0, const std::vector<double>& y0t,
        const VectorXd& z, const MatrixXd& yhat, const VectorXd& lam_T);

private:
    double a, b, T0, T, alpha, gamma, mu, dx, dt;
    int n, m;
    std::vector<double> x, t, bl, br;
};

CoupledWaveFVM::CoupledWaveFVM(double a, double b, double T0, double T, int n, int m, double alpha, double gamma, double mu,
    const std::vector<double>& bl, const std::vector<double>& br)
    : a(a), b(b), T0(T0), T(T), n(n), m(m), alpha(alpha), gamma(gamma), mu(mu), bl(bl), br(br) {
    dx = (b - a) / (n - 1);
    dt = (T - T0) / (m - 1);

    // Generate x and t grids
    x.resize(n);
    t.resize(m);
    for (int i = 0; i < n; ++i) x[i] = a + i * dx;
    for (int i = 0; i < m; ++i) t[i] = T0 + i * dt;
}

std::pair<MatrixXd, MatrixXd> CoupledWaveFVM::solve(const std::vector<double>& y0, const std::vector<double>& y0t,
    const VectorXd& z, const MatrixXd& yhat, const VectorXd& lam_T) {
    double s = std::pow(dt / dx, 2);
    double s1 = dt / dx;

    SparseMatrix<double> A(n - 2, n - 2);
    SparseMatrix<double> C(n - 2, n - 2);
    VectorXd main_diag = VectorXd::Constant(n - 2, 2 * (1 - s));
    VectorXd off_diag = VectorXd::Constant(n - 3, s);
    for (int i = 0; i < n - 3; ++i) {
        A.insert(i, i) = main_diag[i];
        A.insert(i, i + 1) = off_diag[i];
        A.insert(i + 1, i) = off_diag[i];
    }

    // Initialize Yhat
    MatrixXd Yhat = MatrixXd::Zero(2 * (n - 2) * m, 1);
    for (int i = 0; i < n - 2; ++i) {
        Yhat(i, 0) = y0[i + 1];
        Yhat(i + n - 2, 0) = dt * y0t[i + 1];
    }

    // Initialize Z
    VectorXd Z = VectorXd::Zero(2 * (n - 2) * m);
    Z.segment((2 * m - 1) * (n - 2), n - 2) = lam_T.segment(1, n - 2);

    // Solve the sparse linear system
    SparseMatrix<double> M(2 * (n - 2) * m, 2 * (n - 2) * m);
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(M);
    VectorXd GrandVector = solver.solve(Yhat + Z);

    // Extract Y and Lambda
    MatrixXd Y = MatrixXd::Zero(n, m);
    MatrixXd Lambda = MatrixXd::Zero(n, m);
    for (int i = 0; i < m; ++i) {
        Y.block(1, i, n - 2, 1) = GrandVector.segment(i * (n - 2), n - 2);
        Lambda.block(1, i, n - 2, 1) = GrandVector.segment((i + m) * (n - 2), n - 2);
    }

    return { Y, Lambda };
}

// Class for error analysis
class ErrorAnalysis {
public:
    static std::tuple<double, double, double> compute_error(const MatrixXd& Y, const MatrixXd& y_exact,
        const MatrixXd& Lambda, const MatrixXd& L_exact,
        double dt, double dx);
};

std::tuple<double, double, double> ErrorAnalysis::compute_error(const MatrixXd& Y, const MatrixXd& y_exact,
    const MatrixXd& Lambda, const MatrixXd& L_exact,
    double dt, double dx) {
    double error_y = std::sqrt(dt * dx) * (Y - y_exact).norm();
    double error_max = (Y - y_exact).cwiseAbs().maxCoeff();
    double error_control = std::sqrt(dt * dx) * (Lambda - L_exact).norm();
    return { error_y, error_max, error_control };
}

// Main function
int main() {
    double a = 0, b = 1, T0 = 0;
    std::vector<int> N = { 51, 101, 201, 401 };
    std::vector<int> M = { 51, 101, 201, 401 };
    double alpha = 1, gamma = 0, mu = 1;

    std::vector<double> error, error_max, error_control;

    for (size_t l = 0; l < N.size(); ++l) {
        int n = N[l];
        int m = M[l];
        double dx = (b - a) / (n - 1);
        double dt = 1.0 / (m - 1);
        std::vector<double> bl(m, 0), br(m, 0);
        std::vector<double> x(n), t(m);
        for (int i = 0; i < n; ++i) x[i] = a + i * dx;
        for (int i = 0; i < m; ++i) t[i] = T0 + i * dt;

        // Initial conditions
        std::vector<double> y0(n), y0t(n);
        for (int i = 0; i < n; ++i) {
            y0[i] = -std::sin(M_PI * x[i]);
            y0t[i] = std::sin(M_PI * x[i]);
        }

        // Placeholder arrays for yhat and z
        MatrixXd yhat = MatrixXd::Zero(n, m);
        VectorXd z = VectorXd::Zero(n);
        VectorXd lam_T = VectorXd::Zero(n);

        CoupledWaveFVM solver(a, b, T0, 1.0, n, m, alpha, gamma, mu, bl, br);
        auto [Y, Lambda] = solver.solve(y0, y0t, z, yhat, lam_T);

        // Error analysis
        MatrixXd y_exact = MatrixXd::Zero(n, m);
        MatrixXd L_exact = MatrixXd::Zero(n, m);
        auto [err, err_max_val, err_control_val] =
            ErrorAnalysis::compute_error(Y, y_exact, Lambda, L_exact, dt, dx);

        error.push_back(err);
        error_max.push_back(err_max_val);
        error_control.push_back(err_control_val);
    }

    return 0;
}
