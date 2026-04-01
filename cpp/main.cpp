#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Problem {
  int n{};
  int max_iter{};
  double tol{};
  double omega{};
  std::vector<std::vector<double>> A;
  std::vector<double> b;
  std::vector<double> x0;
};

std::vector<std::string> read_nonempty_lines(const std::string &path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open input file: " + path);
  }
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(file, line)) {
    auto trimmed = line;
    trimmed.erase(trimmed.begin(),
                  std::find_if(trimmed.begin(), trimmed.end(),
                               [](unsigned char ch) { return !std::isspace(ch); }));
    trimmed.erase(
        std::find_if(trimmed.rbegin(), trimmed.rend(),
                     [](unsigned char ch) { return !std::isspace(ch); })
            .base(),
        trimmed.end());
    if (!trimmed.empty()) {
      lines.push_back(trimmed);
    }
  }
  return lines;
}

std::vector<double> parse_numbers(const std::string &line) {
  std::istringstream iss(line);
  std::vector<double> values;
  double value;
  while (iss >> value) {
    values.push_back(value);
  }
  return values;
}

Problem read_problem(const std::string &path) {
  auto lines = read_nonempty_lines(path);
  Problem problem;
  std::size_t cursor = 0;
  try {
    problem.n = std::stoi(lines.at(cursor++));
    problem.max_iter = std::stoi(lines.at(cursor++));
    problem.tol = std::stod(lines.at(cursor++));
    problem.A.reserve(problem.n);
    for (int i = 0; i < problem.n; ++i) {
      problem.A.push_back(parse_numbers(lines.at(cursor++)));
    }
    problem.b = parse_numbers(lines.at(cursor++));
    problem.x0 = parse_numbers(lines.at(cursor++));
    problem.omega = std::stod(lines.at(cursor++));
  } catch (const std::exception &e) {
    throw std::runtime_error("Invalid or incomplete input file: " + std::string(e.what()));
  }

  if (static_cast<int>(problem.b.size()) != problem.n ||
      static_cast<int>(problem.x0.size()) != problem.n) {
    throw std::runtime_error("Vectors b and x0 must match matrix size.");
  }
  for (const auto &row : problem.A) {
    if (static_cast<int>(row.size()) != problem.n) {
      throw std::runtime_error("Matrix A must be square.");
    }
  }
  if (problem.max_iter <= 0) {
    throw std::runtime_error("max_iter must be positive.");
  }
  if (problem.tol <= 0.0) {
    throw std::runtime_error("Tolerance must be positive.");
  }
  if (!(0.0 < problem.omega && problem.omega < 2.0)) {
    throw std::runtime_error("Relaxation factor w must be in the interval (0, 2).");
  }

  return problem;
}

double relative_error(const std::vector<double> &current,
                      const std::vector<double> &previous) {
  double err = 0.0;
  for (std::size_t i = 0; i < current.size(); ++i) {
    double denom = std::max(std::abs(current[i]), 1e-15);
    err = std::max(err, std::abs(current[i] - previous[i]) / denom);
  }
  return err;
}

std::vector<std::vector<double>> jacobi(const Problem &p) {
  std::vector<double> x_old = p.x0;
  std::vector<std::vector<double>> history;

  for (int iter = 0; iter < p.max_iter; ++iter) {
    std::vector<double> x(p.n, 0.0);
    for (int i = 0; i < p.n; ++i) {
      double diag = p.A[i][i];
      if (std::abs(diag) < 1e-15) {
        throw std::runtime_error("Zero diagonal encountered in Jacobi method.");
      }
      double sum = 0.0;
      for (int j = 0; j < p.n; ++j) {
        if (j == i) continue;
        sum += p.A[i][j] * x_old[j];
      }
      x[i] = (p.b[i] - sum) / diag;
    }
    history.push_back(x);
    if (relative_error(x, x_old) < p.tol) {
      break;
    }
    x_old = x;
  }
  return history;
}

std::vector<std::vector<double>> gauss_seidel(const Problem &p) {
  std::vector<double> x = p.x0;
  std::vector<std::vector<double>> history;

  for (int iter = 0; iter < p.max_iter; ++iter) {
    std::vector<double> previous = x;
    for (int i = 0; i < p.n; ++i) {
      double diag = p.A[i][i];
      if (std::abs(diag) < 1e-15) {
        throw std::runtime_error("Zero diagonal encountered in Gauss-Seidel method.");
      }
      double sum = 0.0;
      for (int j = 0; j < p.n; ++j) {
        if (j < i) {
          sum += p.A[i][j] * x[j];
        } else if (j > i) {
          sum += p.A[i][j] * previous[j];
        }
      }
      x[i] = (p.b[i] - sum) / diag;
    }
    history.push_back(x);
    if (relative_error(x, previous) < p.tol) {
      break;
    }
  }
  return history;
}

std::vector<std::vector<double>> sor(const Problem &p) {
  std::vector<double> x = p.x0;
  std::vector<std::vector<double>> history;

  for (int iter = 0; iter < p.max_iter; ++iter) {
    std::vector<double> previous = x;
    for (int i = 0; i < p.n; ++i) {
      double diag = p.A[i][i];
      if (std::abs(diag) < 1e-15) {
        throw std::runtime_error("Zero diagonal encountered in SOR method.");
      }
      double sum = 0.0;
      for (int j = 0; j < p.n; ++j) {
        if (j < i) {
          sum += p.A[i][j] * x[j];
        } else if (j > i) {
          sum += p.A[i][j] * previous[j];
        }
      }
      x[i] = (1.0 - p.omega) * previous[i] + p.omega * (p.b[i] - sum) / diag;
    }
    history.push_back(x);
    if (relative_error(x, previous) < p.tol) {
      break;
    }
  }
  return history;
}

void write_results(const std::string &path, const std::vector<std::vector<double>> &jacobi_hist,
                   const std::vector<std::vector<double>> &gs_hist,
                   const std::vector<std::vector<double>> &sor_hist) {
  std::ofstream out(path);
  if (!out.is_open()) {
    throw std::runtime_error("Unable to open output file: " + path);
  }
  out << std::fixed << std::setprecision(6);

  out << "Jacobi Method Result\n";
  for (std::size_t i = 0; i < jacobi_hist.size(); ++i) {
    out << std::setw(4) << (i + 1) << " ";
    for (double v : jacobi_hist[i]) {
      out << std::setw(12) << v << " ";
    }
    out << "\n";
  }

  out << "\nGauss-Seidel Method Result\n";
  for (std::size_t i = 0; i < gs_hist.size(); ++i) {
    out << std::setw(4) << (i + 1) << " ";
    for (double v : gs_hist[i]) {
      out << std::setw(12) << v << " ";
    }
    out << "\n";
  }

  out << "\nSuccessive Over-Relaxation (SOR) Result\n";
  for (std::size_t i = 0; i < sor_hist.size(); ++i) {
    out << std::setw(4) << (i + 1) << " ";
    for (double v : sor_hist[i]) {
      out << std::setw(12) << v << " ";
    }
    out << "\n";
  }
}
}  // namespace

int main(int argc, char *argv[]) {
  try {
    std::string input_path = argc > 1 ? argv[1] : "data/input.txt";
    std::string output_path = argc > 2 ? argv[2] : "data/output_cpp.txt";

    Problem problem = read_problem(input_path);
    auto jacobi_hist = jacobi(problem);
    auto gs_hist = gauss_seidel(problem);
    auto sor_hist = sor(problem);
    write_results(output_path, jacobi_hist, gs_hist, sor_hist);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
