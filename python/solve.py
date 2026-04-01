from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple


def _read_nonempty_lines(path: Path) -> List[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def read_problem(path: Path) -> Tuple[int, float, List[List[float]], List[float], List[float], float]:
    lines = _read_nonempty_lines(path)
    cursor = 0

    try:
        n = int(lines[cursor])
        cursor += 1
        max_iterations = int(lines[cursor])
        cursor += 1
        tol = float(lines[cursor])
        cursor += 1

        matrix: List[List[float]] = []
        for _ in range(n):
            matrix.append([float(val) for val in lines[cursor].split()])
            cursor += 1

        b = [float(val) for val in lines[cursor].split()]
        cursor += 1
        x0 = [float(val) for val in lines[cursor].split()]
        cursor += 1
        omega = float(lines[cursor])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Invalid or incomplete input file: {path}") from exc

    if any(len(row) != n for row in matrix):
        raise ValueError("Matrix A must be square and match the size of vector b.")
    if len(b) != n or len(x0) != n:
        raise ValueError("Vectors b and x0 must be the same length as matrix rows.")
    if max_iterations <= 0:
        raise ValueError("max_iterations must be positive.")
    if tol <= 0.0:
        raise ValueError("Tolerance must be positive.")
    if not (0.0 < omega < 2.0):
        raise ValueError("Relaxation factor w must be in the interval (0, 2).")

    return max_iterations, tol, matrix, b, x0, omega


def relative_error(new_vals: Sequence[float], old_vals: Sequence[float]) -> float:
    err = 0.0
    for new, old in zip(new_vals, old_vals):
        denom = max(abs(new), 1e-15)
        err = max(err, abs(new - old) / denom)
    return err


def jacobi(
    a: Sequence[Sequence[float]],
    b: Sequence[float],
    x0: Sequence[float],
    max_iterations: int,
    tol: float,
) -> List[List[float]]:
    x_old = list(x0)
    n = len(b)
    history: List[List[float]] = []

    for _ in range(max_iterations):
        x = []
        for i in range(n):
            diagonal = a[i][i]
            if abs(diagonal) < 1e-15:
                raise ZeroDivisionError(f"Zero diagonal encountered at row {i + 1} in Jacobi method.")
            summation = sum(a[i][j] * x_old[j] for j in range(n) if j != i)
            x.append((b[i] - summation) / diagonal)
        history.append(x)
        if relative_error(x, x_old) < tol:
            break
        x_old = x
    return history


def gauss_seidel(
    a: Sequence[Sequence[float]],
    b: Sequence[float],
    x0: Sequence[float],
    max_iterations: int,
    tol: float,
) -> List[List[float]]:
    x = list(x0)
    n = len(b)
    history: List[List[float]] = []

    for _ in range(max_iterations):
        x_old = x.copy()
        for i in range(n):
            diagonal = a[i][i]
            if abs(diagonal) < 1e-15:
                raise ZeroDivisionError(f"Zero diagonal encountered at row {i + 1} in Gauss-Seidel method.")
            summation = 0.0
            for j in range(n):
                if j < i:
                    summation += a[i][j] * x[j]
                elif j > i:
                    summation += a[i][j] * x_old[j]
            x[i] = (b[i] - summation) / diagonal
        history.append(x.copy())
        if relative_error(x, x_old) < tol:
            break
    return history


def sor(
    a: Sequence[Sequence[float]],
    b: Sequence[float],
    x0: Sequence[float],
    max_iterations: int,
    tol: float,
    omega: float,
) -> List[List[float]]:
    x = list(x0)
    n = len(b)
    history: List[List[float]] = []

    for _ in range(max_iterations):
        x_old = x.copy()
        for i in range(n):
            diagonal = a[i][i]
            if abs(diagonal) < 1e-15:
                raise ZeroDivisionError(f"Zero diagonal encountered at row {i + 1} in SOR method.")
            summation = 0.0
            for j in range(n):
                if j < i:
                    summation += a[i][j] * x[j]
                elif j > i:
                    summation += a[i][j] * x_old[j]
            x[i] = (1 - omega) * x_old[i] + omega * (b[i] - summation) / diagonal
        history.append(x.copy())
        if relative_error(x, x_old) < tol:
            break
    return history


def write_results(
    output_path: Path,
    jacobi_hist: Sequence[Sequence[float]],
    gs_hist: Sequence[Sequence[float]],
    sor_hist: Sequence[Sequence[float]],
) -> None:
    with output_path.open("w", encoding="utf-8") as f:
        f.write("Jacobi Method Result\n")
        for idx, values in enumerate(jacobi_hist, start=1):
            f.write(f"{idx:4d} " + " ".join(f"{val:12.6f}" for val in values) + "\n")

        f.write("\nGauss-Seidel Method Result\n")
        for idx, values in enumerate(gs_hist, start=1):
            f.write(f"{idx:4d} " + " ".join(f"{val:12.6f}" for val in values) + "\n")

        f.write("\nSuccessive Over-Relaxation (SOR) Result\n")
        for idx, values in enumerate(sor_hist, start=1):
            f.write(f"{idx:4d} " + " ".join(f"{val:12.6f}" for val in values) + "\n")


def solve_file(input_path: Path, output_path: Path) -> None:
    max_iterations, tol, a, b, x0, omega = read_problem(input_path)
    jacobi_hist = jacobi(a, b, x0, max_iterations, tol)
    gs_hist = gauss_seidel(a, b, x0, max_iterations, tol)
    sor_hist = sor(a, b, x0, max_iterations, tol, omega)
    write_results(output_path, jacobi_hist, gs_hist, sor_hist)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Solve linear systems with Jacobi, Gauss-Seidel, and SOR methods.")
    parser.add_argument(
        "input",
        nargs="?",
        default="data/input.txt",
        type=Path,
        help="Path to the input file (default: data/input.txt)",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="data/output_python.txt",
        type=Path,
        help="Path for the output file (default: data/output_python.txt)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    solve_file(args.input, args.output)


if __name__ == "__main__":
    main()
