# Jacobi, Gauss-Seidel, and SOR Methods

Modern, editor-agnostic implementations of the Jacobi, Gauss-Seidel, and Successive Over-Relaxation (SOR) iterative methods.

## Repository layout

```
data/
  input.txt            # sample input
  sample_output.txt    # example output
fortran/
  iterative_methods.f90
  main.f90
python/
  solve.py
cpp/
  main.cpp             # optional C++17 implementation
```

All Visual Studio specific files and build outputs have been removed. The code is designed to work in VS Code with common compilers/interpreters.

## Input format

The solver reads from a text file (`data/input.txt` by default) with the following layout:

```
n
max_iterations
tolerance
A matrix (n lines, n values each)
b vector (1 line, n values)
initial guess x0 (1 line, n values)
relaxation factor w (for SOR, 0 < w < 2)
```

`data/input.txt` contains the example from the original project.

## Fortran (Fortran 2023)

The sources use the `.f90` free-form extension, which remains the standard for
Fortran 2003/2008/2018/2023 code. Build with `gfortran` in Fortran 2023 mode
(GCC 14+ adds `-std=f2023`):

```bash
mkdir -p build
gfortran -std=f2023 -Wall -Wextra -O2 fortran/iterative_methods.f90 fortran/main.f90 -o build/iterative_solvers
./build/iterative_solvers data/input.txt data/output_fortran.txt
```

If your compiler does not yet recognize `-std=f2023` (e.g., GCC 13), use
`-std=f2018` as a drop-in fallback; the source is written to be Fortran 2023
compliant.

## Python

Run directly with Python 3 (no external dependencies):

```bash
python python/solve.py data/input.txt data/output_python.txt
```

## C++17 (optional)

```bash
mkdir -p build
g++ -std=c++17 -O2 -Wall -Wextra cpp/main.cpp -o build/iterative_solvers_cpp
./build/iterative_solvers_cpp data/input.txt data/output_cpp.txt
```

## Sample output

`data/sample_output.txt` shows the expected shape of the output produced by each implementation.
