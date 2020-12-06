#ifndef LaserTherm_Utils_TriDiagonalSolver_hpp
#define LaserTherm_Utils_TriDiagonalSolver_hpp

/** @file TriDiagonalSolver.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/14/19
 */

enum SolverAlgorithm { Thomas };

template<SolverAlgorithm algorithm>
struct TriDiagonalSolver {
};

/**
 * @brief Thomas Algorithm solver
 *
 * Solves the matrix equation
 *
 * \f$A x = b\f$
 *
 * when \f$A\f$ is tridiagonal using the Thomas algorithm.
 *
 * All three diagonals are assumed to be of length N. The first element of aSub
 * and the last element of aSup are not used.
 *
 *  @param	aSub   The subdiagonal (aSub(1) is the first sub diagonal
 * element)
 *  @param	aDiag  The diagonal.
 *  @param 	aSup   The superdiagonal (aSup(N-2) is the last sub diagonal
 * element)
 *  @param	_b      Know part of the system.
 *  @param	_x      Solution, this vector will be overwritten.
 *
 * Note: Inputs are destoyed
 */
template<>
struct TriDiagonalSolver<Thomas> {
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  static void Solve(T1& aSub, T2& aDiag, T3& aSup, T4& b, T5& x)
  {
    int N = x.size();

    // FORWARD SWEEP
    for (int i = 1; i < N;
         ++i)  // can't parallelize, each iteration depends on the last
    {
      aDiag(i) = aDiag(i) - aSub(i) * aSup(i - 1) / aDiag(i - 1);
      b(i)     = b(i) - aSub(i) * b(i - 1) / aDiag(i - 1);
    }

    // BACK SUBSTITUTE
    x(N - 1) = b(N - 1) / aDiag(N - 1);
    for (int i = (N - 2); i >= 0;
         --i)  // can't parallelize, each iteration depends on the last
    {
      x(i) = (b(i) - aSup(i) * x(i + 1)) / aDiag(i);
    }

    return;
  }
};

#endif  // include protector
