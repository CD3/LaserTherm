#pragma once
#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"
#include "LaserTherm/HeatSolvers/_2D/Cylindrical/ExplicitBase.hpp"
#include "LaserTherm/HeatSolvers/_2D/Cylindrical/FiniteDifferenceHeatSolver.hpp"
#include <cmath>
#include <exception>
#include <libField/Field.hpp>

namespace BC = HeatSolvers::BoundaryConditions;
namespace FDHS = HeatSolvers::_2D::Cylindrical;

template <class REAL>
class UniformExplicit : public ExplicitBase<UniformExplicit<REAL>, REAL> {
public:
  // -------------------- PUBLIC CONSTRUCTORS --------------------
  UniformExplicit(size_t _zN, size_t _rN) {
    this->rN = _rN;
    this->zN = _zN;
    this->T.reset(_zN, _rN);
    this->A.reset(this->T.getCoordinateSystemPtr());
    this->VHC.reset(this->T.getCoordinateSystemPtr());
    this->k.reset(this->T.getCoordinateSystemPtr());
    this->T_prime.reset(this->T.getCoordinateSystemPtr());
  }
  UniformExplicit() {
    // Default constructor to be safe
  }

  // get forward r spacing (backwards at rmax)
  REAL get_dr(int i, int j) {
    REAL dr;
    if (j == this->rN - 1) {
      dr = this->k.getCoord(i, j)[1] - this->k.getCoord(i, j - 1)[1];
    } else {
      dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
    }
    return dr;
  }

  // get forward z spacing (backwards at zmax)
  REAL get_dz(int i, int j) {
    REAL dz;
    if (i == this->zN - 1) {
      dz = this->k.getCoord(i, j)[0] - this->k.getCoord(i - 1, j)[0];
    } else {
      dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
    }
    return dz;
  }

  REAL forward_finite_r(Field<REAL, 2> &f, int i, int j) {
    return f[i][j + 1] - f[i][j] / this->get_dr(i, j);
  }
  REAL backward_finite_r(Field<REAL, 2> &f, int i, int j) {
    return f[i][j] - f[i][j - 1] / this->get_dr(i, j);
  }
  REAL central_finite_r(Field<REAL, 2> &f, int i, int j) {
    return f[i][j + 1] - f[i][j - 1] / (2 * this->get_dr(i, j));
  }
  REAL forward_finite_z(Field<REAL, 2> &f, int i, int j) {
    return f[i + 1][j] - f[i][j] / this->get_dz(i, j);
  }
  REAL backward_finite_z(Field<REAL, 2> &f, int i, int j) {
    return f[i][j] - f[i - 1][j] / this->get_dz(i, j);
  }
  REAL central_finite_z(Field<REAL, 2> &f, int i, int j) {
    return f[i + 1][j] - f[i - 1][j] / (2 * this->get_dz(i, j));
  }
  // Finite Difference deriv for kappa wrt r
  REAL dk_dr(int i, int j) {
    if (j == 0) {
      return forward_finite_r(this->k, i, j);
    } else if (j == this->rN - 1) {
      return backward_finite_r(this->k, i, j);
    } else {
      return central_finite_r(this->k, i, j);
    }
  }
  // Finite Difference deriv for kappa wrt z
  REAL dk_dz(int i, int j) {
    if (i == 0) {
      return forward_finite_z(this->k, i, j);
    } else if (i == this->zN - 1) {
      return backward_finite_z(this->k, i, j);
    } else {
      return central_finite_z(this->k, i, j);
    }
  }
  // Calculate Coefficent for T^n_(z, r)
  REAL A_n(int i, int j) {
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    // -2 * kappa / dr**2
    REAL T1 = -2 * this->k[i][j] / pow(dr, 2);
    // -2 * kappa / dz**2
    REAL T2 = -2 * this->k[i][j] / pow(dz, 2);
    return T1 + T2;
  }
  // Calculate Coefficent for T^n_(z, r+1)
  REAL C_n(int i, int j) {
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    // kappa / 2 r dr
    REAL T1 = this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
    // kappa / dr**2
    REAL T2 = this->k[i][j] / pow(dr, 2);
    // dk_dr / 2 dr
    REAL T3 = dk_dr(i, j) / (2 * dr);
    return T1 + T2 + T3;
  }

  // Calculate Coefficent for T^n_(z, r-1)
  REAL B_n(int i, int j) {
    // return (N) 0;
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    // - kappa / 2r * dr
    REAL T1 = -this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
    // kappa / dr**2
    REAL T2 = this->k[i][j] / pow(dr, 2);
    // dk_dr / 2 dr
    REAL T3 = -this->dk_dr(i, j) / (2 * dr);
    return T1 + T2 + T3;
  }

  // Calculate Coefficent for T^n_(z-1, r)
  REAL D_n(int i, int j) {
    // return (N) 0;
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    // kappa / dz**2
    REAL T1 = this->k[i][j] / pow(dz, 2);
    // - dz_dr / 2 dz
    REAL T2 = -dk_dz(i, j) / (2 * dz);
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, r)
  REAL E_n(int i, int j) {
    // return (N) 0;
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    // kappa / dz**2
    REAL T1 = this->k[i][j] / pow(dz, 2);
    // - dz_dr / 2 dz
    REAL T2 = dk_dz(i, j) / (2 * dz);
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL A_nR0(int i, int j) {
    REAL dr = this->get_dr(i, j);
    REAL dz = this->get_dz(i, j);
    REAL T1 = (-4 * this->k[i][j]) / (dr * dr);
    REAL T2 = (-2 * this->k[i][j]) / (dz * dz);
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL C_nR0(int i, int j) {
    REAL dr = this->get_dr(i, j);
    // 2 * kappa / dr**2
    REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
    // 1 / 2 * dr
    REAL T2 = 2 * dr;
    // dk / dr
    REAL T3 = this->dk_dr(i, j) / dr;
    return T1 + T2 * T3;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL B_nR0(int i, int j) {
    REAL dr = this->get_dr(i, j);
    // 2 * kappa / dr**2
    REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
    // 1 / 2 * dr
    REAL T2 = 2 * dr;
    // dk / dr
    REAL T3 = this->dk_dr(i, j) / dr;
    return T1 - T2 * T3;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL D_nR0(int i, int j) { return D_n(i, j); }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL E_nR0(int i, int j) { return E_n(i, j); }

protected:
private:
};
