#pragma once
#include "../../BoundaryConditions.hpp"
#include "./ExplicitBase.hpp"
#include "./FiniteDifferenceHeatSolver.hpp"
#include <cmath>
#include <exception>
#include <libField/Field.hpp>

namespace BC = HeatSolvers::BoundaryConditions;
namespace FDHS = HeatSolvers::_2D::Cylindrical;

template <class REAL>
class NonUniformExplicit : public ExplicitBase<NonUniformExplicit<REAL>, REAL> {
public:
  // -------------------- PUBLIC CONSTRUCTORS --------------------
  NonUniformExplicit(size_t _zN, size_t _rN) {
    this->rN = _rN;
    this->zN = _zN;
    this->T.reset(_zN, _rN);
    this->A.reset(this->T.getCoordinateSystemPtr());
    this->VHC.reset(this->T.getCoordinateSystemPtr());
    this->k.reset(this->T.getCoordinateSystemPtr());
    this->T_prime.reset(this->T.getCoordinateSystemPtr());
  }
  NonUniformExplicit() {
    // Default constructor to be safe
  }
  // ---------------------- PUBLIC METHODS -----------------------
  int get_rdims() { return this->rN; }

  int get_zdims() { return this->zN; }

  /**
   * get the spacing at (i, j) ito r with either forward, backward, or central
   * differences
   */
  REAL get_dr(int i, int j, int o = -1) {
    REAL dr = 0.0;
    if (j == 0) {
      o = 1;
    } else if (j == this->rN - 1) {
      o = -1;
    }
    if (o == -1) {
      dr = this->k.getCoord(i, j)[1] - this->k.getCoord(i, j - 1)[1];
    } else if (o == 0) {
      dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j - 1)[1];
    } else if (o == 1) {
      dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
    }
    return dr;
  }

  /**
   * get the spacing at (i, j) ito z with either forward, backward, or central
   * differences
   */
  REAL get_dz(int i, int j, int o = -1) {
    REAL dz = 0.0;
    if (i == 0) {
      o = 1;
    } else if (i == this->zN - 1) {
      o = -1;
    }
    if (o == -1) {
      dz = this->k.getCoord(i, j)[0] - this->k.getCoord(i - 1, j)[0];
    } else if (o == 0) {
      dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i - 1, j)[0];
    } else if (o == 1) {
      dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
    }
    return dz;
  }

  /**
   * forward finite difference derivative of f wrt r
   */
  REAL forward_finite_r(Field<REAL, 2> &f, int i, int j) {
    return f[i][j + 1] - f[i][j] / get_dr(i, j, 1);
  }

  /**
   * backward finite difference derivative of f wrt r
   */
  REAL backward_finite_r(Field<REAL, 2> &f, int i, int j) {
    return f[i][j] - f[i][j - 1] / get_dr(i, j, -1);
  }

  /**
   * central finite difference derivative of f wrt r
   */
  REAL central_finite_r(Field<REAL, 2> &f, int i, int j) {
    REAL C1, C2, C3;
    C1 = this->get_dr(i, j, -1) /
         (this->get_dr(i, j, 0) * this->get_dr(i, j, 1));
    C2 = (this->get_dr(i, j, 1) - this->get_dr(i, j, -1)) /
         (this->get_dr(i, j, -1) * this->get_dr(i, j, 1));
    C3 = -1 * this->get_dr(i, j, 1) /
         (this->get_dr(i, j, 0) * this->get_dr(i, j, -1));
    return C1 * f[i][j + 1] + C2 * f[i][j] + C3 * f[i][j - 1];
  }

  /**
   * forward finite difference derivative of f wrt z
   */
  REAL forward_finite_z(Field<REAL, 2> &f, int i, int j) {
    return f[i + 1][j] - f[i][j] / get_dz(i, j, 1);
  }

  /**
   * backward finite difference derivative of f wrt z
   */
  REAL backward_finite_z(Field<REAL, 2> &f, int i, int j) {
    return f[i][j] - f[i - 1][j] / get_dz(i, j, -1);
  }

  /**
   * central finite difference derivative of f wrt z
   */
  REAL central_finite_z(Field<REAL, 2> &f, int i, int j) {
    REAL C1, C2, C3;
    C1 = this->get_dz(i, j, -1) /
         (this->get_dz(i, j, 0) * this->get_dz(i, j, 1));
    C2 = (this->get_dz(i, j, 1) - this->get_dz(i, j, -1)) /
         (this->get_dz(i, j, -1) * this->get_dz(i, j, 1));
    C3 = -1 * this->get_dz(i, j, 1) /
         (this->get_dz(i, j, 0) * this->get_dz(i, j, -1));
    return C1 * f[i + 1][j] + C2 * f[i][j] + C3 * f[i - 1][j];
  }

  /**
   * finite difference derivative of kappa wrt r
   */
  REAL dk_dr(int i, int j) {
    if (j == 0) {
      return forward_finite_r(this->k, i, j);
    } else if (j == this->rN - 1) {
      return backward_finite_r(this->k, i, j);
    } else {
      return central_finite_r(this->k, i, j);
    }
  }

  /**
   * finite difference derivative of kappa wrt z
   */
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
    REAL dr_l, dr_c, dr_r;
    REAL dz_l, dz_c, dz_r;
    REAL dk_z, dk_r;
    REAL r;
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dr_l = get_dr(i, j, -1);
    dr_c = get_dr(i, j, 0);
    dr_r = get_dr(i, j, 1);
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dz_l = get_dz(i, j, -1);
    dz_c = get_dz(i, j, 0);
    dz_r = get_dz(i, j, 1);
    // kappa derivatives
    dk_z = this->dk_dz(i, j);
    dk_r = this->dk_dr(i, j);

    r = this->T.getCoord(i, j)[1];

    REAL T1 = (1 / r) * this->k[i][j] * ((dr_r - dr_l) / (dr_r * dr_l));
    REAL T2 = -2 * this->k[i][j] * (1 / (dr_r * dr_l));
    REAL T3 = dk_r * ((dr_r - dr_l) / (dr_r * dr_l));
    REAL T4 = -2 * this->k[i][j] * (1 / (dz_r * dz_l));
    REAL T5 = dk_z * ((dz_r - dz_l) / (dz_r * dz_l));
    return T1 + T2 + T3 + T4 + T5;
  }

  // Calculate Coefficent for T^n_(z, r-1)
  REAL B_n(int i, int j) {
    REAL dr_l, dr_c, dr_r;
    REAL dk_r;
    REAL r;
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dr_l = get_dr(i, j, -1);
    dr_c = get_dr(i, j, 0);
    dr_r = get_dr(i, j, 1);
    // kappa derivatives
    dk_r = this->dk_dr(i, j);

    r = this->T.getCoord(i, j)[1];

    REAL T1 = (-1.0 / r) * this->k[i][j] * dr_r / (dr_c * dr_l);
    REAL T2 = (2.0 / (dr_c * dr_l)) * this->k[i][j];
    REAL T3 = -1.0 * dk_r * dr_r / (dr_c * dr_l);
    return T1 + T2 + T3;
  }

  // Calculate Coefficent for T^n_(z, r+1)
  REAL C_n(int i, int j) {
    REAL dr_l, dr_c, dr_r;
    REAL dk_r;
    REAL r;
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dr_l = get_dr(i, j, -1);
    dr_c = get_dr(i, j, 0);
    dr_r = get_dr(i, j, 1);
    // kappa derivatives
    dk_r = this->dk_dr(i, j);

    r = this->T.getCoord(i, j)[1];

    REAL T1 = (1.0 / r) * this->k[i][j] * dr_l / (dr_c * dr_r);
    REAL T2 = (2.0 / (dr_c * dr_r)) * this->k[i][j];
    REAL T3 = dk_r * dr_l / (dr_c * dr_r);
    return T1 + T2 + T3;
  }

  // Calculate Coefficent for T^n_(z-1, r)
  REAL D_n(int i, int j) {
    REAL dz_l, dz_c, dz_r;
    REAL dk_z;
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dz_l = get_dz(i, j, -1);
    dz_c = get_dz(i, j, 0);
    dz_r = get_dz(i, j, 1);
    // kappa derivatives
    dk_z = this->dk_dz(i, j);

    REAL T1 = 2 * this->k[i][j] * (1 / (dz_c * dz_l));
    REAL T2 = -1 * dk_z * (dz_r / (dz_c * dz_l));
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, r)
  REAL E_n(int i, int j) {
    REAL dz_l, dz_c, dz_r;
    REAL dk_z;
    // the subscripts 'l', 'c', 'r' are left center in place of '-' ' ' '+'
    dz_l = get_dz(i, j, -1);
    dz_c = get_dz(i, j, 0);
    dz_r = get_dz(i, j, 1);
    // kappa derivatives
    dk_z = this->dk_dz(i, j);

    REAL T1 = 2 * this->k[i][j] * (1 / (dz_c * dz_r));
    REAL T2 = dk_z * (dz_l / (dz_c * dz_r));
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL A_nR0(int i, int j) {
    REAL dr = get_dr(i, j, -1);
    REAL dz = get_dz(i, j);
    REAL T1 = (-4 * this->k[i][j]) / (dr * dr);
    REAL T2 = (-2 * this->k[i][j]) / (dz * dz);
    return T1 + T2;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL C_nR0(int i, int j) {
    REAL dr = get_dr(i, j, -1);
    // 2 * kappa / dr**2
    REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
    // 1 / 2 * dr
    REAL T2 = 2 * dr;
    // dk / dr
    REAL T3 = dk_dr(i, j) / dr;
    return T1 + T2 * T3;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL B_nR0(int i, int j) {
    REAL dr = get_dr(i, j, -1);
    // 2 * kappa / dr**2
    REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
    // 1 / 2 * dr
    REAL T2 = 2 * dr;
    // dk / dr
    REAL T3 = dk_dr(i, j) / dr;
    return T1 - T2 * T3;
  }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL D_nR0(int i, int j) { return D_n(i, j); }

  // Calculate Coefficent for T^n_(z+1, 0)
  REAL E_nR0(int i, int j) { return E_n(i, j); }

protected:
private:
};
