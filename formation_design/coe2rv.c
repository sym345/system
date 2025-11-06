/*
 * File: coe2rv.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

/* Include Files */
#include "coe2rv.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * COE2RV Converts classical orbital elements to position and velocity vectors
 *
 *  Inputs:
 *    k     - Gravitational parameter (km^3/s^2)
 *    p     - Semi-latus rectum (km)
 *    ecc   - Eccentricity
 *    inc   - Inclination (rad)
 *    raan  - Right ascension of ascending node (rad)
 *    argp  - Argument of perigee (rad)
 *    nu    - True anomaly (rad)
 *
 *  Outputs:
 *    r_ijk - Position vector in IJK frame (km)
 *    v_ijk - Velocity vector in IJK frame (km/s)
 *
 * Arguments    : double p
 *                double ecc
 *                double inc
 *                double raan
 *                double argp
 *                double nu
 *                double r_ijk[3]
 *                double v_ijk[3]
 * Return Type  : void
 */
void coe2rv(double p, double ecc, double inc, double raan, double argp,
            double nu, double r_ijk[3], double v_ijk[3])
{
  double Q_pqw_to_ijk[9];
  double c_R3_W_tmp[9];
  double d_R3_W_tmp[9];
  double c_a[3];
  double R1_i_tmp;
  double R3_W_tmp;
  double R3_w_tmp;
  double a;
  double b_R1_i_tmp;
  double b_R3_W_tmp;
  double b_R3_w_tmp;
  double b_a;
  double b_r_pqw_tmp;
  double d;
  double r_pqw_tmp;
  int i;
  int i1;
  int i2;
  /*  --- Step 1: Compute r and v in the PQW (perifocal) frame */
  r_pqw_tmp = cos(nu);
  b_r_pqw_tmp = sin(nu);
  a = p / (ecc * r_pqw_tmp + 1.0);
  b_a = sqrt(3.986004418E+14 / p);
  /*  --- Step 2: Build rotation matrix from PQW to IJK */
  /*  This is a 3-1-3 rotation: Rz(-raan) * Rx(-inc) * Rz(-argp) */
  R3_W_tmp = sin(-raan);
  b_R3_W_tmp = cos(-raan);
  R1_i_tmp = sin(-inc);
  b_R1_i_tmp = cos(-inc);
  R3_w_tmp = sin(-argp);
  b_R3_w_tmp = cos(-argp);
  c_R3_W_tmp[0] = b_R3_W_tmp;
  c_R3_W_tmp[3] = -R3_W_tmp;
  c_R3_W_tmp[6] = 0.0;
  c_R3_W_tmp[1] = R3_W_tmp;
  c_R3_W_tmp[4] = b_R3_W_tmp;
  c_R3_W_tmp[7] = 0.0;
  c_R3_W_tmp[2] = 0.0;
  Q_pqw_to_ijk[0] = 1.0;
  c_R3_W_tmp[5] = 0.0;
  Q_pqw_to_ijk[3] = 0.0;
  c_R3_W_tmp[8] = 1.0;
  Q_pqw_to_ijk[6] = 0.0;
  Q_pqw_to_ijk[1] = 0.0;
  Q_pqw_to_ijk[4] = b_R1_i_tmp;
  Q_pqw_to_ijk[7] = -R1_i_tmp;
  Q_pqw_to_ijk[2] = 0.0;
  Q_pqw_to_ijk[5] = R1_i_tmp;
  Q_pqw_to_ijk[8] = b_R1_i_tmp;
  for (i = 0; i < 3; i++) {
    b_R1_i_tmp = c_R3_W_tmp[i];
    d = c_R3_W_tmp[i + 3];
    i1 = (int)c_R3_W_tmp[i + 6];
    for (i2 = 0; i2 < 3; i2++) {
      d_R3_W_tmp[i + 3 * i2] =
          (b_R1_i_tmp * Q_pqw_to_ijk[3 * i2] + d * Q_pqw_to_ijk[3 * i2 + 1]) +
          (double)i1 * Q_pqw_to_ijk[3 * i2 + 2];
    }
  }
  c_R3_W_tmp[0] = b_R3_w_tmp;
  c_R3_W_tmp[3] = -R3_w_tmp;
  c_R3_W_tmp[6] = 0.0;
  c_R3_W_tmp[1] = R3_w_tmp;
  c_R3_W_tmp[4] = b_R3_w_tmp;
  c_R3_W_tmp[7] = 0.0;
  c_R3_W_tmp[2] = 0.0;
  c_R3_W_tmp[5] = 0.0;
  c_R3_W_tmp[8] = 1.0;
  /*  --- Step 3: Rotate r and v into IJK frame */
  c_a[0] = a * r_pqw_tmp;
  c_a[1] = a * b_r_pqw_tmp;
  c_a[2] = a * 0.0;
  for (i = 0; i < 3; i++) {
    b_R1_i_tmp = d_R3_W_tmp[i];
    d = d_R3_W_tmp[i + 3];
    R3_W_tmp = d_R3_W_tmp[i + 6];
    b_R3_W_tmp = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      R1_i_tmp =
          (b_R1_i_tmp * c_R3_W_tmp[3 * i1] + d * c_R3_W_tmp[3 * i1 + 1]) +
          R3_W_tmp * c_R3_W_tmp[3 * i1 + 2];
      Q_pqw_to_ijk[i + 3 * i1] = R1_i_tmp;
      b_R3_W_tmp += R1_i_tmp * c_a[i1];
    }
    r_ijk[i] = b_R3_W_tmp;
  }
  b_R1_i_tmp = b_a * -b_r_pqw_tmp;
  d = b_a * (ecc + r_pqw_tmp);
  R3_W_tmp = b_a * 0.0;
  for (i = 0; i < 3; i++) {
    v_ijk[i] = (Q_pqw_to_ijk[i] * b_R1_i_tmp + Q_pqw_to_ijk[i + 3] * d) +
               Q_pqw_to_ijk[i + 6] * R3_W_tmp;
  }
}

/*
 * File trailer for coe2rv.c
 *
 * [EOF]
 */
