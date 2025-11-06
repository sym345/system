/*
 * File: propagate_formation.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

/* Include Files */
#include "propagate_formation.h"
#include "coe2rv.h"
#include "formation_design_total_internal_types.h"
#include "formation_design_total_rtwutil.h"
#include "formation_design_total_types.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Declarations */
static double rt_atan2d_snf(double u0, double u1);

static double solve_kepler_equation(double M, double e);

/* Function Definitions */
/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int i;
    int i1;
    if (u0 > 0.0) {
      i = 1;
    } else {
      i = -1;
    }
    if (u1 > 0.0) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    y = atan2(i, i1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }
  return y;
}

/*
 * SOLVE_KEPLER_EQUATION 使用牛顿-拉夫逊方法求解开普勒方程
 *
 *  输入：
 *    M   - 平近点角 (rad)
 *    e   - 偏心率
 *    tol - 容差 (可选，默认1e-12)
 *
 *  输出：
 *    E   - 偏近点角 (rad)
 *
 * Arguments    : double M
 *                double e
 * Return Type  : double
 */
static double solve_kepler_equation(double M, double e)
{
  double E;
  double delta_E;
  int iter;
  boolean_T exitg1;
  /*  辅助函数：求解开普勒方程 */
  /*  归一化平近点角到 [0, 2π] */
  delta_E = M;
  if (rtIsNaN(M) || rtIsInf(M)) {
    M = rtNaN;
  } else if (M == 0.0) {
    M = 0.0;
  } else {
    boolean_T rEQ0;
    M = fmod(M, 6.2831853071795862);
    rEQ0 = (M == 0.0);
    if (!rEQ0) {
      double q;
      q = fabs(delta_E / 6.2831853071795862);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }
    if (rEQ0) {
      M = 0.0;
    } else if (delta_E < 0.0) {
      M += 6.2831853071795862;
    }
  }
  /*  初始猜测 */
  if (e < 0.8) {
    E = M;
  } else {
    E = 3.1415926535897931;
  }
  /*  牛顿-拉夫逊迭代 */
  iter = 0;
  exitg1 = false;
  while ((!exitg1) && (iter < 50)) {
    delta_E = -((E - e * sin(E)) - M) / (1.0 - e * cos(E));
    E += delta_E;
    if (fabs(delta_E) < 1.0E-12) {
      exitg1 = true;
    } else {
      iter++;
    }
  }
  return E;
}

/*
 * PROPAGATE_FORMATION 编队轨道传播函数
 *
 *  输入参数：
 *    satellites_coe    - 所有卫星轨道根数 (从formation_design函数输出)
 *    epoch_jd          - 初始历元时间 (儒略日)
 *    target_time_jd    - 目标时间 (儒略日)
 *    formation_config  - 编队构型参数结构体
 *
 *  输出参数：
 *    satellites_lvlh_t   - 目标时间的LVLH坐标系状态
 *    satellites_j2000_t  - 目标时间的J2000坐标系状态
 *
 * Arguments    : const struct2_T satellites_coe[11]
 *                double epoch_jd
 *                double target_time_jd
 *                double formation_config_p
 *                double formation_config_s
 *                double formation_config_l
 *                double formation_config_alpha
 *                double formation_config_theta
 *                struct_T satellites_lvlh_t[11]
 *                b_struct_T satellites_j2000_t[11]
 * Return Type  : void
 */
void propagate_formation(const struct2_T satellites_coe[11], double epoch_jd,
                         double target_time_jd, double formation_config_p,
                         double formation_config_s, double formation_config_l,
                         double formation_config_alpha,
                         double formation_config_theta,
                         struct_T satellites_lvlh_t[11],
                         b_struct_T satellites_j2000_t[11])
{
  double E_t;
  double E_t_tmp;
  double alpha;
  double dt_sec;
  double n_tmp;
  double theta;
  double x_lvlh_tmp;
  double y_lvlh_tmp;
  /*  常量定义 */
  /*  地球引力常数 (m^3/s^2) */
  /*  每天的秒数 */
  /*  计算时间差 */
  dt_sec = (target_time_jd - epoch_jd) * 86400.0;
  /*  主星轨道参数 */
  n_tmp = sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[0].a, 3.0));
  /*  平均角速度 */
  /*  编队构型参数 */
  alpha = 0.017453292519943295 * formation_config_alpha;
  theta = 0.017453292519943295 * formation_config_theta;
  /*  初始化输出 */
  /*  传播每颗卫星 */
  /*  主星：简单的开普勒传播 */
  /*  计算目标时间的平近点角 */
  /*  求解偏近点角 (牛顿-拉夫逊方法) */
  E_t_tmp = n_tmp * dt_sec;
  E_t = solve_kepler_equation(0.017453292519943295 * satellites_coe[0].M +
                                  E_t_tmp,
                              satellites_coe[0].e);
  /*  计算真近点角 */
  /*  转换为位置速度 */
  coe2rv(satellites_coe[0].a *
             (1.0 - satellites_coe[0].e * satellites_coe[0].e),
         satellites_coe[0].e, 0.017453292519943295 * satellites_coe[0].i,
         0.017453292519943295 * satellites_coe[0].Omega,
         0.017453292519943295 * satellites_coe[0].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[0].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[0].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  /*  主星在LVLH中的位置为原点 */
  satellites_lvlh_t[0].x = 0.0;
  satellites_lvlh_t[0].y = 0.0;
  satellites_lvlh_t[0].z = 0.0;
  satellites_j2000_t[0].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[0].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[0].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[0].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[0].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[0].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + theta;
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[1].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[1].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[1].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[1].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[1].a, 3.0)) *
              dt_sec,
      satellites_coe[1].e);
  coe2rv(satellites_coe[1].a *
             (1.0 - satellites_coe[1].e * satellites_coe[1].e),
         satellites_coe[1].e, 0.017453292519943295 * satellites_coe[1].i,
         0.017453292519943295 * satellites_coe[1].Omega,
         0.017453292519943295 * satellites_coe[1].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[1].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[1].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[1].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[1].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[1].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[1].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[1].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[1].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 0.62831853071795862);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[2].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[2].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[2].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[2].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[2].a, 3.0)) *
              dt_sec,
      satellites_coe[2].e);
  coe2rv(satellites_coe[2].a *
             (1.0 - satellites_coe[2].e * satellites_coe[2].e),
         satellites_coe[2].e, 0.017453292519943295 * satellites_coe[2].i,
         0.017453292519943295 * satellites_coe[2].Omega,
         0.017453292519943295 * satellites_coe[2].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[2].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[2].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[2].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[2].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[2].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[2].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[2].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[2].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 1.2566370614359172);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[3].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[3].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[3].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[3].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[3].a, 3.0)) *
              dt_sec,
      satellites_coe[3].e);
  coe2rv(satellites_coe[3].a *
             (1.0 - satellites_coe[3].e * satellites_coe[3].e),
         satellites_coe[3].e, 0.017453292519943295 * satellites_coe[3].i,
         0.017453292519943295 * satellites_coe[3].Omega,
         0.017453292519943295 * satellites_coe[3].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[3].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[3].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[3].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[3].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[3].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[3].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[3].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[3].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 1.8849555921538759);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[4].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[4].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[4].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[4].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[4].a, 3.0)) *
              dt_sec,
      satellites_coe[4].e);
  coe2rv(satellites_coe[4].a *
             (1.0 - satellites_coe[4].e * satellites_coe[4].e),
         satellites_coe[4].e, 0.017453292519943295 * satellites_coe[4].i,
         0.017453292519943295 * satellites_coe[4].Omega,
         0.017453292519943295 * satellites_coe[4].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[4].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[4].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[4].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[4].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[4].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[4].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[4].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[4].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 2.5132741228718345);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[5].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[5].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[5].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[5].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[5].a, 3.0)) *
              dt_sec,
      satellites_coe[5].e);
  coe2rv(satellites_coe[5].a *
             (1.0 - satellites_coe[5].e * satellites_coe[5].e),
         satellites_coe[5].e, 0.017453292519943295 * satellites_coe[5].i,
         0.017453292519943295 * satellites_coe[5].Omega,
         0.017453292519943295 * satellites_coe[5].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[5].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[5].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[5].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[5].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[5].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[5].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[5].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[5].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 3.1415926535897931);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[6].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[6].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[6].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[6].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[6].a, 3.0)) *
              dt_sec,
      satellites_coe[6].e);
  coe2rv(satellites_coe[6].a *
             (1.0 - satellites_coe[6].e * satellites_coe[6].e),
         satellites_coe[6].e, 0.017453292519943295 * satellites_coe[6].i,
         0.017453292519943295 * satellites_coe[6].Omega,
         0.017453292519943295 * satellites_coe[6].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[6].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[6].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[6].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[6].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[6].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[6].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[6].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[6].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 3.7699111843077517);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[7].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[7].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[7].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[7].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[7].a, 3.0)) *
              dt_sec,
      satellites_coe[7].e);
  coe2rv(satellites_coe[7].a *
             (1.0 - satellites_coe[7].e * satellites_coe[7].e),
         satellites_coe[7].e, 0.017453292519943295 * satellites_coe[7].i,
         0.017453292519943295 * satellites_coe[7].Omega,
         0.017453292519943295 * satellites_coe[7].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[7].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[7].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[7].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[7].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[7].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[7].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[7].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[7].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 4.39822971502571);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[8].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[8].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[8].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[8].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[8].a, 3.0)) *
              dt_sec,
      satellites_coe[8].e);
  coe2rv(satellites_coe[8].a *
             (1.0 - satellites_coe[8].e * satellites_coe[8].e),
         satellites_coe[8].e, 0.017453292519943295 * satellites_coe[8].i,
         0.017453292519943295 * satellites_coe[8].Omega,
         0.017453292519943295 * satellites_coe[8].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[8].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[8].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[8].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[8].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[8].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[8].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[8].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[8].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 5.026548245743669);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[9].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[9].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[9].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[9].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[9].a, 3.0)) *
              dt_sec,
      satellites_coe[9].e);
  coe2rv(satellites_coe[9].a *
             (1.0 - satellites_coe[9].e * satellites_coe[9].e),
         satellites_coe[9].e, 0.017453292519943295 * satellites_coe[9].i,
         0.017453292519943295 * satellites_coe[9].Omega,
         0.017453292519943295 * satellites_coe[9].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[9].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[9].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
  satellites_j2000_t[9].r[0] = satellites_j2000_t[10].r[0];
  satellites_j2000_t[9].v[0] = satellites_j2000_t[10].v[0];
  satellites_j2000_t[9].r[1] = satellites_j2000_t[10].r[1];
  satellites_j2000_t[9].v[1] = satellites_j2000_t[10].v[1];
  satellites_j2000_t[9].r[2] = satellites_j2000_t[10].r[2];
  satellites_j2000_t[9].v[2] = satellites_j2000_t[10].v[2];
  /*  子星：使用相对运动方程 */
  /*  为每颗子星分配相位 */
  /*  当前时刻的相对位置 (LVLH坐标系) */
  E_t = E_t_tmp + (theta + 5.6548667764616276);
  x_lvlh_tmp = cos(E_t);
  satellites_lvlh_t[10].x = -formation_config_p * x_lvlh_tmp;
  y_lvlh_tmp = sin(E_t);
  satellites_lvlh_t[10].y =
      2.0 * formation_config_p * y_lvlh_tmp + formation_config_l;
  E_t -= alpha;
  satellites_lvlh_t[10].z = formation_config_s * sin(E_t);
  /*  相对速度 (LVLH坐标系) */
  /*  子星在J2000坐标系中的状态需要通过坐标变换得到 */
  /*  这里使用简化的开普勒传播 */
  E_t = solve_kepler_equation(
      0.017453292519943295 * satellites_coe[10].M +
          sqrt(3.986004418E+14 / rt_powd_snf(satellites_coe[10].a, 3.0)) *
              dt_sec,
      satellites_coe[10].e);
  coe2rv(satellites_coe[10].a *
             (1.0 - satellites_coe[10].e * satellites_coe[10].e),
         satellites_coe[10].e, 0.017453292519943295 * satellites_coe[10].i,
         0.017453292519943295 * satellites_coe[10].Omega,
         0.017453292519943295 * satellites_coe[10].omega,
         2.0 * rt_atan2d_snf(sqrt(satellites_coe[10].e + 1.0) * sin(E_t / 2.0),
                             sqrt(1.0 - satellites_coe[10].e) * cos(E_t / 2.0)),
         satellites_j2000_t[10].r, satellites_j2000_t[10].v);
}

/*
 * File trailer for propagate_formation.c
 *
 * [EOF]
 */
