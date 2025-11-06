/*
 * File: formation_design_total.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

/* Include Files */
#include "formation_design_total.h"
#include "coe2rv.h"
#include "formation_design_total_emxutil.h"
#include "formation_design_total_internal_types.h"
#include "formation_design_total_rtwutil.h"
#include "formation_design_total_types.h"
#include "propagate_formation.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Type Definitions */
#ifndef typedef_cell_wrap_0
#define typedef_cell_wrap_0
typedef struct {
  emxArray_char_T_1x24 f1;
} cell_wrap_0;
#endif /* typedef_cell_wrap_0 */

/* Function Definitions */
/*
 * 输入参数：
 *    chief_coe        - 主星轨道根数
 *    epoch_jd         - 初始历元时间 (儒略日)
 *    formation_config - 编队构型参数
 *    num_satellites   - 卫星数量
 *    formation_type   - 编队类型
 *    duration_hours   - 演化时长 (小时)
 *
 * Arguments    : const struct0_T *chief_coe
 *                double epoch_jd
 *                const struct1_T *formation_config
 *                double num_satellites
 *                double formation_type
 *                double duration_hours
 *                struct2_T satellites_coe[11]
 *                struct3_T satellites_lvlh[11]
 *                struct4_T satellites_j2000[11]
 *                struct5_T *formation_params
 *                struct6_T *trajectory
 * Return Type  : void
 */
void formation_design_total(const struct0_T *chief_coe, double epoch_jd,
                            const struct1_T *formation_config,
                            double num_satellites, double formation_type,
                            double duration_hours, struct2_T satellites_coe[11],
                            struct3_T satellites_lvlh[11],
                            struct4_T satellites_j2000[11],
                            struct5_T *formation_params, struct6_T *trajectory)
{
  static const char cv2[24] = {'S', 'P', 'A', 'C', 'E', '_', 'C', 'I',
                               'R', 'C', 'U', 'L', 'A', 'R', '_', 'F',
                               'O', 'R', 'M', 'A', 'T', 'I', 'O', 'N'};
  static const char cv1[16] = {'D', 'U', 'A', 'L', '_', 'L', 'A', 'Y',
                               'E', 'R', '_', 'O', 'R', 'B', 'I', 'T'};
  static const char cv[14] = {'C', 'O', 'P', 'L', 'A', 'N', 'A',
                              'R', '_', 'O', 'R', 'B', 'I', 'T'};
  static const char cv3[12] = {'U', 'N', 'K', 'N', 'O', 'W',
                               'N', '_', 'T', 'Y', 'P', 'E'};
  cell_wrap_0 r;
  cell_wrap_0 r1;
  cell_wrap_0 r2;
  emxArray_real_T *phase_inner;
  emxArray_real_T *phase_outer;
  emxArray_real_T *y;
  struct2_T expl_temp;
  struct3_T b_expl_temp;
  double r_chief_j2000[3];
  double v_chief_j2000[3];
  double Omega_ref;
  double T;
  double a_ref;
  double alpha;
  double d;
  double i_ref;
  double idx_outer;
  double n;
  double num_inner;
  double num_outer;
  double omega_ref;
  double p_inner;
  double p_outer;
  double s;
  double theta;
  double *phase_inner_data;
  double *phase_outer_data;
  double *y_data;
  int i;
  int loop_ub;
  int sat_idx;
  /*  输出参数： */
  /*    satellites_coe    - 所有卫星轨道根数 */
  /*    satellites_lvlh   - LVLH坐标系下状态 */
  /*    satellites_j2000  - J2000坐标系下状态 */
  /*    formation_params  - 编队参数 */
  /*    trajectory        - 轨迹数据结构（X/Y/Z） */
  /*  ====== 初始编队设计 ====== */
  /*  常量定义 */
  /*  地球引力常数 (m^3/s^2) */
  /*  输入参数转换 */
  a_ref = chief_coe->a;
  i_ref = 0.017453292519943295 * chief_coe->i;
  Omega_ref = 0.017453292519943295 * chief_coe->Omega;
  omega_ref = 0.017453292519943295 * chief_coe->omega;
  /*  编队构型参数 */
  s = formation_config->s;
  alpha = 0.017453292519943295 * formation_config->alpha;
  theta = 0.017453292519943295 * formation_config->theta;
  /*  计算轨道参数 */
  n = sqrt(3.986004418E+14 / rt_powd_snf(chief_coe->a, 3.0));
  /*  平均角速度 */
  T = 6.2831853071795862 / n;
  /*  轨道周期 */
  /*  初始化变量 */
  num_inner = 0.0;
  num_outer = 0.0;
  p_inner = 0.0;
  p_outer = 0.0;
  emxInit_real_T(&phase_inner);
  phase_inner_data = phase_inner->data;
  phase_inner->size[0] = 0;
  phase_inner->size[1] = 0;
  emxInit_real_T(&phase_outer);
  phase_outer_data = phase_outer->data;
  phase_outer->size[0] = 0;
  phase_outer->size[1] = 0;
  /*  根据编队类型调整参数 */
  if ((int)formation_type == 2) {
    /*  双层绕飞 */
    p_inner = formation_config->p_inner;
    p_outer = formation_config->p_outer;
    num_inner = formation_config->num_inner;
    num_outer = formation_config->num_outer;
    /*  内层相位 */
    emxInit_real_T(&y);
    y_data = y->data;
    if (rtIsNaN(formation_config->num_inner - 1.0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      y_data[0] = rtNaN;
    } else if (formation_config->num_inner - 1.0 < 0.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)(formation_config->num_inner - 1.0) + 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      loop_ub = (int)(formation_config->num_inner - 1.0);
      for (i = 0; i <= loop_ub; i++) {
        y_data[i] = i;
      }
    }
    idx_outer = 6.2831853071795862 / formation_config->num_inner;
    i = phase_inner->size[0] * phase_inner->size[1];
    phase_inner->size[0] = 1;
    phase_inner->size[1] = y->size[1];
    emxEnsureCapacity_real_T(phase_inner, i);
    phase_inner_data = phase_inner->data;
    loop_ub = y->size[1];
    for (i = 0; i < loop_ub; i++) {
      phase_inner_data[i] = y_data[i] * idx_outer;
    }
    /*  外层相位 */
    if (rtIsNaN(formation_config->num_outer - 1.0)) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      y_data[0] = rtNaN;
    } else if (formation_config->num_outer - 1.0 < 0.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)(formation_config->num_outer - 1.0) + 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      loop_ub = (int)(formation_config->num_outer - 1.0);
      for (i = 0; i <= loop_ub; i++) {
        y_data[i] = i;
      }
    }
    idx_outer = 6.2831853071795862 / formation_config->num_outer;
    i = phase_outer->size[0] * phase_outer->size[1];
    phase_outer->size[0] = 1;
    phase_outer->size[1] = y->size[1];
    emxEnsureCapacity_real_T(phase_outer, i);
    phase_outer_data = phase_outer->data;
    d = 3.1415926535897931 / formation_config->num_outer;
    loop_ub = y->size[1];
    for (i = 0; i < loop_ub; i++) {
      phase_outer_data[i] = y_data[i] * idx_outer + d;
    }
    emxFree_real_T(&y);
  }
  /*  主星在J2000坐标系中的状态 */
  d = chief_coe->e * chief_coe->e;
  coe2rv(chief_coe->a * (1.0 - d), chief_coe->e, i_ref, Omega_ref, omega_ref,
         0.0, r_chief_j2000, v_chief_j2000);
  /*  初始化输出 */
  /* 要固定住数组大小 */
  /*  默认结构体 */
  /*  初始化所有元素 */
  satellites_coe[1].sat_id = -1.0;
  satellites_coe[1].a = 0.0;
  satellites_coe[1].e = 0.0;
  satellites_coe[1].i = 0.0;
  satellites_coe[1].Omega = 0.0;
  satellites_coe[1].omega = 0.0;
  satellites_coe[1].M = 0.0;
  satellites_coe[1].epoch_jd = epoch_jd;
  satellites_lvlh[1].sat_id = -1.0;
  satellites_lvlh[1].x = 0.0;
  satellites_lvlh[1].y = 0.0;
  satellites_lvlh[1].z = 0.0;
  satellites_lvlh[1].vx = 0.0;
  satellites_lvlh[1].vy = 0.0;
  satellites_lvlh[1].vz = 0.0;
  satellites_lvlh[1].epoch_jd = epoch_jd;
  satellites_j2000[1].sat_id = -1.0;
  satellites_j2000[1].r.size[0] = 3;
  satellites_j2000[1].r.size[1] = 1;
  satellites_j2000[1].v.size[0] = 3;
  satellites_j2000[1].v.size[1] = 1;
  satellites_j2000[1].epoch_jd = epoch_jd;
  satellites_coe[2].sat_id = -1.0;
  satellites_coe[2].a = 0.0;
  satellites_coe[2].e = 0.0;
  satellites_coe[2].i = 0.0;
  satellites_coe[2].Omega = 0.0;
  satellites_coe[2].omega = 0.0;
  satellites_coe[2].M = 0.0;
  satellites_coe[2].epoch_jd = epoch_jd;
  satellites_lvlh[2].sat_id = -1.0;
  satellites_lvlh[2].x = 0.0;
  satellites_lvlh[2].y = 0.0;
  satellites_lvlh[2].z = 0.0;
  satellites_lvlh[2].vx = 0.0;
  satellites_lvlh[2].vy = 0.0;
  satellites_lvlh[2].vz = 0.0;
  satellites_lvlh[2].epoch_jd = epoch_jd;
  satellites_j2000[2].sat_id = -1.0;
  satellites_j2000[2].r.size[0] = 3;
  satellites_j2000[2].r.size[1] = 1;
  satellites_j2000[2].v.size[0] = 3;
  satellites_j2000[2].v.size[1] = 1;
  satellites_j2000[2].epoch_jd = epoch_jd;
  satellites_coe[3].sat_id = -1.0;
  satellites_coe[3].a = 0.0;
  satellites_coe[3].e = 0.0;
  satellites_coe[3].i = 0.0;
  satellites_coe[3].Omega = 0.0;
  satellites_coe[3].omega = 0.0;
  satellites_coe[3].M = 0.0;
  satellites_coe[3].epoch_jd = epoch_jd;
  satellites_lvlh[3].sat_id = -1.0;
  satellites_lvlh[3].x = 0.0;
  satellites_lvlh[3].y = 0.0;
  satellites_lvlh[3].z = 0.0;
  satellites_lvlh[3].vx = 0.0;
  satellites_lvlh[3].vy = 0.0;
  satellites_lvlh[3].vz = 0.0;
  satellites_lvlh[3].epoch_jd = epoch_jd;
  satellites_j2000[3].sat_id = -1.0;
  satellites_j2000[3].r.size[0] = 3;
  satellites_j2000[3].r.size[1] = 1;
  satellites_j2000[3].v.size[0] = 3;
  satellites_j2000[3].v.size[1] = 1;
  satellites_j2000[3].epoch_jd = epoch_jd;
  satellites_coe[4].sat_id = -1.0;
  satellites_coe[4].a = 0.0;
  satellites_coe[4].e = 0.0;
  satellites_coe[4].i = 0.0;
  satellites_coe[4].Omega = 0.0;
  satellites_coe[4].omega = 0.0;
  satellites_coe[4].M = 0.0;
  satellites_coe[4].epoch_jd = epoch_jd;
  satellites_lvlh[4].sat_id = -1.0;
  satellites_lvlh[4].x = 0.0;
  satellites_lvlh[4].y = 0.0;
  satellites_lvlh[4].z = 0.0;
  satellites_lvlh[4].vx = 0.0;
  satellites_lvlh[4].vy = 0.0;
  satellites_lvlh[4].vz = 0.0;
  satellites_lvlh[4].epoch_jd = epoch_jd;
  satellites_j2000[4].sat_id = -1.0;
  satellites_j2000[4].r.size[0] = 3;
  satellites_j2000[4].r.size[1] = 1;
  satellites_j2000[4].v.size[0] = 3;
  satellites_j2000[4].v.size[1] = 1;
  satellites_j2000[4].epoch_jd = epoch_jd;
  satellites_coe[5].sat_id = -1.0;
  satellites_coe[5].a = 0.0;
  satellites_coe[5].e = 0.0;
  satellites_coe[5].i = 0.0;
  satellites_coe[5].Omega = 0.0;
  satellites_coe[5].omega = 0.0;
  satellites_coe[5].M = 0.0;
  satellites_coe[5].epoch_jd = epoch_jd;
  satellites_lvlh[5].sat_id = -1.0;
  satellites_lvlh[5].x = 0.0;
  satellites_lvlh[5].y = 0.0;
  satellites_lvlh[5].z = 0.0;
  satellites_lvlh[5].vx = 0.0;
  satellites_lvlh[5].vy = 0.0;
  satellites_lvlh[5].vz = 0.0;
  satellites_lvlh[5].epoch_jd = epoch_jd;
  satellites_j2000[5].sat_id = -1.0;
  satellites_j2000[5].r.size[0] = 3;
  satellites_j2000[5].r.size[1] = 1;
  satellites_j2000[5].v.size[0] = 3;
  satellites_j2000[5].v.size[1] = 1;
  satellites_j2000[5].epoch_jd = epoch_jd;
  satellites_coe[6].sat_id = -1.0;
  satellites_coe[6].a = 0.0;
  satellites_coe[6].e = 0.0;
  satellites_coe[6].i = 0.0;
  satellites_coe[6].Omega = 0.0;
  satellites_coe[6].omega = 0.0;
  satellites_coe[6].M = 0.0;
  satellites_coe[6].epoch_jd = epoch_jd;
  satellites_lvlh[6].sat_id = -1.0;
  satellites_lvlh[6].x = 0.0;
  satellites_lvlh[6].y = 0.0;
  satellites_lvlh[6].z = 0.0;
  satellites_lvlh[6].vx = 0.0;
  satellites_lvlh[6].vy = 0.0;
  satellites_lvlh[6].vz = 0.0;
  satellites_lvlh[6].epoch_jd = epoch_jd;
  satellites_j2000[6].sat_id = -1.0;
  satellites_j2000[6].r.size[0] = 3;
  satellites_j2000[6].r.size[1] = 1;
  satellites_j2000[6].v.size[0] = 3;
  satellites_j2000[6].v.size[1] = 1;
  satellites_j2000[6].epoch_jd = epoch_jd;
  satellites_coe[7].sat_id = -1.0;
  satellites_coe[7].a = 0.0;
  satellites_coe[7].e = 0.0;
  satellites_coe[7].i = 0.0;
  satellites_coe[7].Omega = 0.0;
  satellites_coe[7].omega = 0.0;
  satellites_coe[7].M = 0.0;
  satellites_coe[7].epoch_jd = epoch_jd;
  satellites_lvlh[7].sat_id = -1.0;
  satellites_lvlh[7].x = 0.0;
  satellites_lvlh[7].y = 0.0;
  satellites_lvlh[7].z = 0.0;
  satellites_lvlh[7].vx = 0.0;
  satellites_lvlh[7].vy = 0.0;
  satellites_lvlh[7].vz = 0.0;
  satellites_lvlh[7].epoch_jd = epoch_jd;
  satellites_j2000[7].sat_id = -1.0;
  satellites_j2000[7].r.size[0] = 3;
  satellites_j2000[7].r.size[1] = 1;
  satellites_j2000[7].v.size[0] = 3;
  satellites_j2000[7].v.size[1] = 1;
  satellites_j2000[7].epoch_jd = epoch_jd;
  satellites_coe[8].sat_id = -1.0;
  satellites_coe[8].a = 0.0;
  satellites_coe[8].e = 0.0;
  satellites_coe[8].i = 0.0;
  satellites_coe[8].Omega = 0.0;
  satellites_coe[8].omega = 0.0;
  satellites_coe[8].M = 0.0;
  satellites_coe[8].epoch_jd = epoch_jd;
  satellites_lvlh[8].sat_id = -1.0;
  satellites_lvlh[8].x = 0.0;
  satellites_lvlh[8].y = 0.0;
  satellites_lvlh[8].z = 0.0;
  satellites_lvlh[8].vx = 0.0;
  satellites_lvlh[8].vy = 0.0;
  satellites_lvlh[8].vz = 0.0;
  satellites_lvlh[8].epoch_jd = epoch_jd;
  satellites_j2000[8].sat_id = -1.0;
  satellites_j2000[8].r.size[0] = 3;
  satellites_j2000[8].r.size[1] = 1;
  satellites_j2000[8].v.size[0] = 3;
  satellites_j2000[8].v.size[1] = 1;
  satellites_j2000[8].epoch_jd = epoch_jd;
  satellites_coe[9].sat_id = -1.0;
  satellites_coe[9].a = 0.0;
  satellites_coe[9].e = 0.0;
  satellites_coe[9].i = 0.0;
  satellites_coe[9].Omega = 0.0;
  satellites_coe[9].omega = 0.0;
  satellites_coe[9].M = 0.0;
  satellites_coe[9].epoch_jd = epoch_jd;
  satellites_lvlh[9].sat_id = -1.0;
  satellites_lvlh[9].x = 0.0;
  satellites_lvlh[9].y = 0.0;
  satellites_lvlh[9].z = 0.0;
  satellites_lvlh[9].vx = 0.0;
  satellites_lvlh[9].vy = 0.0;
  satellites_lvlh[9].vz = 0.0;
  satellites_lvlh[9].epoch_jd = epoch_jd;
  satellites_j2000[9].sat_id = -1.0;
  satellites_j2000[9].r.size[0] = 3;
  satellites_j2000[9].r.size[1] = 1;
  satellites_j2000[9].v.size[0] = 3;
  satellites_j2000[9].v.size[1] = 1;
  satellites_j2000[9].epoch_jd = epoch_jd;
  satellites_coe[10].sat_id = -1.0;
  satellites_coe[10].a = 0.0;
  satellites_coe[10].e = 0.0;
  satellites_coe[10].i = 0.0;
  satellites_coe[10].Omega = 0.0;
  satellites_coe[10].omega = 0.0;
  satellites_coe[10].M = 0.0;
  satellites_coe[10].epoch_jd = epoch_jd;
  satellites_lvlh[10].sat_id = -1.0;
  satellites_lvlh[10].x = 0.0;
  satellites_lvlh[10].y = 0.0;
  satellites_lvlh[10].z = 0.0;
  satellites_lvlh[10].vx = 0.0;
  satellites_lvlh[10].vy = 0.0;
  satellites_lvlh[10].vz = 0.0;
  satellites_lvlh[10].epoch_jd = epoch_jd;
  satellites_j2000[10].sat_id = -1.0;
  satellites_j2000[10].r.size[0] = 3;
  satellites_j2000[10].r.size[1] = 1;
  satellites_j2000[10].v.size[0] = 3;
  satellites_j2000[10].v.size[1] = 1;
  for (i = 0; i < 3; i++) {
    satellites_j2000[1].r.data[i] = 0.0;
    satellites_j2000[1].v.data[i] = 0.0;
    satellites_j2000[2].r.data[i] = 0.0;
    satellites_j2000[2].v.data[i] = 0.0;
    satellites_j2000[3].r.data[i] = 0.0;
    satellites_j2000[3].v.data[i] = 0.0;
    satellites_j2000[4].r.data[i] = 0.0;
    satellites_j2000[4].v.data[i] = 0.0;
    satellites_j2000[5].r.data[i] = 0.0;
    satellites_j2000[5].v.data[i] = 0.0;
    satellites_j2000[6].r.data[i] = 0.0;
    satellites_j2000[6].v.data[i] = 0.0;
    satellites_j2000[7].r.data[i] = 0.0;
    satellites_j2000[7].v.data[i] = 0.0;
    satellites_j2000[8].r.data[i] = 0.0;
    satellites_j2000[8].v.data[i] = 0.0;
    satellites_j2000[9].r.data[i] = 0.0;
    satellites_j2000[9].v.data[i] = 0.0;
    satellites_j2000[10].r.data[i] = 0.0;
    satellites_j2000[10].v.data[i] = 0.0;
  }
  satellites_j2000[10].epoch_jd = epoch_jd;
  /*  主星数据 */
  satellites_coe[0].sat_id = 0.0;
  satellites_coe[0].a = chief_coe->a;
  satellites_coe[0].e = chief_coe->e;
  satellites_coe[0].i = chief_coe->i;
  satellites_coe[0].Omega = chief_coe->Omega;
  satellites_coe[0].omega = chief_coe->omega;
  satellites_coe[0].M = chief_coe->M;
  satellites_coe[0].epoch_jd = epoch_jd;
  satellites_lvlh[0].sat_id = 0.0;
  satellites_lvlh[0].x = 0.0;
  satellites_lvlh[0].y = 0.0;
  satellites_lvlh[0].z = 0.0;
  satellites_lvlh[0].vx = 0.0;
  satellites_lvlh[0].vy = 0.0;
  satellites_lvlh[0].vz = 0.0;
  satellites_lvlh[0].epoch_jd = epoch_jd;
  satellites_j2000[0].sat_id = 0.0;
  satellites_j2000[0].r.size[0] = 1;
  satellites_j2000[0].r.size[1] = 3;
  satellites_j2000[0].v.size[0] = 1;
  satellites_j2000[0].v.size[1] = 3;
  satellites_j2000[0].r.data[0] = r_chief_j2000[0];
  satellites_j2000[0].v.data[0] = v_chief_j2000[0];
  satellites_j2000[0].r.data[1] = r_chief_j2000[1];
  satellites_j2000[0].v.data[1] = v_chief_j2000[1];
  satellites_j2000[0].r.data[2] = r_chief_j2000[2];
  satellites_j2000[0].v.data[2] = v_chief_j2000[2];
  satellites_j2000[0].epoch_jd = epoch_jd;
  /*  子星计算 */
  if (num_satellites == 0.0) {
    /*  如果没有子星，直接返回主星信息 */
    formation_params->chief_coe = *chief_coe;
    formation_params->formation_config = *formation_config;
    formation_params->num_satellites = 0.0;
    formation_params->formation_type = formation_type;
    formation_params->epoch_jd = epoch_jd;
    formation_params->orbital_period = T;
    formation_params->mean_motion = n;
    r.f1.size[0] = 1;
    r.f1.size[1] = 14;
    for (i = 0; i < 14; i++) {
      r.f1.data[i] = cv[i];
    }
    r1.f1.size[0] = 1;
    r1.f1.size[1] = 16;
    for (i = 0; i < 16; i++) {
      r1.f1.data[i] = cv1[i];
    }
    r2.f1.size[0] = 1;
    r2.f1.size[1] = 24;
    for (i = 0; i < 24; i++) {
      r2.f1.data[i] = cv2[i];
    }
    if ((formation_type >= 1.0) && (formation_type <= 3.0)) {
      cell_wrap_0 rv[3];
      cell_wrap_0 rv1[3];
      formation_params->type_name.size[0] = 1;
      rv[0] = r;
      rv[1] = r1;
      rv[2] = r2;
      formation_params->type_name.size[1] =
          rv[(int)formation_type - 1].f1.size[1];
      rv[0] = r;
      rv[1] = r1;
      rv[2] = r2;
      rv1[0] = r;
      rv1[1] = r1;
      rv1[2] = r2;
      loop_ub = rv1[(int)formation_type - 1].f1.size[1];
      for (i = 0; i < loop_ub; i++) {
        formation_params->type_name.data[i] =
            rv[(int)formation_type - 1].f1.data[i];
      }
    } else {
      formation_params->type_name.size[0] = 1;
      formation_params->type_name.size[1] = 12;
      for (i = 0; i < 12; i++) {
        formation_params->type_name.data[i] = cv3[i];
      }
    }
  } else {
    /*  子星计算 */
    i = (int)num_satellites;
    if ((int)num_satellites - 1 >= 0) {
      expl_temp.a = chief_coe->a;
      expl_temp.Omega = 57.295779513082323 * Omega_ref;
      expl_temp.omega = 57.295779513082323 * omega_ref;
      expl_temp.M = 57.295779513082323 * (0.017453292519943295 * chief_coe->M);
      expl_temp.epoch_jd = epoch_jd;
      b_expl_temp.epoch_jd = epoch_jd;
    }
    for (sat_idx = 0; sat_idx < i; sat_idx++) {
      double current_p;
      double current_s;
      double current_theta;
      double i_cir;
      double x_lvlh_tmp_tmp;
      double y_lvlh_tmp;
      /*  初始化当前卫星的参数 */
      current_p = formation_config->p;
      /*  默认值 */
      current_s = s;
      /*  默认值 */
      current_theta = theta;
      /*  默认值 */
      switch ((int)formation_type) {
      case 1:
        /*  共面绕飞 */
        current_s = 0.0;
        current_theta = theta + 6.2831853071795862 *
                                    (((double)sat_idx + 1.0) - 1.0) /
                                    num_satellites;
        break;
      case 2:
        /*  双层绕飞 */
        if ((double)sat_idx + 1.0 <= num_inner) {
          /*  内层卫星 */
          current_p = p_inner;
          current_theta = theta + phase_inner_data[sat_idx];
        } else if (((double)sat_idx + 1.0 > num_inner) && (num_outer > 0.0)) {
          /*  外层卫星 */
          idx_outer = ((double)sat_idx + 1.0) - num_inner;
          if ((phase_outer->size[0] == 0) || (phase_outer->size[1] == 0)) {
            loop_ub = 0;
          } else {
            loop_ub = phase_outer->size[1];
          }
          if (idx_outer <= loop_ub) {
            /*  边界检查 */
            current_p = p_outer;
            current_s = -s;
            current_theta = theta + phase_outer_data[(int)idx_outer - 1];
          }
        }
        break;
      case 3:
        /*  空间圆构型 */
        current_theta = theta + 6.2831853071795862 *
                                    (((double)sat_idx + 1.0) - 1.0) /
                                    num_satellites;
        break;
      }
      /*     %% LVLH位置速度 */
      x_lvlh_tmp_tmp = cos(current_theta);
      y_lvlh_tmp = sin(current_theta);
      /*     %% 子星轨道根 */
      idx_outer = current_p / a_ref;
      idx_outer = sqrt((idx_outer * idx_outer + d) +
                       2.0 * current_p / a_ref * chief_coe->e * x_lvlh_tmp_tmp);
      i_cir = i_ref + current_s / a_ref;
      coe2rv(a_ref * (1.0 - idx_outer * idx_outer), idx_outer, i_cir, Omega_ref,
             omega_ref, 0.0, r_chief_j2000, v_chief_j2000);
      /*     %% 存储 */
      expl_temp.sat_id = (double)sat_idx + 1.0;
      expl_temp.e = idx_outer;
      expl_temp.i = 57.295779513082323 * i_cir;
      satellites_coe[(int)(((double)sat_idx + 1.0) + 1.0) - 1] = expl_temp;
      b_expl_temp.sat_id = (double)sat_idx + 1.0;
      b_expl_temp.x = -current_p * x_lvlh_tmp_tmp;
      b_expl_temp.y = 2.0 * current_p * y_lvlh_tmp + formation_config->l;
      idx_outer = current_theta - alpha;
      b_expl_temp.z = current_s * sin(idx_outer);
      b_expl_temp.vx = current_p * n * y_lvlh_tmp;
      b_expl_temp.vy = 2.0 * current_p * n * x_lvlh_tmp_tmp;
      b_expl_temp.vz = current_s * n * cos(idx_outer);
      satellites_lvlh[(int)(((double)sat_idx + 1.0) + 1.0) - 1] = b_expl_temp;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].sat_id =
          (double)sat_idx + 1.0;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].r.size[0] = 1;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].r.size[1] = 3;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].v.size[0] = 1;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].v.size[1] = 3;
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].r.data[0] =
          r_chief_j2000[0];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].v.data[0] =
          v_chief_j2000[0];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].r.data[1] =
          r_chief_j2000[1];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].v.data[1] =
          v_chief_j2000[1];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].r.data[2] =
          r_chief_j2000[2];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].v.data[2] =
          v_chief_j2000[2];
      satellites_j2000[(int)(((double)sat_idx + 1.0) + 1.0) - 1].epoch_jd =
          epoch_jd;
    }
    /*  编队参数信息 */
    formation_params->chief_coe = *chief_coe;
    formation_params->formation_config = *formation_config;
    formation_params->num_satellites = num_satellites;
    formation_params->formation_type = formation_type;
    formation_params->epoch_jd = epoch_jd;
    formation_params->orbital_period = T;
    formation_params->mean_motion = n;
    r.f1.size[0] = 1;
    r.f1.size[1] = 14;
    for (i = 0; i < 14; i++) {
      r.f1.data[i] = cv[i];
    }
    r1.f1.size[0] = 1;
    r1.f1.size[1] = 16;
    for (i = 0; i < 16; i++) {
      r1.f1.data[i] = cv1[i];
    }
    r2.f1.size[0] = 1;
    r2.f1.size[1] = 24;
    for (i = 0; i < 24; i++) {
      r2.f1.data[i] = cv2[i];
    }
    if ((formation_type >= 1.0) && (formation_type <= 3.0)) {
      cell_wrap_0 rv[3];
      cell_wrap_0 rv1[3];
      formation_params->type_name.size[0] = 1;
      rv[0] = r;
      rv[1] = r1;
      rv[2] = r2;
      formation_params->type_name.size[1] =
          rv[(int)formation_type - 1].f1.size[1];
      rv[0] = r;
      rv[1] = r1;
      rv[2] = r2;
      rv1[0] = r;
      rv1[1] = r1;
      rv1[2] = r2;
      loop_ub = rv1[(int)formation_type - 1].f1.size[1];
      for (i = 0; i < loop_ub; i++) {
        formation_params->type_name.data[i] =
            rv[(int)formation_type - 1].f1.data[i];
      }
    } else {
      formation_params->type_name.size[0] = 1;
      formation_params->type_name.size[1] = 12;
      for (i = 0; i < 12; i++) {
        formation_params->type_name.data[i] = cv3[i];
      }
    }
  }
  emxFree_real_T(&phase_outer);
  emxFree_real_T(&phase_inner);
  /*  ====== 参数设置 ====== */
  /*  时间步数 */
  /*  创建时间序列 */
  trajectory->time_array_hours[99] = duration_hours;
  trajectory->time_array_hours[0] = 0.0;
  if (-duration_hours == 0.0) {
    idx_outer = duration_hours / 99.0;
    for (loop_ub = 0; loop_ub < 98; loop_ub++) {
      trajectory->time_array_hours[loop_ub + 1] =
          (2.0 * ((double)loop_ub + 2.0) - 101.0) * idx_outer;
    }
  } else if ((duration_hours < 0.0) &&
             (fabs(duration_hours) > 8.9884656743115785E+307)) {
    idx_outer = duration_hours / 99.0;
    for (loop_ub = 0; loop_ub < 98; loop_ub++) {
      trajectory->time_array_hours[loop_ub + 1] =
          idx_outer * ((double)loop_ub + 1.0);
    }
  } else {
    idx_outer = duration_hours / 99.0;
    for (loop_ub = 0; loop_ub < 98; loop_ub++) {
      trajectory->time_array_hours[loop_ub + 1] =
          ((double)loop_ub + 1.0) * idx_outer;
    }
  }
  /*  ====== 初始化轨迹矩阵 ====== */
  /*  ====== 计算编队演化 ====== */
  for (loop_ub = 0; loop_ub < 100; loop_ub++) {
    b_struct_T a__1[11];
    struct_T satellites_lvlh_t[11];
    d = epoch_jd + trajectory->time_array_hours[loop_ub] / 24.0;
    trajectory->time_array_jd[loop_ub] = d;
    /*  传播到目标时间 */
    propagate_formation(satellites_coe, epoch_jd, d, formation_config->p,
                        formation_config->s, formation_config->l,
                        formation_config->alpha, formation_config->theta,
                        satellites_lvlh_t, a__1);
    /*  记录轨迹 */
    for (sat_idx = 0; sat_idx < 11; sat_idx++) {
      i = sat_idx + 11 * loop_ub;
      trajectory->x[i] = satellites_lvlh_t[sat_idx].x;
      trajectory->y[i] = satellites_lvlh_t[sat_idx].y;
      trajectory->z[i] = satellites_lvlh_t[sat_idx].z;
    }
  }
  /*  ====== 打包轨迹输出 ====== */
}

/*
 * File trailer for formation_design_total.c
 *
 * [EOF]
 */
