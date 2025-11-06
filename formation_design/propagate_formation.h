/*
 * File: propagate_formation.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

#ifndef PROPAGATE_FORMATION_H
#define PROPAGATE_FORMATION_H

/* Include Files */
#include "formation_design_total_internal_types.h"
#include "formation_design_total_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void propagate_formation(const struct2_T satellites_coe[11], double epoch_jd,
                         double target_time_jd, double formation_config_p,
                         double formation_config_s, double formation_config_l,
                         double formation_config_alpha,
                         double formation_config_theta,
                         struct_T satellites_lvlh_t[11],
                         b_struct_T satellites_j2000_t[11]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for propagate_formation.h
 *
 * [EOF]
 */
