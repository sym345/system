/*
 * File: formation_design_total.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

#ifndef FORMATION_DESIGN_TOTAL_H
#define FORMATION_DESIGN_TOTAL_H

/* Include Files */
#include "formation_design_total_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern __declspec(dllexport) void formation_design_total(
    const struct0_T *chief_coe, double epoch_jd,
    const struct1_T *formation_config, double num_satellites,
    double formation_type, double duration_hours, struct2_T satellites_coe[11],
    struct3_T satellites_lvlh[11], struct4_T satellites_j2000[11],
    struct5_T *formation_params, struct6_T *trajectory);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for formation_design_total.h
 *
 * [EOF]
 */
