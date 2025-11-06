/*
 * File: formation_design_total_types.h
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

#ifndef FORMATION_DESIGN_TOTAL_TYPES_H
#define FORMATION_DESIGN_TOTAL_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  double a;
  double e;
  double i;
  double Omega;
  double omega;
  double M;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  double p;
  double s;
  double l;
  double alpha;
  double theta;
  double pphi;
  double p_inner;
  double p_outer;
  double num_inner;
  double num_outer;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  double sat_id;
  double a;
  double e;
  double i;
  double Omega;
  double omega;
  double M;
  double epoch_jd;
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  double sat_id;
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double epoch_jd;
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_struct6_T
#define typedef_struct6_T
typedef struct {
  double x[1100];
  double y[1100];
  double z[1100];
  double time_array_hours[100];
  double time_array_jd[100];
} struct6_T;
#endif /* typedef_struct6_T */

#ifndef struct_emxArray_real_T_3x3
#define struct_emxArray_real_T_3x3
struct emxArray_real_T_3x3 {
  double data[9];
  int size[2];
};
#endif /* struct_emxArray_real_T_3x3 */
#ifndef typedef_emxArray_real_T_3x3
#define typedef_emxArray_real_T_3x3
typedef struct emxArray_real_T_3x3 emxArray_real_T_3x3;
#endif /* typedef_emxArray_real_T_3x3 */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  double sat_id;
  emxArray_real_T_3x3 r;
  emxArray_real_T_3x3 v;
  double epoch_jd;
} struct4_T;
#endif /* typedef_struct4_T */

#ifndef struct_emxArray_char_T_1x24
#define struct_emxArray_char_T_1x24
struct emxArray_char_T_1x24 {
  char data[24];
  int size[2];
};
#endif /* struct_emxArray_char_T_1x24 */
#ifndef typedef_emxArray_char_T_1x24
#define typedef_emxArray_char_T_1x24
typedef struct emxArray_char_T_1x24 emxArray_char_T_1x24;
#endif /* typedef_emxArray_char_T_1x24 */

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  struct0_T chief_coe;
  struct1_T formation_config;
  double num_satellites;
  double formation_type;
  double epoch_jd;
  double orbital_period;
  double mean_motion;
  emxArray_char_T_1x24 type_name;
} struct5_T;
#endif /* typedef_struct5_T */

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#endif
/*
 * File trailer for formation_design_total_types.h
 *
 * [EOF]
 */
