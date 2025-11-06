/*
 * File: main.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 11-Aug-2025 16:57:37
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include Files */
#include "main.h"
#include "formation_design_total.h"
#include "formation_design_total_terminate.h"
#include "formation_design_total_types.h"
#include "rt_nonfinite.h"
#include <stdio.h>    /* printf */
#include <math.h>     /* sqrt */
#include <stdlib.h>   /* malloc, free */

/* 函数原型声明 */
static double argInit_real_T(void);
static struct0_T argInit_struct0_T(void);
static void argInit_struct1_T(struct1_T* result);
static void save_results(struct2_T satellites_coe[],
    struct3_T satellites_lvlh[],
    struct4_T satellites_j2000[],
    struct5_T* formation_params,
    struct6_T* trajectory,
    int num_satellites);

/* 初始化 double */
static double argInit_real_T(void)
{
    return 0.0;
}

/* 初始化 struct0_T */
static struct0_T argInit_struct0_T(void)
{
    struct0_T result;
    double result_tmp = argInit_real_T();
    result.e = result_tmp;
    result.i = result_tmp;
    result.Omega = result_tmp;
    result.omega = result_tmp;
    result.M = result_tmp;
    result.a = result_tmp;
    return result;
}

/* 初始化 struct1_T */
static void argInit_struct1_T(struct1_T* result)
{
    double result_tmp = argInit_real_T();
    result->s = result_tmp;
    result->l = result_tmp;
    result->alpha = result_tmp;
    result->theta = result_tmp;
    result->pphi = result_tmp;
    result->p_inner = result_tmp;
    result->p_outer = result_tmp;
    result->num_inner = result_tmp;
    result->num_outer = result_tmp;
    result->p = result_tmp;
}

/* main */
int main(int argc, char** argv)
{
    (void)argc;
    (void)argv;
    main_formation_design_total();
    formation_design_total_terminate();
    return 0;
}

/* 主函数 */
void main_formation_design_total(void)
{
    struct0_T chief_coe;
    struct1_T cfg;
    struct2_T satellites_coe[11];
    struct3_T satellites_lvlh[11];
    struct4_T satellites_j2000[11];
    struct5_T formation_params;
    struct6_T trajectory;

    /* 主星轨道要素 */
    chief_coe.a = 6899807;
    chief_coe.e = 0.0011;
    chief_coe.i = 97.5008;
    chief_coe.Omega = 347.0438;
    chief_coe.omega = 90;
    chief_coe.M = 0;

    double epoch_jd = 2460463.0;

    /* 编队构型 */
    cfg.p = 50000;
    cfg.s = cfg.p * sqrt(3.0) / 2.0;
    cfg.l = 0;
    cfg.alpha = 90;
    cfg.theta = 0;
    cfg.pphi = 0;
    cfg.p_inner = 0.5 * cfg.p;
    cfg.p_outer = cfg.p;
    cfg.num_inner = 3;
    cfg.num_outer = 6;

    double num_satellites = 9;
    double formation_type = 3;
    double duration_hours = 6;

    printf("=========== 测试参数 ===========\n");
    printf("主星轨道要素:\n");
    printf("  a = %.2f m\n", chief_coe.a);
    printf("  e = %.4f\n", chief_coe.e);
    printf("  i = %.4f deg\n", chief_coe.i);
    printf("  Omega = %.4f deg\n", chief_coe.Omega);
    printf("  omega = %.2f deg\n", chief_coe.omega);
    printf("  M = %.2f deg\n", chief_coe.M);
    printf("\n");

    printf("编队构型参数:\n");
    printf("  p = %.2f m\n", cfg.p);
    printf("  s = %.2f m\n", cfg.s);
    printf("  l = %.2f m\n", cfg.l);
    printf("  alpha = %.2f deg\n", cfg.alpha);
    printf("  theta = %.2f deg\n", cfg.theta);
    printf("  pphi = %.2f deg\n", cfg.pphi);
    printf("  p_inner = %.2f m\n", cfg.p_inner);
    printf("  p_outer = %.2f m\n", cfg.p_outer);
    printf("  num_inner = %.0f\n", cfg.num_inner);
    printf("  num_outer = %.0f\n", cfg.num_outer);
    printf("\n");

    printf("其他参数:\n");
    printf("  epoch_jd = %.1f\n", epoch_jd);
    printf("  num_satellites = %.0f\n", num_satellites);
    printf("  formation_type = %.0f\n", formation_type);
    printf("  duration_hours = %.1f\n", duration_hours);
    printf("================================\n\n");

    /* 调用生成函数 */
    formation_design_total(&chief_coe, epoch_jd, &cfg,
        num_satellites, formation_type, duration_hours,
        satellites_coe, satellites_lvlh, satellites_j2000,
        &formation_params, &trajectory);

    /* 保存结果 */
    save_results(satellites_coe, satellites_lvlh, satellites_j2000,
        &formation_params, &trajectory, (int)formation_params.num_satellites);
}

/* 保存 CSV */
static void save_results(struct2_T satellites_coe[], struct3_T satellites_lvlh[],
    struct4_T satellites_j2000[], struct5_T* formation_params,
    struct6_T* trajectory, int num_satellites)
{
    FILE* f;
    int i;

    /* 1. COE */
    errno_t err = fopen_s(&f, "satellites_coe.csv", "w");
    if (err != 0 || !f) { perror("satellites_coe.csv"); return; }
    fprintf(f, "sat_id,a,e,i,Omega,omega,M,epoch_jd\n");
    for (i = 0; i < num_satellites; i++) {
        fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%.10f\n",
            satellites_coe[i].sat_id, satellites_coe[i].a, satellites_coe[i].e,
            satellites_coe[i].i, satellites_coe[i].Omega, satellites_coe[i].omega,
            satellites_coe[i].M, satellites_coe[i].epoch_jd);
    }
    fclose(f);

    /* 2. LVLH */
    err = fopen_s(&f, "satellites_lvlh.csv", "w");
    if (err != 0 || !f) { perror("satellites_lvlh.csv"); return; }
    fprintf(f, "sat_id,x,y,z,vx,vy,vz,epoch_jd\n");
    for (i = 0; i < num_satellites; i++) {
        fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%.10f\n",
            satellites_lvlh[i].sat_id, satellites_lvlh[i].x, satellites_lvlh[i].y,
            satellites_lvlh[i].z, satellites_lvlh[i].vx, satellites_lvlh[i].vy,
            satellites_lvlh[i].vz, satellites_lvlh[i].epoch_jd);
    }
    fclose(f);


    /* 4. Formation Params */
     err = fopen_s(&f, "formation_params.csv", "w");
    if (err != 0 || !f) { perror("formation_params.csv"); return; }
    fprintf(f, "chief_a,%g\n", formation_params->chief_coe.a);
    fprintf(f, "chief_e,%g\n", formation_params->chief_coe.e);
    fprintf(f, "chief_i,%g\n", formation_params->chief_coe.i);
    fprintf(f, "chief_Omega,%g\n", formation_params->chief_coe.Omega);
    fprintf(f, "chief_omega,%g\n", formation_params->chief_coe.omega);
    fprintf(f, "chief_M,%g\n", formation_params->chief_coe.M);
    fprintf(f, "p,%g\n", formation_params->formation_config.p);
    fprintf(f, "s,%g\n", formation_params->formation_config.s);
    fprintf(f, "l,%g\n", formation_params->formation_config.l);
    fprintf(f, "alpha,%g\n", formation_params->formation_config.alpha);
    fprintf(f, "theta,%g\n", formation_params->formation_config.theta);
    fprintf(f, "pphi,%g\n", formation_params->formation_config.pphi);
    fprintf(f, "p_inner,%g\n", formation_params->formation_config.p_inner);
    fprintf(f, "p_outer,%g\n", formation_params->formation_config.p_outer);
    fprintf(f, "num_inner,%g\n", formation_params->formation_config.num_inner);
    fprintf(f, "num_outer,%g\n", formation_params->formation_config.num_outer);
    fprintf(f, "num_satellites,%g\n", formation_params->num_satellites);
    fprintf(f, "formation_type,%g\n", formation_params->formation_type);
    fprintf(f, "epoch_jd,%.10f\n", formation_params->epoch_jd);
    fprintf(f, "orbital_period,%g\n", formation_params->orbital_period);
    fprintf(f, "mean_motion,%g\n", formation_params->mean_motion);
    int str_len = formation_params->type_name.size[1];
    if (str_len > 0 && formation_params->type_name.data[str_len - 1] == '\0') {
        str_len--; /* 去掉 '\0' */
    }
    fprintf(f, "type_name,%.*s\n", str_len, formation_params->type_name.data);
    fclose(f);


    fclose(f);
}
