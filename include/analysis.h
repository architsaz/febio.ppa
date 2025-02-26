#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "mesh_types.h"
int analz_double(mesh *M, double *area, int *Melem, int *relems, int *bleb, double *Eval_max, double *field, char *casename, char *study, char *filename, char *field_name);
int analz_int(mesh *M, double *area, int *Melem, int *relems, int *bleb, double *Eval_max, int *field, char *casename, char *study, char *filename, char *field_name);
#endif