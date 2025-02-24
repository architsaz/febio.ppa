#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "mesh_types.h"
int analzs(mesh *M, double *area, int *Melem, int *relems, double *Eval_max, double *field, char *casename, char *study, char *filename);
#endif