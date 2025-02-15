#ifndef FEBIOLOG_H
#define FEBIOLOG_H
    #include "input_type.h"
    #include "logfile_types.h"
    int read_regionmask(char *path, int npoin,int nelem,int *elems, input *inp, int **region_id2, int **region_idp2);
    int read_wallmask(char *path,int nelem, input *inp, int **Melem2);
    int readfebiolog(char *path,int nelem2, double **st2, read_time logtime);
    void unibimask(int nelem, double *smax, double *smin, double critical_ratio, double threshold, int **sdir2);
#endif 