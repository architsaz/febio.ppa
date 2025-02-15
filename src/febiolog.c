#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "logfile_types.h"
#include "input_type.h"
#include "common.h"
#include "mesh.h"
// read stress tensor in the log file of febio
int readfebiolog(char *path,int nelem2, double **st2, read_time logtime)
{

    int e = 0;
    double *st;
    // Open log file :
    FILE *fptr = fopen(path, "r");
    char str[256];
    if (fptr == NULL)
    {
        fprintf(stderr, "* ERROR: there is no file in this path : %s.\n", path);
        exit(EXIT_FAILURE);
    }
    // find the number of shell element
    int nelem = 0;
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "	Number of shell elements ....................... : %d", &nelem) == 1)
            break;
    }
    fclose(fptr);
    printf("* number of shell element is %d\n", nelem);
    if (nelem2 != nelem)
    {
        fprintf(stderr, "ERROR: Number of shell element in the file : %s does not match with what zfem file on data directory", path);
        exit(EXIT_FAILURE);
    }
    // find the maximum time in the file
    double time_value, max_time;
    max_time = 0.0;
    fptr = fopen(path, "r");
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "Time = %lf", &time_value) == 1)
            max_time = MAX(time_value, max_time);
    }
    fclose(fptr);
    printf("* %s study execute till time = %lf\n", path, max_time);
    // save the stress tensor
    st = calloc(9 * (size_t)nelem, sizeof(*st));
    double logtime_value = 0.0;
    if (logtime == end_first_step)
        logtime_value = 1.0;
    if (logtime == end_second_step)
        logtime_value = 2.0;
    if (logtime == time_max)
        logtime_value = max_time;
    // st [9] = [sxx,sxy,sxz;syx,syy,syz;szx,szy,szz]
    fptr = fopen(path, "r");
    int junk, nscan = 0;
    int find_time = 0;
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "Time = %lf", &time_value) == 1)
        {
            if (time_value == logtime_value)
            {
                printf("* start to read stress tensor at time: %lf\n",time_value);
                find_time++;
                fgets(str, 256, fptr);
                for (int ele = 0; ele < nelem; ele++)
                {
                    fgets(str, 256, fptr);
                    nscan = 0;
                    // log_st[6]= [sxx(0),syy(4),szz(8),sxy(1)(3),syz(5)(7),sxz(2)(6)]
                    nscan = sscanf(str, "%d %lf %lf %lf %lf %lf %lf", &junk, &st[9 * ele], &st[9 * ele + 4], &st[9 * ele + 8],
                                   &st[9 * ele + 1], &st[9 * ele + 5], &st[9 * ele + 2]);
                    if (nscan != 7)
                    {
                        fprintf(stderr, "there is error on number of element in line %d", ele);
                        exit(EXIT_FAILURE);
                    }
                    st[9 * ele + 3] = st[9 * ele + 1];
                    st[9 * ele + 7] = st[9 * ele + 5];
                    st[9 * ele + 6] = st[9 * ele + 2];
                }
            }
        }
    }
    if (find_time == 0){
        fprintf(stderr, "ERROR: can not find time %lf in the log file!\n", max_time);
        exit(EXIT_FAILURE);
    }
    // for(int ele=0;ele<10;ele++){
    //     printf("ele: %d ",ele);
    //     for(int i =0;i<9;i++) printf("%lf ",st[9*ele+i]);
    //     printf("\n");
    // }
    fclose(fptr);
    printf("* the stress tensor saved !\n");

    *st2 = st;
    return e;
}
int read_wallmask(char *path,int nelem, input *inp, int **Melem2)
{
    int e = 0;
    int *Melem = calloc((size_t)nelem, sizeof(*Melem));
    if (Melem == NULL)
    {
        fprintf(stderr, "Memory allocation failed for Melem.\n");
        return -1;
    }

    /* Open the file */
    FILE *fptr = fopen(path, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
        free(Melem); // Free allocated memory before returning
        return -1;
    }

    /* Read all lines of the file */
    int buffer = 100;
    char *str, line[buffer];
    int nscan, iline;

    /* Read labels of all elements */
    int temp = 0;
    for (iline = 0; iline < nelem; iline++)
    {
        str = edit_endline_character(line, buffer, fptr);
        nscan = sscanf(str, "%d", &temp);
        if (nscan != 1)
        {
            fprintf(stderr, "ERROR: Incorrect number of entries on POINTS line.\n");
            free(Melem);  // Free before returning
            fclose(fptr); // Close the file before returning
            return -1;
        }

        // Check the value of the Melem
        for (int k = 0; k < inp->label_num; k++)
        {
            if (temp == inp->label[k])
            {
                Melem[iline] = inp->label[k];
            }
        }
    }

    /* Close the file */
    if (fclose(fptr) == EOF)
    {
        printf("Error closing %s\n", path);
        free(Melem); // Free before returning in case of error
        return -1;
    }

    /* Assign the allocated array to the output pointer */
    *Melem2 = Melem;
    printf("* Done reading .WALL mask file!\n");

    return e;
}
int read_regionmask(char *path, int npoin,int nelem,int *elems, input *inp, int **region_id2, int **region_idp2)
{
    int e = 0;

    // Allocate memory for region_id
    int *region_id = malloc((size_t)npoin * sizeof(*region_id));
    if (region_id == NULL)
    {
        fprintf(stderr, "Memory allocation failed for region_id.\n");
        return -1;
    }

    // Open the file
    FILE *fptr = fopen(path, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
        free(region_id); // Free allocated memory before returning
        return -1;
    }

    /* Read all lines of the file */
    int buffer = 100;
    char *str;
    char line[buffer];
    int endcount = 0;
    int nscan, iline;
    char test[20];

    // Reading region labels
    while (1)
    {
        str = edit_endline_character(line, buffer, fptr);
        nscan = sscanf(str, "%s", test);

        if (!strcmp(test, "regions"))
        {
            str = edit_endline_character(line, buffer, fptr);
            nscan = sscanf(str, "%s", test);

            if (!strcmp(test, "1"))
            {
                for (iline = 0; iline < npoin; iline++)
                {
                    str = edit_endline_character(line, buffer, fptr);
                    nscan = sscanf(str, "%d", &region_id[iline]);

                    // Error check
                    if (nscan != 1)
                    {
                        fprintf(stderr, "ERROR: Incorrect data on line %d.\n", iline + 1);
                        free(region_id);
                        fclose(fptr);
                        return -1;
                    }

                    // Check if value is in colorid array
                    int checkflag = 0;
                    for (int k = 0; k < inp->colorid_num; k++)
                    {
                        if (region_id[iline] == inp->colorid[k])
                        {
                            checkflag++;
                            break;
                        }
                    }
                    if (checkflag != 1)
                    {
                        printf("ERROR: Value of element %d (value=%d) is not in colorId array\n", iline, region_id[iline]);
                        free(region_id);
                        fclose(fptr);
                        return -1;
                    }
                }
                endcount += 1;
            }
        }

        if (endcount == 1)
        {
            printf("* Done Reading region mask file!\n");
            break;
        }
    }

    // Close the file
    if (fclose(fptr) == EOF)
    {
        printf("Error closing %s\n", path);
        free(region_id);
        return -1;
    }

    // Allocate memory for region_id_ele
    int *region_id_ele = malloc((size_t)nelem * sizeof(*region_id_ele));
    if (region_id_ele == NULL)
    {
        fprintf(stderr, "Memory allocation failed for region_id_ele.\n");
        free(region_id);
        return -1;
    }

    // Fill region_id_ele based on region_id
    int points[3] = {0, 0, 0};
    for (int ele = 0; ele < nelem; ele++)
    {
        points[0] = elems[3 * ele];
        points[1] = elems[3 * ele + 1];
        points[2] = elems[3 * ele + 2];

        // Assign region ID for the element based on one point (e.g., points[0] - 1)
        region_id_ele[ele] = region_id[points[0] - 1];
    }

    /* Return results */
    *region_id2 = region_id_ele;
    *region_idp2 = region_id;

    return e;
}
// find unidirectional or bidirectional stress region mask:
void unibimask(int nelem, double *smax, double *smin, double critical_ratio, double threshold, int **sdir2)
{
    int *sdir = (int *)calloc((size_t)nelem, sizeof(int));
    for (int ele = 0; ele < nelem; ele++)
    {
        if (fabs(smin[ele]) < 1)
            continue; // avoid unpressurized region
        double ratio = fabs(smin[ele]) / fabs(smax[ele]);
        sdir[ele] = (ratio < critical_ratio) ? 1 : 2; // unidiretional : 1 bidirectional : 2
        if (fabs(smax[ele]) <= threshold)
        {
            // sdir[ele] = (sdir[ele] == 4) ? 2 : 1;
            sdir[ele] = 0;
        }
    }
    *sdir2 = sdir;
}