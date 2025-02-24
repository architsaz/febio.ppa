#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mesh_types.h"
#include "analysis_type.h"
// Function to sum the array 
double sumarr(double *arr, int size)
{
	double sum = 0.0;
	if (size == 0)
		return sum;
	for (int ele = 0; ele < size; ele++)
	{
		sum += arr[ele];
	}
	return sum;
}
// Function to sort the array (used for calculating the median)
void sort_array(double arr[], int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        for (int j = i + 1; j < size; j++)
        {
            if (arr[i] > arr[j])
            {
                double temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
    }
}
// Function to find the maximum size among arrays
int findMaxSize(int numArrays, int sizes[])
{
    int max = sizes[0];
    for (int i = 1; i < numArrays; i++)
    {
        if (sizes[i] > max)
        {
            max = sizes[i];
        }
    }
    return max;
}
// Function to write multiple arrays with headers (names) to a text file
void saveMultipleArraysToFile(const char *path, int numArrays, void *arrays[], int sizes[], DataType types[], const char *headers[])
{
    FILE *file = fopen(path, "w");

    if (file == NULL)
    {
        printf("Error opening file!\n");
        return;
    }

    // Print headers (array names)
    for (int j = 0; j < numArrays; j++)
    {
        fprintf(file, "%s\t", headers[j]); // Print each header name
    }
    fprintf(file, "\n"); // Newline after headers

    int maxRows = findMaxSize(numArrays, sizes); // Find the maximum length of the arrays

    // Loop through each row of the data
    for (int i = 0; i < maxRows; i++)
    {
        // Loop through each array (column)
        for (int j = 0; j < numArrays; j++)
        {
            if (i < sizes[j])
            {
                // Print based on the data type of the array
                if (types[j] == INT_TYPE)
                {
                    fprintf(file, "%d\t", ((int *)arrays[j])[i]); // Print int
                }
                else if (types[j] == FLOAT_TYPE)
                {
                    fprintf(file, "%8.2f\t", ((float *)arrays[j])[i]); // Print float
                }
                else if (types[j] == DOUBLE_TYPE)
                {
                    fprintf(file, "%9.2lf\t", ((double *)arrays[j])[i]); // Print double
                }
                else if (types[j] == CHAR_TYPE)
                {
                    fprintf(file, "%c\t", ((char *)arrays[j])[i]); // Print single char
                }
                else if (types[j] == STRING_TYPE)
                {
                    fprintf(file, "%s\t", ((char **)arrays[j])[i]); // Print string (array of char*)
                }
            }
            else
            {
                fprintf(file, "\t"); // Print an empty tab for missing elements
            }
        }
        fprintf(file, "\n"); // Newline after each row
    }

    fclose(file);
    printf("* Data written successfully to %s!\n", path);
}
// Function to calculate the mean
double calculate_mean(double *arr, int size, double *weight)
{
    double sum = 0.0;
    double sum_weight = 0.0;
    for (int i = 0; i < size; i++)
    {
        sum += arr[i] * weight[i];
        sum_weight += weight[i];
    }
    return sum / sum_weight;
}
// Function to calculate the median
double calculate_median(double arr[], int size)
{
    sort_array(arr, size);

    if (size % 2 == 0)
    {
        return (arr[size / 2 - 1] + arr[size / 2]) / 2.0;
    }
    else
    {
        return arr[size / 2];
    }
}
// Function to find the maximum value
double find_max(double arr[], int size)
{
    double max = arr[0];
    for (int i = 1; i < size; i++)
    {
        if (arr[i] > max)
        {
            max = arr[i];
        }
    }
    return max;
}
// Function to find the minimum value
double find_min(double arr[], int size)
{
    double min = arr[0];
    for (int i = 1; i < size; i++)
    {
        if (arr[i] < min)
        {
            min = arr[i];
        }
    }
    return min;
}
// Function to calculate the standard deviation
double calculate_stddev(double arr[], int size, double mean, double *weight)
{
    double sum = 0.0;
    double sum_weight = 0.0;
    for (int i = 0; i < size; i++)
    {
        sum += weight[i] * pow(arr[i] - mean, 2);
        sum_weight += weight[i];
    }
    return sqrt(sum / sum_weight);
}
// Calling defined statictic analysis function in ordered
void mystat(double *arr, int n, double *area, double **output1)
{
    static double *output;
    output = calloc((size_t)5, sizeof(*output));
    output[0] = calculate_mean(arr, n, area);
    output[1] = find_max(arr, n);
    output[2] = find_min(arr, n);
    output[3] = calculate_stddev(arr, n, output[0], area);
    *output1 = output;
}
// analysis the max stress @ color region mask and region mask
int analzs(mesh *M, double *area, int *Melem, int *relems, double *Eval_max, double *field, char *casename, char *study, char *filename)
{
    int e = 0;

    int nele_red, nele_yel, nele_wht, nele_rupt, nele_press, nele_dom, nele_bod, nele_nek, nele_part, nele_aneu;
    nele_red = nele_yel = nele_wht = nele_rupt = nele_press = nele_dom = nele_bod = nele_nek = nele_part = nele_aneu = 0;
    double *field_red, *field_yel, *field_wht, *field_rupt, *field_press, *field_aneu, *field_dom, *field_bod, *field_nek, *field_part;
    double *area_red, *area_yel, *area_wht, *area_rupt, *area_press, *area_aneu, *area_dom, *area_bod, *area_nek, *area_part;
    // find the size of each domain
    for (int ele = 0; ele < M->nelem; ele++)
    {
        if (Eval_max[ele] < 1)
            continue;
        nele_press++;
        if (Melem[ele] == 1)
            nele_red++;
        if (Melem[ele] == 4)
            nele_yel++;
        if (Melem[ele] == 7)
            nele_wht++;
        if (Melem[ele] == 9)
            nele_rupt++;    
        if (relems[ele] == 16)
            nele_dom++;
        if (relems[ele] == 8)
            nele_bod++;
        if (relems[ele] == 4)
            nele_nek++;
        if (relems[ele] == 1)
            nele_part++;
        if (relems[ele] == 4 || relems[ele] == 8 || relems[ele] == 16)
            nele_aneu++;
    }
    printf("nele_press: %d \n", nele_press);
    printf("nele_red: %d \n", nele_red);
    printf("nele_yel: %d \n", nele_yel);
    printf("nele_wht: %d \n", nele_wht);
    printf("nele_rupt: %d \n", nele_rupt);
    printf("nele_dom: %d \n", nele_dom);
    printf("nele_bod: %d \n", nele_bod);
    printf("nele_nek: %d \n", nele_nek);
    printf("nele_part: %d \n", nele_part);
    printf("nele_aneu: %d \n", nele_aneu);
    field_red = calloc((size_t)nele_red, sizeof(double));
    field_yel = calloc((size_t)nele_yel, sizeof(double));
    field_wht = calloc((size_t)nele_wht, sizeof(double));
    field_rupt = calloc((size_t)nele_rupt, sizeof(double));
    field_aneu = calloc((size_t)nele_aneu, sizeof(double));
    field_press = calloc((size_t)nele_press, sizeof(double));
    field_dom = calloc((size_t)nele_dom, sizeof(double));
    field_bod = calloc((size_t)nele_bod, sizeof(double));
    field_nek = calloc((size_t)nele_nek, sizeof(double));
    field_part = calloc((size_t)nele_part, sizeof(double));
    area_red = calloc((size_t)nele_red, sizeof(double));
    area_yel = calloc((size_t)nele_yel, sizeof(double));
    area_wht = calloc((size_t)nele_wht, sizeof(double));
    area_rupt = calloc((size_t)nele_rupt, sizeof(double));
    area_aneu = calloc((size_t)nele_aneu, sizeof(double));
    area_press = calloc((size_t)nele_press, sizeof(double));
    area_dom = calloc((size_t)nele_dom, sizeof(double));
    area_bod = calloc((size_t)nele_bod, sizeof(double));
    area_nek = calloc((size_t)nele_nek, sizeof(double));
    area_part = calloc((size_t)nele_part, sizeof(double));
    nele_red = nele_yel = nele_wht = nele_rupt = nele_press = nele_dom = nele_bod = nele_nek = nele_part = nele_aneu = 0;
    for (int ele = 0; ele < M->nelem; ele++)
    {
        if (Eval_max[ele] < 1)
            continue;
        field_press[nele_press] = fabs(field[ele]);
        area_press[nele_press] = area[ele];
        nele_press++;
        if (Melem[ele] == 1)
        {
            field_red[nele_red] = fabs(field[ele]);
            area_red[nele_red] = area[ele];
            nele_red++;
        }

        if (Melem[ele] == 4)
        {
            field_yel[nele_yel] = fabs(field[ele]);
            area_yel[nele_yel] = area[ele];
            nele_yel++;
        }

        if (Melem[ele] == 7)
        {
            field_wht[nele_wht] = fabs(field[ele]);
            area_wht[nele_wht] = area[ele];
            nele_wht++;
        }

        if (Melem[ele] == 9)
        {
            field_rupt[nele_rupt] = fabs(field[ele]);
            area_rupt[nele_rupt] = area[ele];
            nele_rupt++;
        }

        if (relems[ele] == 16)
        {
            field_dom[nele_dom] = fabs(field[ele]);
            area_dom[nele_dom] = area[ele];
            nele_dom++;
        }

        if (relems[ele] == 8)
        {
            field_bod[nele_bod] = fabs(field[ele]);
            area_bod[nele_bod] = area[ele];
            nele_bod++;
        }

        if (relems[ele] == 4)
        {
            field_nek[nele_nek] = fabs(field[ele]);
            area_nek[nele_nek] = area[ele];
            nele_nek++;
        }

        if (relems[ele] == 1)
        {
            field_part[nele_part] = fabs(field[ele]);
            area_part[nele_part] = area[ele];
            nele_part++;
        }

        if (relems[ele] == 4 || relems[ele] == 8 || relems[ele] == 16)
        {
            field_aneu[nele_aneu] = fabs(field[ele]);
            area_aneu[nele_aneu] = area[ele];
            nele_aneu++;
        }
    }
    // Calculate statistics
    double *stat_red, *stat_yel, *stat_wht, *stat_rupt, *stat_aneu, *stat_dom, *stat_bod, *stat_nek, *stat_part, *stat_press;
    double stat_empty[4] = {0, 0, 0, 0};
    if (nele_red != 0)
    {
        mystat(field_red, nele_red, area_red, &stat_red);
    }
    else
    {
        stat_red = stat_empty;
    }
    if (nele_yel != 0)
    {
        mystat(field_yel, nele_yel, area_yel, &stat_yel);
    }
    else
    {
        stat_yel = stat_empty;
    }
    if (nele_wht != 0)
    {
        mystat(field_wht, nele_wht, area_wht, &stat_wht);
    }
    else
    {
        stat_wht = stat_empty;
    }
    if (nele_rupt != 0)
    {
        mystat(field_rupt, nele_rupt, area_rupt, &stat_rupt);
    }
    else
    {
        stat_rupt = stat_empty;
    }
    mystat(field_aneu, nele_aneu, area_aneu, &stat_aneu);
    mystat(field_dom, nele_dom, area_dom, &stat_dom);
    mystat(field_bod, nele_bod, area_bod, &stat_bod);
    mystat(field_nek, nele_nek, area_nek, &stat_nek);
    mystat(field_part, nele_part, area_part, &stat_part);
    mystat(field_press, nele_press, area_press, &stat_press);
    printf("red: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_red, nele_red), stat_red[0], stat_red[1], stat_red[2], stat_red[3]);
    printf("yel: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_yel, nele_yel), stat_yel[0], stat_yel[1], stat_yel[2], stat_yel[3]);
    printf("wht: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_wht, nele_wht), stat_wht[0], stat_wht[1], stat_wht[2], stat_wht[3]);
    printf("rupt: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_rupt, nele_rupt), stat_rupt[0], stat_rupt[1], stat_rupt[2], stat_rupt[3]);
    printf("dom: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_dom, nele_dom), stat_dom[0], stat_dom[1], stat_dom[2], stat_dom[3]);
    printf("bod: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_bod, nele_bod), stat_bod[0], stat_bod[1], stat_bod[2], stat_bod[3]);
    printf("nek: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_nek, nele_nek), stat_nek[0], stat_nek[1], stat_nek[2], stat_nek[3]);
    printf("ane: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_aneu, nele_aneu), stat_aneu[0], stat_aneu[1], stat_aneu[2], stat_aneu[3]);
    printf("par: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_part, nele_part), stat_part[0], stat_part[1], stat_part[2], stat_part[3]);
    printf("pre: area: %0.5lf mean: %9.2lf max: %9.2lf min: %9.2lf stddev: %9.2lf\n", sumarr(area_press, nele_press), stat_press[0], stat_press[1], stat_press[2], stat_press[3]);

    // save in a table txt
    // arrays of different types and sizes
    char *arr1[4];
    for (int i = 0; i < 4; i++)
    {
        arr1[i] = (char *)malloc(strlen(casename) + 1); // +1 for the null terminator
        if (arr1[i] != NULL)
        {
            strcpy(arr1[i], casename);
        }
        else
        {
            // Handle memory allocation failure
            fprintf(stderr, "ERROR: Memory allocation failed for arr1[%d]\n", i);
            exit(EXIT_FAILURE);
        }
    }
    char *arr2[4];
    for (int i = 0; i < 4; i++)
    {
        arr2[i] = (char *)malloc(strlen(study) + 1); // +1 for the null terminator
        if (arr2[i] != NULL)
        {
            strcpy(arr2[i], study);
        }
        else
        {
            // Handle memory allocation failure
            fprintf(stderr, "ERROR: Memory allocation failed for arr2[%d]\n", i);
            exit(EXIT_FAILURE);
        }
    }
    char *arr3[4] = {"mean", "max", "min", "stddev"};

    // Number of arrays and their sizes
    int numArrays = 13;
    int sizes[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}; // Sizes of the arrays

    // Array of pointers to the arrays
    void *arrays[] = {(void *)arr1, (void *)arr2, (void *)arr3, (void *)stat_aneu, (void *)stat_red, (void *)stat_yel, (void *)stat_wht, (void *)stat_rupt, (void *)stat_dom,
                      (void *)stat_bod, (void *)stat_nek, (void *)stat_part, (void *)stat_press};

    // Data types of each array (int, float, char, string)
    DataType types[] = {STRING_TYPE, STRING_TYPE, STRING_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE};

    // Array of column headers (names)
    const char *headers[] = {"Casename", "Study", "stat_para", "stat_aneu", "stat_red", "stat_yel", "stat_wht", "stat_rupt", "stat_dom", "stat_bod", "stat_nek", "stat_part", "stat_press"};

    // Save the arrays to the file with headers
    saveMultipleArraysToFile(filename, numArrays, arrays, sizes, types, headers);

    // free memory
    for (int i = 0; i < 4; i++)
    {
        free(arr1[i]); // Free memory allocated for arr1
        free(arr2[i]); // Free memory allocated for arr2
    }

    // Free the dynamically allocated smax and area arrays
    free(field_red);
    free(field_yel);
    free(field_wht);
    free(field_rupt);
    free(field_aneu);
    free(field_press);
    free(field_dom);
    free(field_bod);
    free(field_nek);
    free(field_part);

    free(area_red);
    free(area_yel);
    free(area_wht);
    free(area_rupt);
    free(area_aneu);
    free(area_press);
    free(area_dom);
    free(area_bod);
    free(area_nek);
    free(area_part);

    // Free dynamically allocated stat arrays if they are not pointing to stat_empty
    if (stat_red != stat_empty)
        free(stat_red);
    if (stat_yel != stat_empty)
        free(stat_yel);
    if (stat_wht != stat_empty)
        free(stat_wht);
    if (stat_rupt != stat_empty)
        free(stat_rupt);    
    // Other stat arrays that were dynamically allocated
    free(stat_aneu);
    free(stat_dom);
    free(stat_bod);
    free(stat_nek);
    free(stat_part);
    free(stat_press);

    return e;
}
