#ifndef INPUT_TYPE_H
#define INPUT_TYPE_H
    typedef struct
    {
        // Solver
        char nonlinear_FE[50];
        int symetric_stiff;
        // curvature mask
        double norm_ang, bad_ang;
        double young_highcurv;
        // Neo-Hooken model
        char Mmodel[50]; // isotropic elastic  or  neo-Hookean or coupled Mooney-Rivlin
        double *young_r, *young_l, NJyoung, incyoung;
        double pois, ro;
        // mask status
        int used_cmask, used_rmask;
        // label : <red, yellow, white, cyan, rupture, remain>
        int *label, label_num;
        // boundary condition:
        int *colorid, colorid_num, *fix_region, fix_region_num, *load_region, load_region_num; // region { remain(another aneu),diastal,parent,neck,body,dome}
        int used_BCmask;
        // thickness
        double *thick_r, *thick_l; // {red, yellow, white} [cm]
        int thick_r_num, thick_l_num;
        // pre
        double pres, ultipres; //[dyne/cm^2]             120 mmHg
                            // filename
        char filename[50];
        // logfile
        int print_st;
    } input;
#endif