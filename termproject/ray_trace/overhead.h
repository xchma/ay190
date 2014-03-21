#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define YES 1
#define NO 0

/* Define ray structure */
struct ray_structure
{
    int ray_num;                 // ray number
    double x;                    // Cartesian coordinate x
    double y;                    // Cartesian coordinate y
    double r;                    // polar coordinate r
    double theta;                // polar coordinate theta
    double nx;                   // photon propagate direction n_x
    double ny;                   // photon propagate direction n_y
    int r_cell;                  // which r cell is the photon in
    int theta_cell;              // which theta cell is the photon in
    int out_flag;                // a flag indicates if the photon escapes
};

/* claim the functions in use */
void allocate_ray(void);
double rand_float(void);
void initialize_ray(void);
double find_step_size(int);
void move_one_step(int);
void free_ray(void);

/*----------------------------------------------------------------------------*/
/* general purpose macros and definitions (never modified) */
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define SQR(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define STR(x) #x
#define AND &&
#define OR ||
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#ifndef PI
#define PI       3.14159265358979323846
#endif
#define ONE_3RD  0.3333333333333333
#define TWO_3RDS 0.6666666666666667
#define SMALL_NUMBER 1.0e-10
#define LARGE_NUMBER 1.0e+10
#define KAPPA 2.0               // opacity
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
/*----------------------------------------------------------------------------*/
/* DMAX,DMIN,IMAX,IMIN defs */
#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)