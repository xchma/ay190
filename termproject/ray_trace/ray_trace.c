#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "overhead.h"

/* global variables */
struct ray_structure **rays;

int N_r, N_theta, N_ray;
double r_min, r_max, theta_min, theta_max;
double dr, dtheta;
float *density;


int ray(int N_r_cells, int N_theta_cells,
        float radius_min, float radius_max,
        float angle_min, float angle_max,
        int Num_rays, float *density_map,
        float *r_final, float *theta_final){
    
    /* assign global variables */
    N_r = N_r_cells;
    N_theta = N_theta_cells;
    N_ray = Num_rays;
    r_min = radius_min;
    r_max = radius_max;
    theta_min = angle_min;
    theta_max = angle_max;
    dr = (r_max - r_min) / N_r;
    dtheta = (theta_max - theta_min) / N_theta;
    density = density_map;
    
    srand(time(NULL));
    
    int i,j;
    
    /* allocate memory, initialize the ray */
    allocate_ray();
    initialize_ray();
    
    /* start the mail loop */
    for (i=0; i<N_ray; i++) {
        
        if ( !(i % 10000) ) {
            printf("...ray=%d...\n",i);
        }
        
        while (rays[i]->out_flag == NO) {
            move_one_step(i);
        }
        
    }
    
    /* write the final results */
    for (i=0; i<N_ray; i++) {
        
        *(r_final+i) = rays[i]->r;
        *(theta_final+i) = rays[i]->theta;
        
    }
    
    /* free the memory */
    free_ray();
    
    return 1;
    
}

void allocate_ray(void){
    
    int i;
    
    if(N_ray>0){
        if (!(rays = malloc(N_ray*sizeof(struct ray_structure*)))) {
            printf("Fail to allocate memory. (A)\n");
            exit(0);
        }
    }
    
    if (!(rays[0] = malloc(N_ray*sizeof(struct ray_structure)))) {
        printf("Fail to allocate memory. (B)\n");
        exit(0);
    }
    
    for (i=1; i<N_ray; i++) {
        rays[i] = rays[i-1] + 1;
    }
    
}

double rand_float(void){
    
    return (double)rand() / ((double)(RAND_MAX) * (1+SMALL_NUMBER));
    
}

void initialize_ray(void){
    
    int i;
    double theta, theta_min_alter;
    
    for (i=0; i<N_ray; i++) {
        
        theta_min_alter = theta_min - dtheta*SMALL_NUMBER; // to avoid boundary
        theta = theta_min_alter + (theta_max-theta_min_alter) * rand_float();
        
        rays[i]->ray_num = i;
        rays[i]->r = r_min + dr*SMALL_NUMBER; // to avoid cell boundary
        rays[i]->theta = theta;
        rays[i]->r_cell = 0;
        rays[i]->theta_cell = (int)((theta-theta_min) / dtheta);
        rays[i]->nx = cos(theta);
        rays[i]->ny = sin(theta);
        rays[i]->x = r_min*cos(theta);
        rays[i]->y = r_min*sin(theta);
        rays[i]->out_flag = NO;
        
    }
    
}

double find_step_size(int i){
    
    int r_cell, theta_cell;
    double x, y, nx, ny;
    
    x = rays[i]->x;
    y = rays[i]->y;
    nx = rays[i]->nx;
    ny = rays[i]->ny;
    r_cell = rays[i]->r_cell;
    theta_cell = rays[i]->theta_cell;
    
    double s, s_temp;
    double theta, denom, numer;
    double radius, coeff_B, coeff_C, delta;
    
    s = LARGE_NUMBER;
    
    /* lower theta boundary */
    theta = theta_min + theta_cell*dtheta;
    denom = nx*sin(theta) - ny*cos(theta);
    if (denom != 0) {
        
        numer = y*cos(theta) - x*sin(theta);
        s_temp = numer/denom;
        
        if (s_temp > 0) {
            s = MIN(s, s_temp);
        }
        
    }
    
    /* upper theta boundary */
    theta = theta_min + (theta_cell+1)*dtheta;
    denom = nx*sin(theta) - ny*cos(theta);
    if (denom != 0) {
        
        numer = y*cos(theta) - x*sin(theta);
        s_temp = numer/denom;
        
        if (s_temp > 0) {
            s = MIN(s, s_temp);
        }
        
    }
    
    /* inner radius boundary */
    radius = r_min + r_cell*dr;
    coeff_B = 2 * (x*nx + y*ny);
    coeff_C = x*x + y*y - radius*radius;
    delta = coeff_B*coeff_B - 4*coeff_C;
    if (delta > 0 && coeff_B < 0) {
        
        s_temp = (-coeff_B - sqrt(delta)) / 2.0;
        s = MIN(s, s_temp);
        
    }
    
    /* outer radius boundary */
    radius = r_min + (r_cell+1)*dr;
    coeff_B = 2 * (x*nx + y*ny);
    coeff_C = x*x + y*y - radius*radius;
    delta = coeff_B*coeff_B - 4*coeff_C;
    s_temp = (-coeff_B + sqrt(delta)) / 2.0;
    s = MIN(s, s_temp);
    
    return s;
    
}

void move_one_step(int i){
    
    double x, y, r, theta, nx, ny, phi;
    int r_cell, theta_cell;
    double s, den_cell, tau, zeta, tau_rand;
    
    x = rays[i]->x;
    y = rays[i]->y;
    nx = rays[i]->nx;
    ny = rays[i]->ny;
    r = rays[i]->r;
    theta = rays[i]->theta;
    r_cell = rays[i]->r_cell;
    theta_cell = rays[i]->theta_cell;
    den_cell = *(density+r_cell*N_theta+theta_cell);

    s = find_step_size(i);
    tau = den_cell * KAPPA * s;
    
    zeta = rand_float();
    tau_rand = -log(1-zeta);
    
    /* able to left this cell */
    if (tau_rand >= tau) {
        
        s = s + dr*SMALL_NUMBER; /* to avoid boundary */
        x = x + s*nx;
        y = y + s*ny;
        r = sqrt(x*x + y*y);
        theta = acos(x/r);
        
        if (y<0) {
            theta = 2*PI - theta;
        }
        
        if (r<r_min OR r>r_max OR theta<theta_min OR theta>theta_max) {
            rays[i]->out_flag = YES;
        }
        
        r_cell = (int)((r - r_min) / dr);
        theta_cell = (int)((theta - theta_min) / dtheta);
        
        /* up date the info */
        rays[i]->x = x;
        rays[i]->y = y;
        rays[i]->r = r;
        rays[i]->theta = theta;
        rays[i]->r_cell = r_cell;
        rays[i]->theta_cell = theta_cell;
        
    }
    
    /* a scatter happens */
    if (tau_rand < tau) {
        
        s = tau_rand / (KAPPA * den_cell);
        x = x + s*nx;
        y = y + s*ny;
        r = sqrt(x*x + y*y);
        theta = acos(x/r);
        
        phi = 2 * PI * rand_float();
        nx = cos(phi);
        ny = sin(phi);
        
        /* update the info */
        rays[i]->x = x;
        rays[i]->y = y;
        rays[i]->r = r;
        rays[i]->theta = theta;
        rays[i]->nx = nx;
        rays[i]->ny = ny;
        
    }
    
    /*
    printf("x=%f, y=%f, r=%f, theta=%f, nx=%f, ny=%f\n",
           rays[i]->x, rays[i]->y, rays[i]->r, rays[i]->theta,
           rays[i]->nx, rays[i]->ny);
    */
    
}

void free_ray(void){
    
    int i;
    free(rays[0]);
    free(rays);
    
}
