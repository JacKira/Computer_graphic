#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
void swap_f(double *a, double *b);
double **create_z_buffer(size_t N, size_t M); //Allocate z_buffer and fill
void z_buffer_free(double **A, size_t N); //free z_buffer memory
void BuildNormal (Vec3 p[3], Vec3 normal);
void transform_normal(Vec3 normal, int phi);
void rotate(Vec3 p, int phi);
void perspec(Vec3 p, double c);

// Triangle rasterization algorithm using a scan line
void triangle(tgaImage *image,
              int x0, int y0, double z0, double u0, double v0,
              int x1, int y1, double z1, double u1, double v1,
              int x2, int y2, double z2, double u2, double v2, double **zbuffer,
              double intensity, Model *model);

// Drawing with simple lighting model
void render(tgaImage *image, Model *model);

int main(int argc, char **argv)
{
    int rv = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <objfile> <outfile>\n", argv[0]);
        return -1;
    }

    Model *model = NULL;
    tgaImage *image = NULL;

    do {
        model = loadFromObj(argv[1]);
        if (!model) {
            perror("loadFromObj");
            rv = -1;
            break;
        }

        image = tgaNewImage(4000, 4000, RGB);
        if (!image) {
            perror("tgaNewImage");
            rv = -1;
            break;
        }
        if(loadDiffuseMap(model, argv[3]) == 0){ //load texture
            perror("loadDiffuseMap");
            rv = -1;
            break;
        }
        // Drawing
        render(image, model);

        if (-1 == tgaSaveToFile(image, argv[2])) {
            perror("tgaSaveToFile");
            rv = -1;
            break;
        }
    } while (0);

    if (model) {
        freeModel(model);
    }
    if (image) {
        tgaFreeImage(image);
    } 
    return rv;
}

void render(tgaImage *image, Model *model)
{
    int face, vert, phi = 20;
    double c = 7;
    int h = image->height;
    int w = image->width;
    Vec3 *p[3], *vt[3], normal, light = {0, 0, -1};
    Vec3 screen_coords[3];
    double I;
    int scale = 1;
    double **z_buffer = create_z_buffer(w, h);

    for (face = 0; face < model->nface; ++face) {
        // Translating into screen coordinates
        for (vert = 0; vert < 3; ++vert) {
            p[vert] = getVertex(model, face, vert);
            screen_coords[vert][0] = (*p[vert])[0];
            screen_coords[vert][1] = (*p[vert])[1];
            screen_coords[vert][2] = (*p[vert])[2];
            rotate(screen_coords[vert], phi);
            perspec(screen_coords[vert], c);
        }
        BuildNormal(screen_coords, normal);
        for (vert = 0; vert < 3; ++vert) {
            screen_coords[vert][0] = (screen_coords[vert][0] / scale + 1.0) * w / 2;
            screen_coords[vert][1] = (screen_coords[vert][1] / scale + 1.0) * h / 2;
            vt[vert] = getDiffuseUV(model, face, vert);
        }
        // Computing light intensity
        I = (light[0] * normal[0] + light[1] * normal[1] + light[2] * normal[2]);
        if(I <= 0){
            // Fill the triangle
            triangle(image, screen_coords[0][0], screen_coords[0][1],  screen_coords[0][2], (*vt[0])[0], (*vt[0])[1],
                            screen_coords[1][0], screen_coords[1][1],  screen_coords[1][2], (*vt[1])[0], (*vt[1])[1],
                            screen_coords[2][0], screen_coords[2][1],  screen_coords[2][2], (*vt[2])[0], (*vt[2])[1], 
                            z_buffer, fabs(I), model);
        }
    }
    tgaFlipVertically(image); 
    z_buffer_free(z_buffer, h);
}

void triangle(tgaImage *image,
              int x0, int y0, double z0, double u0, double v0,
              int x1, int y1, double z1, double u1, double v1,
              int x2, int y2, double z2, double u2, double v2,
              double **z_buffer, double intensity, Model *model)
{
    tgaColor color;  //define color varible
    // Sorting vertices by 'y' coord
    if (y0 > y1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
        swap_f(&z0, &z1);
        swap_f(&u0, &u1);
        swap_f(&v0, &v1);
    }
    if (y0 > y2) {
        swap(&x0, &x2);
        swap(&y0, &y2);
        swap_f(&z0, &z2);
        swap_f(&u0, &u2);
        swap_f(&v0, &v2);
    }
    if (y1 > y2) {
        swap(&x1, &x2);
        swap(&y1, &y2);
        swap_f(&z1, &z2);
        swap_f(&u1, &u2);
        swap_f(&v1, &v2);
    }

    double t_less, t_more, tz, z, z_a, z_b;
    double x_left, x_right, x;
    double u_a, u_b, v_a, v_b;
    int r, g, b;
    Vec3 uv;
    for(int y = y0; y <= y2; y++){
        // solve parametric equalation [x(t), y(t), z(t), u(t), v(t)] and find t by y(t)
        t_more = (double)(y - y0) / (y2 - y0);
        if(y > y1){                                    //lower traingle
            t_less = (double)(y - y1) / (y2 - y1);
            x_left = x2 * t_less + x1 * (1 - t_less);
            z_a = z2 * t_less + z1 * (1 - t_less);
            u_a = u2 * t_less + u1 * (1 - t_less);
            v_a = v2 * t_less + v1 * (1 - t_less);
        }
        else{                                         //upper triangle
            t_less = (double)(y - y0) / (y1 - y0);
            x_left = x1 * t_less + x0 * (1 - t_less);
            z_a = z1 * t_less + z0 * (1 - t_less);
            u_a = u1 * t_less + u0 * (1 - t_less);
            v_a = v1 * t_less + v0 * (1 - t_less);
        }
        x_right = x2 * t_more + x0 * (1 - t_more); 
        z_b = z2 * t_more + z0 * (1 - t_more);
        u_b = u2 * t_more + u0 * (1 - t_more);
        v_b = v2 * t_more + v0 * (1 - t_more);
        // interval boundaries
        if (x_left > x_right){
            swap_f(&x_left, &x_right);
        }      
        if (z_a > z_b){
            swap_f(&z_a, &z_b);
        }
        if (u_a > u_b){
            swap_f(&u_a, &u_b);
        }
        if (v_a > v_b){
            swap_f(&v_a, &v_b);
        }
        for(x = x_left; x <= x_right; x++){
            tz = (double)(x - x_left) / (x_right - x_left);
            // calculate value of z_buffer
            z = z_b * tz + z_a * (1 - tz);
            if(z > z_buffer[(int)x][y]){
                z_buffer[(int)x][y] = z;
                uv[0] = u_b * tz + u_a * (1 - tz);
                uv[1] = v_b * tz + v_a * (1 - tz);
                color = getDiffuseColor(model, &uv);
                r = intensity * Red(color);
                g = intensity * Green(color);
                b = intensity * Blue(color);
                tgaSetPixel(image, x, y, tgaRGB(r, g, b));
            }
        }
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

void swap_f(double *a, double *b){
    double t = *a;
    *a = *b;
    *b = t;    
}

int sign(int a) {
    return (a < 0) ? -1 : 1;
}

double ** create_z_buffer(size_t N, size_t M)
{
    double **A = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
        A[i] = (double *)malloc(M*sizeof(double));
        for(int j = 0; j < M; j++){
            A[i][j] = -1e10;
        }
    }
    return A;
}

void z_buffer_free(double **A, size_t N)
{
    for(int i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
}

void BuildNormal(Vec3 p[3], Vec3 normal)
{
    
    Vec3 v1, v2;
    double length, ilength;

    // find v1 (v1 = p2 - p1)
    v1[0] = p[1][0] - p[0][0]; // x
    v1[1] = p[1][1] - p[0][1]; // y
    v1[2] = p[1][2] - p[0][2]; // z
    
    // find v2 (v2 = p3 - p1)
    v2[0] = p[2][0] - p[0][0]; // x
    v2[1] = p[2][1] - p[0][1]; //y
    v2[2] = p[2][2] - p[0][2]; // z
         
    // find normal [v1*v2]
    normal[0] = v1[1] * v2[2] - v1[2] * v2[1]; // y1*z2 - z1*y2 = Xn
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2]; // z1*x2 - x1*z2 = Yn
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0]; // x1*y2 - y1*x2 = Zn

    // normalizing vector
    // find length of the normal
    length = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    length = sqrt(length); 

    ilength = 1/length;
    
    normal[0] *= ilength;
    normal[1] *= ilength;
    normal[2] *= ilength;
}

void rotate(Vec3 p, int angle){
    double phi = angle * 3.14 / 180;
    double T[4][4] = {{cos(phi),      0,   sin(phi), 0},
                     { 0,  1,   0, 0},
                     {  -sin(phi),  0,   cos(phi), 0},
                     {   0,        0,        0, 1}};
    p[0] = T[0][0] * p[0] + T[1][0] * p[1] + T[2][0] * p[2] + T[3][0];
    p[1] = T[0][1] * p[0] + T[1][1] * p[1] + T[2][1] * p[2] + T[3][1];
    p[2] = T[0][2] * p[0] + T[1][2] * p[1] + T[2][2] * p[2] + T[3][2];

}

void perspec(Vec3 p, double c){
    p[0] = p[0] / (1 - p[2] / c);
    p[1] = p[1] / (1 - p[2] / c);
    p[2] = p[2] / (1 - p[2] / c);
}

void transform_normal(Vec3 p, int angle){
    double phi = angle * 3.14 / 180;
    double T[4][4] = {{cos(phi) ,      sin(phi) * sin(phi),   -sin(phi), 0},
                     { 0,  1,    0, 0},
                     {  sin(phi), -cos(phi)*sin(phi),    cos(phi), 0},
                     {   0,        0,        0, 1}};
    p[0] = T[0][0] * p[0] + T[1][0] * p[1] + T[2][0] * p[2] + T[3][0];
    p[1] = T[0][1] * p[0] + T[1][1] * p[1] + T[2][1] * p[2] + T[3][1];
    p[2] = T[0][2] * p[0] + T[1][2] * p[1] + T[2][2] * p[2] + T[3][2];
    double length = 1; // sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);                 
    p[0] *= length;
    p[1] *= length;
    p[2] *= length;
}   
