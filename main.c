#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
void swap_f(double *a, double *b);
int abs(int a);
int sign(float a);
void BuildNormal (Vec3 *p[3], double normal[3]);


void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);


void triangle(tgaImage *image,
              int x0, int y0,
              int x1, int y1,
              int x2, int y2,
              tgaColor color);

void meshgrid(tgaImage *image, Model *model);

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

        image = tgaNewImage(8000, 8000, RGB);
        if (!image) {
            perror("tgaNewImage");
            rv = -1;
            break;
        }

        render(image, model);

        if (-1 == tgaSaveToFile(image, argv[2])) {
            perror("tgaSateToFile");
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

//Fill triangles in meshgrid by random color
void render(tgaImage *image, Model *model)
{
    int face, vert;
    int h = image->height;
    int w = image->width;
    int r, g, b;
    double  I, light[3] = {0, 0, -1};
    Vec3 *p[3];
    double normal[3];
    int screen_coords[3][2];
    for (face = 0; face < model->nface; ++face) {
        for (vert = 0; vert < 3; ++vert) {
            p[vert] = getVertex(model, face, vert);
            screen_coords[vert][0] = ((*p[vert])[2]/1000 + 1.0)*w/2;
            screen_coords[vert][1] = ((*p[vert])[1]/1000 + 1.0)*h/2;
        }
        BuildNormal(p, normal);
        I = light[0] * normal[0] + light[1] * normal[1] + light[2] * normal[2];
        if (I < 0){
            r = abs((int)(I * 255)) % 256;
            g = abs((int)(I * 255)) % 256;
            b = abs((int)(I * 255)) % 256; 
            triangle(image, screen_coords[0][0], screen_coords[0][1],
                            screen_coords[1][0], screen_coords[1][1],
                            screen_coords[2][0], screen_coords[2][1],
                            tgaRGB(r, g, b));
        }
        
    }
    tgaFlipVertically(image);    
}


void triangle(tgaImage *image,
              int x0, int y0,
              int x1, int y1,
              int x2, int y2,
              tgaColor color)
{
    /* Sort vertices by y coord */
    if (y0 > y1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }
    if (y0 > y2) {
        swap(&x0, &x2);
        swap(&y0, &y2);
    }
    if (y1 > y2) {
        swap(&x1, &x2);
        swap(&y1, &y2);
    }

    /* Compute x coord deltas */
    double dx01 = 0, dx02 = 0, dx12 = 0;
    if (y0 != y1) {
        dx01 = x1 - x0;
        dx01 /= y1 - y0;
    }
    if (y0 != y2) {
        dx02 = x2 - x0;
        dx02 /= y2 - y0;
    }
    if (y1 != y2) {
        dx12 = x2 - x1;
        dx12 /= y2 - y1;
    }
    double _dx02 = dx02;

    if (dx01 > dx02) {
        swap_f(&dx01, &dx02);
    }

    if (dx12 > _dx02) {
        swap_f(&dx12, &_dx02);
    }

    /* Fill the triangle from up to down by Y*/
    int y, x;
    double xleft = x0, xright = x0;
    for (y = y0; y <= y2; ++y) {
        for (x = xleft; x <= xright; ++x) {
            tgaSetPixel(image, x, y, color);
        }
        xleft += (y < y1) ? dx01 : _dx02;
        xright += (y < y1) ? dx02 : dx12;
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

int abs(int a) {
    return (a >= 0) ? a : -a;
}

int sign(float a){
    return (a > 0) ? 1 : -1;
}

void BuildNormal(Vec3 *p[3], double normal[3])
{
    
    double v1[3], v2[3];
    double length, ilength;

    // найдем вектор v1 задающий одно ребро треугольника (v1 = p2 - p1)
    v1[0] = (*p[1])[0] - (*p[0])[0]; // x
    v1[1] = (*p[1])[1] - (*p[0])[1]; // y
    v1[2] = (*p[1])[2] - (*p[0])[2]; // z
    
    // найдем вектор v2 задающий другое ребро треугольника (v2 = p3 - p1)
    v2[0] = (*p[2])[0] - (*p[0])[0]; // x
    v2[1] = (*p[2])[1] - (*p[0])[1]; //y
    v2[2] = (*p[2])[2] - (*p[0])[2]; // z
         
    // векторное произведение [v1*v2]
    normal[0] = v1[1] * v2[2] - v1[2] * v2[1]; // y1*z2 - z1*y2 = Xn
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2]; // z1*x2 - x1*z2 = Yn
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0]; // x1*y2 - y1*x2 = Zn
    
    //
    // нормализуем вектор
    //
    
    // найдем длину вектора
    length = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    length = sqrt(length); 

    ilength = 1/length;
    // возьмем обратную величину, чтоб не делать потом лишних делений
    
    normal[0] *= ilength;
    normal[1] *= ilength;
    normal[2] *= ilength;
    
}
