#include <stdio.h>

#include "tga.h"
#include "model.h"

void swap(int *a, int *b);
int abs(int a);
int sign(float a);
/*
* Using Digital Bresenhame algorihm
* to draw interval connecting (x0, y0) with (x1, y1)
* on image using color
*/
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);


/*
* Draw model using meshgrid
*/
void meshgrid(tgaImage *image, Model *model);

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

        image = tgaNewImage(3000, 3000, RGB);
        if (!image) {
            perror("tgaNewImage");
            rv = -1;
            break;
        }

        meshgrid(image, model);

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

//Draw msehgrid
void meshgrid(tgaImage *image, Model *model)
{
    int face, vert, i;
    int h = image->height;
    int w = image->width;
    tgaColor white = tgaRGB(255, 255, 255);
    Vec3 *p[3];
    for (face = 0; face < model->nface; ++face) {
        for (vert = 0; vert < 3; ++vert) {
            p[vert] = getVertex(model, face, vert);
        }
        for (i = 0; i < 3; ++i) {
            line(image,
                ((*p[i])[0]/1000 + 1.0)*w/2,         ((*p[i])[1]/1000 + 1.0)*h/2,
                ((*p[(i + 1)%3])[0]/1000 + 1.0)*w/2, ((*p[(i + 1)%3])[1]/1000 + 1.0)*h/2,
                white);
        }
    }
    tgaFlipVertically(image);
}

//Digital Bresenhame algorihm
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color)
{
    int steep = 0;
    if (abs(y1 - y0) > abs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }
    
    int x;
    int y, dx = x1 - x0, dy = y1 - y0;
    int e, de = 2*abs(dy);
    for (x = x0, y = y0, e = 0; x <= x1; ++x) {
        if (steep) {
            tgaSetPixel(image, y, x, color);
        } else {
            tgaSetPixel(image, x, y, color);
        }
        e += de;
        if (e > dx){
            y += sign(dy);
            e -= 2*dx;
        }
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int abs(int a) {
    return (a >= 0) ? a : -a;
}

int sign(float a){
    return (a > 0) ? 1 : -1;
}
