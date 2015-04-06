/**
* OS specific
**/
#define OS_UNIX 1
#define OS_WIN 2

#ifdef USE_UNIX
#   define OS OS_UNIX
#else
#   define OS OS_WIN
#endif

#if OS == OS_WIN
#   define _CRT_SECURE_NO_WARNINGS
#   include <io.h>
#   include <fcntl.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "bmp.h"


/**
 * prototypes
 **/
static void ShowUsage(void);
static void CheckOption(char *option, int argc, int minargc);
static void SetImageMask(Image *img, Image *mask);
static void SetBinaryIOMode(void);
static void ReadCorrespondences(char *fileName, int& numSamples, Line **L0, Line **L1);


/**
 * main
 **/
int main(int argc, char **argv)
{
    Image *img = NULL;

    // setup
    SetBinaryIOMode();

    // first argument is program name
    argv++, argc--;

    // look for help
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-help")) {
            ShowUsage();
        }
    }

    // no argument case
    if (argc == 0) {
        img = BMPReadImage(stdin);
    }

    // parse arguments
    while (argc > 0)
    {
        if (**argv == '-')
        {
	  /*
            if (!strcmp(*argv, "-noise"))
            {
                double factor;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                factor = atof(argv[1]);
                img->AddNoise(factor);
                argv += 2, argc -= 2;
            }
	  */
            if (!strcmp(*argv, "-brightness"))
            {
                double factor;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                factor = atof(argv[1]);
                img->Brighten(factor);
                argv += 2, argc -=2;
            }

            else if (!strcmp(*argv, "-contrast"))
            {
                double factor;
                CheckOption(*argv, argc, 2);
                    if (img == NULL)
                        img = BMPReadImage(stdin);

                factor = atof(argv[1]);
                img->ChangeContrast(factor);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-saturation"))
            {
                double factor;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                factor = atof(argv[1]);
                img->ChangeSaturation(factor);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-gamma"))
            {
                double factor;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                factor = atof(argv[1]);
                img->ChangeGamma(factor);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-crop"))
            {
                int x, y, w, h;
                CheckOption(*argv, argc, 5);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                x = atoi(argv[1]);
                y = atoi(argv[2]);
                w = atoi(argv[3]);
                h = atoi(argv[4]);

                Image *dst = img->Crop(x, y, w, h);
                delete img;
                img = dst;

                argv += 5, argc -= 5;
            }
	    
	    /*
            else if (!strcmp(*argv, "-extractChannel"))
            {
                int channel;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                channel = atoi(argv[1]);
                img->ExtractChannel(channel);
                argv += 2, argc -= 2;
            }
	    */

            else if (!strcmp(*argv, "-quantize"))
            {
                int nbits;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                nbits = atoi(argv[1]);
                img->Quantize(nbits);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-randomDither"))
            {
                int nbits;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                nbits = atoi(argv[1]);
                img->RandomDither(nbits);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-composite"))
            {
                FILE *fp;
                Image *bottomMask, *top, *topMask;

                CheckOption(*argv, argc, 4);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                fp = fopen(argv[1], "rb");
                assert(fp != NULL);
                bottomMask = BMPReadImage(fp);
                SetImageMask(img, bottomMask);

                fp = fopen(argv[2], "rb");
                assert(fp != NULL);
                top = BMPReadImage(fp);

                fp = fopen(argv[3], "rb");
                assert(fp != NULL);
                topMask = BMPReadImage(fp);
                SetImageMask(top, topMask);

                ImageComposite(img, top, img);
                argv += 4, argc -= 4;
            }

            else if (!strcmp(*argv, "-blur"))
            {
                int n;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                n = atoi(argv[1]);
                img->Blur(n);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-sharpen"))
            {
                if (img == NULL)
                    img = BMPReadImage(stdin);

                img->Sharpen();
                argv++ , argc--;
            }

            else if (!strcmp(*argv, "-edgeDetect"))
            {
                int factor;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

		factor = atoi(argv[1]) ;
                img->EdgeDetect(factor);
                argv+= 2, argc -= 2;
            }
	    /*
            else if (!strcmp(*argv, "-orderedDither"))
            {
                int nbits;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                nbits = atoi(argv[1]);
                img->OrderedDither(nbits);
                argv += 2, argc -= 2;
            }
	    */
            else if (!strcmp(*argv, "-FloydSteinbergDither"))
            {
                int nbits;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                nbits = atoi(argv[1]);
                img->FloydSteinbergDither(nbits);
                argv += 2, argc -= 2;
            }

            else if (!strcmp(*argv, "-size"))
            {
                CheckOption(*argv, argc, 3);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                int sizex = atoi(argv[1]);
                int sizey = atoi(argv[2]);

                Image *dst = img->Scale(sizex, sizey);
                delete img;

                img = dst;
                argv += 3, argc -= 3;
            }

            else if (!strcmp(*argv, "-shift"))
            {
                CheckOption(*argv, argc, 3);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                double sx = atof(argv[1]);
                double sy = atof(argv[2]);

                img->Shift(sx, sy);
                argv += 3, argc -= 3;
            }

	    /*
            else if (!strcmp(*argv, "-rotate"))
            {
                double angle;
                Image *dst;
                CheckOption(*argv, argc, 2);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                angle = atof(argv[1]);
                dst = img->Rotate(angle);
                delete img;
                img = dst;
                dst = NULL;
                argv += 2, argc -= 2;
            }
	    */

            else if (!strcmp(*argv, "-fun"))
            {
                if (img == NULL)
                    img = BMPReadImage(stdin);

                img->Fun();
                argv++, argc--;
            }

            else if (!strcmp(*argv, "-morph"))
            {
                CheckOption(*argv, argc, 4);
                if (img == NULL)
                    img = BMPReadImage(stdin);

                // read destination image
                FILE *fp = fopen(argv[1], "rb");
                assert(fp != NULL);

                Image *dstImage = BMPReadImage(fp);
                fclose(fp);

                // read correspondences
                int numLines;
                Line *L0, *L1;
                ReadCorrespondences(argv[2], numLines, &L0, &L1);

                // read t
                double t = atof(argv[3]);

                // warp
                Image *warp = ImageMorph(img, dstImage, numLines, L0, L1, t);

                delete img;
                delete dstImage;
                img = warp;

                argv += 4, argc -= 4;
            }

            else if (!strcmp(*argv, "-sampling"))
            {
                if (img == NULL)
                    img = BMPReadImage(stdin);

                int method;
                CheckOption(*argv, argc, 2);
                method = atoi(argv[1]);
                img->SetSamplingMethod(method);
                argv += 2, argc -= 2;
            }

            else
            {
                fprintf(stderr, "image: invalid option: %s\n", *argv);
                ShowUsage();
            }
        } else {
            fprintf(stderr, "image: invalid option: %s\n", *argv);
            ShowUsage();
        }
    }

    // write image
    if (img != NULL)
    {
        BMPWriteImage(img, stdout);
        delete img;
    }

    // done
    return EXIT_SUCCESS;
}



/**
 * ShowUsage
 **/
static char options[] =
"-help\n"
"-brightness <factor>\n"
"-contrast <factor>\n"
"-saturation <factor>\n"
"-gamma <gamma>\n"
"-crop <x> <y> <width> <height>\n"
"-quantize <nbits>\n"
"-randomDither <nbits>\n"
"-FloydSteinbergDither <nbits>\n"
"-blur <maskSize>\n"
"-sharpen\n"
"-edgeDetect <threshold>\n"
"-size <sx> <sy>\n"
"-shift <tx> <ty>\n"
"-fun\n"
"-composite <bottomMaskFileName> <topImageFileName> <topMaskFileName>\n"
"-morph <dstImageFileName> <correspondencesFileName> <t>\n"
"-sampling <method no>\n"
;

static void ShowUsage(void)
{
    fprintf(stderr, "Usage: image [-option [arg ...] ...]\n");
    fprintf(stderr, options);
    exit(EXIT_FAILURE);
}



/**
 * CheckOption
 **/
static void CheckOption(char *option, int argc, int minargc)
{
    if (argc < minargc)
    {
        fprintf(stderr, "Too few arguments for %s\n", option);
        ShowUsage();
    }
}



/**
 * SetImageMask
 *
 * Reads the alpha channel of "img", from the blue channel of "mask".
 * This gets around the fact that we are using 24-bit files
 * that don't store alpha.
 **/
static void SetImageMask(Image *img, Image *mask)
{
    assert(img && mask);
    assert(img->width == mask->width);
    assert(img->height == mask->height);

    for (int i=0; i<img->NumPixels(); i++)
        img->pixels[i].a = 255 - mask->pixels[i].b;
}



/**
 * SetBinaryIOMode
 *
 * In WindowsNT by default stdin and stdout are opened as text files.
 * This code opens both as binary files.
 **/
static void SetBinaryIOMode(void)
{
#if OS == OS_WIN
    int result;

    result = _setmode(_fileno(stdin), _O_BINARY);
    if (result == -1)
    {
        perror("Cannot set stdin binary mode");
        exit(EXIT_FAILURE);
    }

    result = _setmode(_fileno(stdout), _O_BINARY);
    if (result == -1)
    {
        perror("Cannot set stdout binary mode");
        exit(EXIT_FAILURE);
    }
#endif
}


/**
 * Read correspondences for morph.
 *
 * File format is as follows:
 *    [numLines]
 *    [line 0] [line 0]'
 *    ....
 *    [line NumLines-1] [line NumLines-1]'
 *
 * where [line] is:
 *    [point P] [point Q]
 *
 * where [point] is:
 *    [x y]
 *
 * whitespace doesn't count.
 **/
static void ReadCorrespondences(char *fileName, int& numSamples, Line **L0, Line **L1)
{
    FILE *fp = fopen(fileName, "rt");
    assert(fp);

    int res = fscanf(fp, "%d", &numSamples);
    assert(res==1 && numSamples>0);

    *L0 = new Line[numSamples];
    *L1 = new Line[numSamples];
    for(int i=0; i<numSamples; ++i)
    {
        float px, py, qx, qy;

        res = fscanf(fp, "%f %f %f %f", &px, &py, &qx, &qy);
        assert(res == 4);
        (*L0)[i].p.Set(px,py);
        (*L0)[i].q.Set(qx,qy);

        res = fscanf(fp, "%f %f %f %f", &px, &py, &qx, &qy);
        assert(res == 4);
        (*L1)[i].p.Set(px,py);
        (*L1)[i].q.Set(qx,qy);
    }
}
