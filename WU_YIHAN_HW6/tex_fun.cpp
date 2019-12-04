/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
GzColor * minMap;
int xs, ys;
int reset = 1;
float max1 = 1;
float max2 = 1;
float max3 = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color, float scaleSizeU, float scaleSizeV, float scaleSizeUm, float scaleSizeVm)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  color[0] = 0;
  color[1] = 0;
  color[2] = 0;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
	  if (max2 > (image[i][RED] + image[i][GREEN] + image[i][BLUE]))
		  max2 = image[i][RED] + image[i][GREEN] + image[i][BLUE];
      }



    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
  if (u < 0) {
	  u = 0;
  }
  else if (u > 1) {
	  u = 1;
  }
  if (v < 0) {
	  v = 0;
  }
  else if (v > 1) {
	  v = 1;
  }

  float Uc, Vc;
  Uc = u * (xs - 1);
  Vc = v * (ys - 1);

  float minUc;
  float minVc;
  float maxUc;
  float maxVc;

  if (scaleSizeU > 0 || scaleSizeUm < 0) {
	  minUc = max(0, Uc + scaleSizeUm);
	  maxUc = min(xs - 1, Uc + scaleSizeU);
  }
  else {
	  minUc = max(0, Uc + scaleSizeU);
	  maxUc = min(xs - 1, Uc + scaleSizeUm);
  }
  if (scaleSizeV > 0 || scaleSizeVm < 0) {
	  minVc = max(0, Vc + scaleSizeVm);
	  maxVc = min(ys - 1, Vc + scaleSizeV);
  }
  else {
	  minVc = max(0, Vc + scaleSizeV);
	  maxVc = min(ys - 1, Vc + scaleSizeVm);
  }

  float a = minUc, b = minVc;
  int cnt = 0;


  /*int Ax = floor(Uc);
  int Ay = floor(Vc);
  int Bx = ceil(Uc);
  int By = floor(Vc);
  int Cx = ceil(Uc);
  int Cy = ceil(Vc);
  int Dx = floor(Uc);
  int Dy = ceil(Vc);
  GzColor pixelA, pixelB, pixelC, pixelD;
  for (int i = 0; i < 3; i++) {
	  pixelA[i] = image[Ax + Ay * (xs)][i];
	  pixelB[i] = image[Bx + By * (xs)][i];
	  pixelC[i] = image[Cx + Cy * (xs)][i];
	  pixelD[i] = image[Dx + Dy * (xs)][i];
  }
  float s = Uc - Ax;
  float t = Vc - Ay;

  for (int i = 0; i < 3; i++) {
	  color[i] = s * t * pixelC[i] + (1 - s) * t * pixelD[i] + s * (1 - t) * pixelB[i] + (1 - s) * (1 - t) * pixelA[i];
  }*/
  while (a <= maxUc) {
	  while (b <= maxVc) {

		  int Ax = floor(a);
		  int Ay = floor(b);
		  int Bx = ceil(a);
		  int By = floor(b);
		  int Cx = ceil(a);
		  int Cy = ceil(b);
		  int Dx = floor(a);
		  int Dy = ceil(b);
		  GzColor pixelA, pixelB, pixelC, pixelD;
		  for (int i = 0; i < 3; i++) {
			  pixelA[i] = image[Ax + Ay * (xs)][i];
			  pixelB[i] = image[Bx + By * (xs)][i];
			  pixelC[i] = image[Cx + Cy * (xs)][i];
			  pixelD[i] = image[Dx + Dy * (xs)][i];
		  }
		  float s = a - Ax;
		  float t = b - Ay;

		  float temp[3];

		  for (int i = 0; i < 3; i++) {
			  temp[i] = s * t * pixelC[i] + (1 - s) * t * pixelD[i] + s * (1 - t) * pixelB[i] + (1 - s) * (1 - t) * pixelA[i];
			  color[i] += s * t * pixelC[i] + (1 - s) * t * pixelD[i] + s * (1 - t) * pixelB[i] + (1 - s) * (1 - t) * pixelA[i];
		  }
		  if (scaleSizeV == 62)
			  color[0] = color[0];
		  if (temp[0] + temp[1] + temp[2] < max3)
			  max3 = temp[0] + temp[1] + temp[2];
		  b++;
		  cnt++;
		  //break;
	  }
	  a++;
	  }
  for (int i = 0; i < 3; i++) {
	  color[i] /= cnt;
	  if (color[0] + color[1] + color[2] < max1)
		  max1 = color[0] + color[1] + color[2];
  }
  if (cnt == 0)
	  cnt = cnt;
  

  
  return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color, int scaleSizeX, int scaleSizeY)
{
	float CX[2];
	CX[0] = u;
	CX[1] = v;
	int N = 100;
	float limit = pow(10, 10);
	for (int i = 0; i < 100; i++) {
		float temp1 = CX[0] * CX[0] - CX[1] * CX[1] - 0.12375;
		float temp2 = 2 * CX[0] * CX[1] + 0.56805;
		if (abs(temp1) > limit || abs(temp2) > limit)
			break;
		CX[0] = temp1;
		CX[1] = temp2;
	}
	
	float length;
	length = sqrt(CX[0] * CX[0] + CX[1] * CX[1]);

	color[0] = color[1] = color[2] = length / pow(10, 9);


	return GZ_SUCCESS;
}

int ptex_fun_checkboard(float u, float v, GzColor color, int scaleSizeX, int scaleSizeY) {
	u *= 8;
	v *= 8;
	int x = int(u);
	int y = int(v);
	if ((x + y) % 2 == 0) {
		color[0] = color[1] = color[2] = 0;
	}
	else {
		color[0] = color[1] = color[2] = 0.5;
	}

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}




/*
if (reset) {         
	xs = 1024;
	ys = 1024;
	image = new GzColor[xs * ys * 3];
	float CX[2];
	CX[0] = 0;
	CX[1] = 0;
	int N = 500;
	float length;
	for (int i = 0; i < N; i++) {
		float temp1 = CX[0] * CX[0] - CX[1] * CX[1] - 0.12375;
		float temp2 = 2 * CX[0] * CX[1] + 0.56805;
		length = sqrt(temp1 * temp1 + temp2 * temp2);
		if (length > 2)
			break;
		CX[0] = temp1;
		CX[1] = temp2;
		int x = round(max(min(CX[0] + 1, 1.0), 0.0) * 1024);
		int y = round(max(min(CX[1], 1.0), 0.0) * 1024);
		image[x + y * xs][0] = image[x + y * xs][2] = image[x + y * xs][1] = 1;
	}
	reset = 0;
}
image[100]

if (u < 0) {
	u = 0;
}
else if (u > 1) {
	u = 1;
}
if (v < 0) {
	v = 0;
}
else if (v > 1) {
	v = 1;
}

float Uc, Vc;
Uc = u * (xs - 1);
Vc = v * (ys - 1);

int Ax = floor(Uc);
int Ay = floor(Vc);
int Bx = ceil(Uc);
int By = floor(Vc);
int Cx = ceil(Uc);
int Cy = ceil(Vc);
int Dx = floor(Uc);
int Dy = ceil(Vc);
GzColor pixelA, pixelB, pixelC, pixelD;
for (int i = 0; i < 3; i++) {
	pixelA[i] = image[Ax + Ay * (xs)][i];
	pixelB[i] = image[Bx + By * (xs)][i];
	pixelC[i] = image[Cx + Cy * (xs)][i];
	pixelD[i] = image[Dx + Dy * (xs)][i];
}
float s = Uc - Ax;
float t = Vc - Ay;

for (int i = 0; i < 3; i++) {
	color[i] = s * t * pixelC[i] + (1 - s) * t * pixelD[i] + s * (1 - t) * pixelB[i] + (1 - s) * (1 - t) * pixelA[i];
}*/