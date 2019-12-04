/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <algorithm>

int maxScale = 0;
#define PI (float) 3.14159265358979323846

float AAFilter[AAKERNEL_SIZE][3]{
-0.52, 0.38, 0.128, 		0.41, 0.56, 0.119,		0.27, 0.08, 0.294,
-0.17, -0.29, 0.249,		0.58, -0.55, 0.104,		-0.31, -0.71, 0.106
};


void copyMatrix(GzMatrix a, GzMatrix b) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			a[i][j] = b[i][j];
}

void multiMatrix(GzMatrix a, GzMatrix c, GzMatrix b) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			a[i][j] = c[i][0] * b[0][j] + c[i][1] * b[1][j] + c[i][2] * b[2][j] + c[i][3] * b[3][j];
}

float multiVector(GzCoord vector1, GzCoord vector2) {
	return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
}

void norm(GzCoord v) {

	float n = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));

	v[0] = v[0] / n;
	v[1] = v[1] / n;
	v[2] = v[2] / n;
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	degree *= (PI / 180);

	mat[0][0] = 1;
	mat[1][1] = cos(degree);
	mat[1][2] = -sin(degree);
	mat[2][1] = sin(degree);
	mat[2][2] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	degree *= (PI / 180);

	mat[1][1] = 1;
	mat[0][0] = cos(degree);
	mat[0][2] = sin(degree);
	mat[2][0] = -sin(degree);
	mat[2][2] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/

	degree *= (PI / 180);

	mat[2][2] = 1;
	mat[0][0] = cos(degree);
	mat[0][1] = -sin(degree);
	mat[1][0] = sin(degree);
	mat[1][1] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	mat[0][0] = 1;
	mat[1][1] = 1;
	mat[2][2] = 1;
	mat[3][3] = 1;
	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */
	xres = xRes;
	yres = yRes;
	for (int k = 0; k < AAKERNEL_SIZE + 1; k++) {
		pixelbuffer[k] = new GzPixel[xres * yres];
	}
	framebuffer = new char[xRes * yRes * 3];


	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/

	matlevel = 0;
	matNormLevel = 0;
	numlights = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			Xsp[i][j] = 0;
	}

	Xsp[0][0] = Xsp[0][3] = xres / 2;
	Xsp[1][1] = -yres / 2;
	Xsp[1][3] = yres / 2;
	Xsp[2][2] = INT_MAX;
	Xsp[3][3] = 1;

	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;

	m_camera.lookat[0] = 0.0;
	m_camera.lookat[1] = 0.0;
	m_camera.lookat[2] = 0.0;

	m_camera.worldup[0] = 0.0;
	m_camera.worldup[1] = 1.0;
	m_camera.worldup[2] = 0.0;

	m_camera.FOV = DEFAULT_FOV;

}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	for (int k = 0; k < AAKERNEL_SIZE + 1; k++) {
		delete[] pixelbuffer[k];
	}
	delete[] framebuffer;
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (int i = 0; i < xres * yres; i++) {
		for (int k = 0; k < AAKERNEL_SIZE + 1; k++) {
			pixelbuffer[k][i] = { 128 * 16, 112 * 16, 96 * 16, 0, INT_MAX };
		}
	}
	return GZ_SUCCESS;

}

int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/

	matlevel = 0;

	GzCoord z = { m_camera.lookat[0] - m_camera.position[0],  m_camera.lookat[1] - m_camera.position[1] , m_camera.lookat[2] - m_camera.position[2] };
	norm(z);
	float productZU = z[0] * m_camera.worldup[0] + z[1] * m_camera.worldup[1] + z[2] * m_camera.worldup[2];
	GzCoord c = { m_camera.worldup[0] - z[0] * productZU, m_camera.worldup[1] - z[1] * productZU, m_camera.worldup[2] - z[2] * productZU };
	norm(c);
	GzCoord x = { c[1] * z[2] - c[2] * z[1], c[2] * z[0] - c[0] * z[2], c[0] * z[1] - c[1] * z[0] };
	norm(x);

	GzMatrix Xiw = {
		{ x[0], x[1], x[2], -(x[0] * m_camera.position[0] + x[1] * m_camera.position[1] + x[2] * m_camera.position[2])},
		{ c[0], c[1], c[2], -(c[0] * m_camera.position[0] + c[1] * m_camera.position[1] + c[2] * m_camera.position[2])},
		{ z[0], z[1], z[2], -(z[0] * m_camera.position[0] + z[1] * m_camera.position[1] + z[2] * m_camera.position[2])},
		{ 0, 0, 0, 1}
	};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			m_camera.Xiw[i][j] = Xiw[i][j];
	}

	GzMatrix Xpi = {
		{ 1,0,0,0 },
		{ 0,1,0,0 },
		{ 0,0,tan(m_camera.FOV * PI / 180 / 2),0 },
		{ 0,0,tan(m_camera.FOV * PI / 180 / 2),1 }
	};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			m_camera.Xpi[i][j] = Xpi[i][j];
	}



	GzPushMatrix(Xsp, false);
	GzPushMatrix(Xpi, false);
	GzPushMatrix(Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/

	m_camera.position[0] = camera.position[0];
	m_camera.position[1] = camera.position[1];
	m_camera.position[2] = camera.position[2];

	m_camera.lookat[0] = camera.lookat[0];
	m_camera.lookat[1] = camera.lookat[1];
	m_camera.lookat[2] = camera.lookat[2];

	m_camera.worldup[0] = camera.worldup[0];
	m_camera.worldup[1] = camera.worldup[1];
	m_camera.worldup[2] = camera.worldup[2];

	m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;
}

int GzRender::GzPushNormMatrix(GzMatrix	matrix) {

	float S = 1 / sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2]);

	matrix[0][0] *= S;
	matrix[0][1] *= S;
	matrix[0][2] *= S;
	matrix[0][3] = 0;

	matrix[1][0] *= S;
	matrix[1][1] *= S;
	matrix[1][2] *= S;
	matrix[1][3] = 0;

	matrix[2][0] *= S;
	matrix[2][1] *= S;
	matrix[2][2] *= S;
	matrix[2][3] = 0;

	matrix[3][0] *= S;
	matrix[3][1] *= S;
	matrix[3][2] *= S;
	matrix[3][3] = 0;

	if (matNormLevel == 0)
		copyMatrix(Xnorm[0], matrix);

	else
		multiMatrix(Xnorm[matNormLevel], Xnorm[matNormLevel - 1], matrix);

	matNormLevel++;

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix, bool judge)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	if (matlevel == 0)
		copyMatrix(Ximage[0], matrix);

	else
		multiMatrix(Ximage[matlevel], Ximage[matlevel - 1], matrix);

	matlevel++;

	if (judge) {
		GzPushNormMatrix(matrix);
	}
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel > 0)
		matlevel--;
	if (matNormLevel > 0)
		matNormLevel--;
	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z, int k)
{
	/* HW1.4 write pixel values into the buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		if (z < pixelbuffer[k][j * xres + i].z)
			pixelbuffer[k][j * xres + i] = { min(4095, max(r, 0)), min(4095, max(g, 0)), min(4095, max(b, 0)), min(4095, max(a, 0)), z };
	}
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z, int k)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		*r = pixelbuffer[k][i * yres + j].red;
		*g = pixelbuffer[k][i * yres + j].green;
		*b = pixelbuffer[k][i * yres + j].blue;
		*a = pixelbuffer[k][i * yres + j].alpha;
		*z = pixelbuffer[k][i * yres + j].z;
	}

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\r", xres, yres);
	for (int i = 0; i < xres * yres; i++) {
		double temp0 = 0, temp1 = 0, temp2 = 0;
		for (int k = 0; k < AAKERNEL_SIZE; k++) {
			temp0 += (pixelbuffer[k][i].red * AAFilter[k][2]);
			temp1 += (pixelbuffer[k][i].green * AAFilter[k][2]);
			temp2 += (pixelbuffer[k][i].blue * AAFilter[k][2]);
		}
		pixelbuffer[AAKERNEL_SIZE][i].red = int(temp0);
		pixelbuffer[AAKERNEL_SIZE][i].green = int(temp1);
		pixelbuffer[AAKERNEL_SIZE][i].blue = int(temp2);
		char r = pixelbuffer[AAKERNEL_SIZE][i].red >> 4;
		char g = pixelbuffer[AAKERNEL_SIZE][i].green >> 4;
		char b = pixelbuffer[AAKERNEL_SIZE][i].blue >> 4;
		fwrite(&r, 1, 1, outfile);
		fwrite(&g, 1, 1, outfile);
		fwrite(&b, 1, 1, outfile);
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	for (int i = 0; i < xres * yres; i++) {
		framebuffer[i * 3 + 2] = (char)(pixelbuffer[AAKERNEL_SIZE][i].red >> 4);
		framebuffer[i * 3 + 1] = (char)(pixelbuffer[AAKERNEL_SIZE][i].green >> 4);
		framebuffer[i * 3] = (char)(pixelbuffer[AAKERNEL_SIZE][i].blue >> 4);
	}

	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int i = 0; i < numAttributes; i++) {
		if (nameList[i] == GZ_RGB_COLOR) {
			flatcolor[0] = ((GzColor*)valueList[i])[0][0];
			flatcolor[1] = ((GzColor*)valueList[i])[0][1];
			flatcolor[2] = ((GzColor*)valueList[i])[0][2];
		}
		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT) {
			lights[numlights] = *(GzLight*)valueList[i];
			numlights++;
		}
		else if (nameList[i] == GZ_AMBIENT_LIGHT) {
			ambientlight = *(GzLight*)valueList[i];
		}
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT) {
			Kd[0] = ((GzColor*)valueList[i])[0][0];
			Kd[1] = ((GzColor*)valueList[i])[0][1];
			Kd[2] = ((GzColor*)valueList[i])[0][2];
		}
		else if (nameList[i] == GZ_INTERPOLATE) {
			interp_mode = *(int*)valueList[i];
		}
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT) {
			Ka[0] = ((GzColor*)valueList[i])[0][0];
			Ka[1] = ((GzColor*)valueList[i])[0][1];
			Ka[2] = ((GzColor*)valueList[i])[0][2];
		}
		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT) {
			Ks[0] = ((GzColor*)valueList[i])[0][0];
			Ks[1] = ((GzColor*)valueList[i])[0][1];
			Ks[2] = ((GzColor*)valueList[i])[0][2];
		}
		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT) {
			spec = *(float*)valueList[i];
		}
		else if (nameList[i] == GZ_TEXTURE_MAP) {
			tex_fun = (GzTexture)valueList[i]; // what am i doing here?! I have no diea
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzShade(point& value, GzColor texK) {
	GzCoord nor;
	nor[0] = value.xn;
	nor[1] = value.yn;
	nor[2] = value.zn;
	GzCoord E = { 0, 0, -1 };
	//GzCoord* R = new GzCoord[numlights];
	GzColor specularLight = { 0, 0, 0 };
	GzColor diffuseLight = { 0, 0, 0 };
	GzColor ambientLight = { 0, 0, 0 };
	float judge = multiVector(E, nor);
	if (judge < 0) {
		nor[0] = -nor[0];
		nor[1] = -nor[1];
		nor[2] = -nor[2];
	}
	for (int i = 0; i < numlights; i++) {
		float dot = multiVector(lights[i].direction, nor);
		if (dot < 0)
			continue;
		GzCoord R;
		R[0] = 2 * dot*nor[0] - lights[i].direction[0];
		R[1] = 2 * dot*nor[1] - lights[i].direction[1];
		R[2] = 2 * dot*nor[2] - lights[i].direction[2];
		norm(R);

		float dotRE = multiVector(R, E);
		if (dotRE > 0) {
			specularLight[0] += lights[i].color[0] * pow(dotRE, spec);
			specularLight[1] += lights[i].color[1] * pow(dotRE, spec);
			specularLight[2] += lights[i].color[2] * pow(dotRE, spec);
		}

		diffuseLight[0] += lights[i].color[0] * dot;
		diffuseLight[1] += lights[i].color[1] * dot;
		diffuseLight[2] += lights[i].color[2] * dot;

	}
	specularLight[0] *= Ks[0];
	specularLight[1] *= Ks[1];
	specularLight[2] *= Ks[2];

	if (tex_fun) {
		diffuseLight[0] *= texK[0];
		diffuseLight[1] *= texK[1];
		diffuseLight[2] *= texK[2];

		ambientLight[0] = texK[0] * ambientlight.color[0];
		ambientLight[1] = texK[1] * ambientlight.color[1];
		ambientLight[2] = texK[2] * ambientlight.color[2];
	}
	else {
		diffuseLight[0] *= Kd[0];
		diffuseLight[1] *= Kd[1];
		diffuseLight[2] *= Kd[2];

		ambientLight[0] = Ka[0] * ambientlight.color[0];
		ambientLight[1] = Ka[1] * ambientlight.color[1];
		ambientLight[2] = Ka[2] * ambientlight.color[2];
	}

	value.red = specularLight[0] + diffuseLight[0] + ambientLight[0];
	value.green = specularLight[1] + diffuseLight[1] + ambientLight[1];
	value.blue = specularLight[2] + diffuseLight[2] + ambientLight[2];

	return GZ_SUCCESS;
}


int GzRender::GzTexShade(point& value) {
	GzCoord nor;
	nor[0] = value.xn;
	nor[1] = value.yn;
	nor[2] = value.zn;
	GzCoord E = { 0, 0, -1 };
	//GzCoord* R = new GzCoord[numlights];
	GzColor specularLight = { 0, 0, 0 };
	GzColor diffuseLight = { 0, 0, 0 };
	GzColor ambientLight = { 0, 0, 0 };
	float judge = multiVector(E, nor);
	if (judge < 0) {
		nor[0] = -nor[0];
		nor[1] = -nor[1];
		nor[2] = -nor[2];
	}
	for (int i = 0; i < numlights; i++) {
		float dot = multiVector(lights[i].direction, nor);
		if (dot < 0)
			continue;
		GzCoord R;
		R[0] = 2 * dot*nor[0] - lights[i].direction[0];
		R[1] = 2 * dot*nor[1] - lights[i].direction[1];
		R[2] = 2 * dot*nor[2] - lights[i].direction[2];
		norm(R);

		float dotRE = multiVector(R, E);
		if (dotRE > 0) {
			specularLight[0] += lights[i].color[0] * pow(dotRE, spec);
			specularLight[1] += lights[i].color[1] * pow(dotRE, spec);
			specularLight[2] += lights[i].color[2] * pow(dotRE, spec);
		}

		diffuseLight[0] += lights[i].color[0] * dot;
		diffuseLight[1] += lights[i].color[1] * dot;
		diffuseLight[2] += lights[i].color[2] * dot;

	}


	ambientLight[0] = ambientlight.color[0];
	ambientLight[1] = ambientlight.color[1];
	ambientLight[2] = ambientlight.color[2];

	value.red = specularLight[0] + diffuseLight[0] + ambientLight[0];
	value.green = specularLight[1] + diffuseLight[1] + ambientLight[1];
	value.blue = specularLight[2] + diffuseLight[2] + ambientLight[2];

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a ((GzCoord*)valueList[i])angle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	point value[3];
	for (int a = 0; a < AAKERNEL_SIZE; a++) {
		for (int i = 0; i < numParts; i++) {

			if (nameList[i] == GZ_POSITION) {

				float W1;
				for (int j = 0; j < 3; ++j) {
					value[j].x = (double)Ximage[matlevel - 1][0][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][0][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][0][2] * ((GzCoord*)valueList[i])[j][2]
						+ Ximage[matlevel - 1][0][3];
					value[j].y = (double)Ximage[matlevel - 1][1][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][1][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][1][2] * ((GzCoord*)valueList[i])[j][2]
						+ Ximage[matlevel - 1][1][3] * 1.0;
					value[j].z = (double)Ximage[matlevel - 1][2][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][2][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][2][2] * ((GzCoord*)valueList[i])[j][2]
						+ Ximage[matlevel - 1][2][3] * 1.0;
					W1 = (double)Ximage[matlevel - 1][3][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][3][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][3][2] * ((GzCoord*)valueList[i])[j][2]
						+ Ximage[matlevel - 1][3][3] * 1.0;
					value[j].x /= W1;
					value[j].y /= W1;
					value[j].z /= W1;

					value[j].x -= AAFilter[a][0];
					value[j].y -= AAFilter[a][1];
				}

				for (int j = 0; j < 3; j++) {
					if (value[j].z < 0)
						return GZ_SUCCESS;
				}
			}
			else if (nameList[i] == GZ_NORMAL) {
				for (int j = 0; j < 3; ++j) {
					value[j].xn = Xnorm[matNormLevel - 1][0][0] * ((GzCoord*)valueList[i])[j][0] + Xnorm[matNormLevel - 1][0][1] * ((GzCoord*)valueList[i])[j][1] + Xnorm[matNormLevel - 1][0][2] * ((GzCoord*)valueList[i])[j][2];
					value[j].yn = Xnorm[matNormLevel - 1][1][0] * ((GzCoord*)valueList[i])[j][0] + Xnorm[matNormLevel - 1][1][1] * ((GzCoord*)valueList[i])[j][1] + Xnorm[matNormLevel - 1][1][2] * ((GzCoord*)valueList[i])[j][2];
					value[j].zn = Xnorm[matNormLevel - 1][2][0] * ((GzCoord*)valueList[i])[j][0] + Xnorm[matNormLevel - 1][2][1] * ((GzCoord*)valueList[i])[j][1] + Xnorm[matNormLevel - 1][2][2] * ((GzCoord*)valueList[i])[j][2];
				}
			}
			if (nameList[i] == GZ_TEXTURE_INDEX) {
				for (int j = 0; j < 3; j++) {
					// transform uv into perspective space UV for each vertex
					float Vz = value[j].z / (INT_MAX - value[j].z);
					value[j].tex[0] = ((GzTextureIndex*)valueList[i])[j][0] / (Vz + 1);
					value[j].tex[1] = ((GzTextureIndex*)valueList[i])[j][1] / (Vz + 1);

				}
			}
		}

		GzTexShade(value[0]);
		GzTexShade(value[1]);
		GzTexShade(value[2]);

		sort_y(value);
		float A, B, C, D, AR, AG, AB, BR, BG, BB, DR, DG, DB, AX, AY, AZ, BX, BY, BZ, DX, DY, DZ, AU, AV, BU, BV, DU, DV;
		float x10, x20, y10, y20, z10, z20, r10, r20, g10, g20, b10, b20, xx10, xx20, yy10, yy20, zz10, zz20, u10, u20, v10, v20;
		x10 = (value[1].x - value[0].x);
		y10 = (value[1].y - value[0].y);
		z10 = (value[1].z - value[0].z);
		r10 = (value[1].red - value[0].red);
		g10 = (value[1].green - value[0].green);
		b10 = (value[1].blue - value[0].blue);
		xx10 = (value[1].xn - value[0].xn);
		yy10 = (value[1].yn - value[0].yn);
		zz10 = (value[1].zn - value[0].zn);
		u10 = (value[1].tex[0] - value[0].tex[0]);
		v10 = (value[1].tex[1] - value[0].tex[1]);


		x20 = (value[2].x - value[0].x);
		y20 = (value[2].y - value[0].y);
		z20 = (value[2].z - value[0].z);
		r20 = (value[2].red - value[0].red);
		g20 = (value[2].green - value[0].green);
		b20 = (value[2].blue - value[0].blue);
		xx20 = (value[2].xn - value[0].xn);
		yy20 = (value[2].yn - value[0].yn);
		zz20 = (value[2].zn - value[0].zn);
		u20 = (value[2].tex[0] - value[0].tex[0]);
		v20 = (value[2].tex[1] - value[0].tex[1]);


		A = ((y10 * z20) - (y20 * z10));
		AR = ((y10 * r20) - (y20 * r10));
		AG = ((y10 * g20) - (y20 * g10));
		AB = ((y10 * b20) - (y20 * b10));
		AX = ((y10 * xx20) - (y20 * xx10));
		AY = ((y10 * yy20) - (y20 * yy10));
		AZ = ((y10 * zz20) - (y20 * zz10));
		AU = ((y10 * u20) - (y20 * u10));
		AV = ((y10 * v20) - (y20 * v10));


		B = ((x20 * z10) - (x10 * z20));
		BR = ((x20 * r10) - (x10 * r20));
		BG = ((x20 * g10) - (x10 * g20));
		BB = ((x20 * b10) - (x10 * b20));
		BX = ((x20 * xx10) - (x10 * xx20));
		BY = ((x20 * yy10) - (x10 * yy20));
		BZ = ((x20 * zz10) - (x10 * zz20));
		BU = ((x20 * u10) - (x10 * u20));
		BV = ((x20 * v10) - (x10 * v20));


		C = ((x10 * y20) - (x20 * y10));

		D = ((value[0].y * (B)) + (value[0].x * (A)) + (value[0].z * (C))) * (-1);
		DR = ((value[0].y * (BR)) + (value[0].x * (AR)) + (value[0].red * (C))) * (-1);
		DG = ((value[0].y * (BG)) + (value[0].x * (AG)) + (value[0].green * (C))) * (-1);
		DB = ((value[0].y * (BB)) + (value[0].x * (AB)) + (value[0].blue * (C))) * (-1);
		DX = ((value[0].y * (BX)) + (value[0].x * (AX)) + (value[0].xn * (C))) * (-1);
		DY = ((value[0].y * (BY)) + (value[0].x * (AY)) + (value[0].yn * (C))) * (-1);
		DZ = ((value[0].y * (BZ)) + (value[0].x * (AZ)) + (value[0].zn * (C))) * (-1);
		DU = ((value[0].y * (BU)) + (value[0].x * (AU)) + (value[0].tex[0] * (C))) * (-1);
		DV = ((value[0].y * (BV)) + (value[0].x * (AV)) + (value[0].tex[1] * (C))) * (-1);
		

		float xt;
		if (value[0].x == value[2].x)
			xt = value[0].x;
		else
			xt = (value[0].x - value[2].x) * (value[1].y - value[2].y) / (value[0].y - value[2].y) + value[2].x;
		float x1 = min(xt, value[1].x);
		float x2 = max(xt, value[1].x);


		for (int i = ((int)(value[0].y)) + 1; i < ((int)(value[1].y)) + 1; i++) {
			int xl = (int)(((float)i - value[0].y) / (value[1].y - value[0].y) * (x1 - value[0].x) + value[0].x) + 1;
			int xr = (int)(((float)i - value[0].y) / (value[1].y - value[0].y) * (x2 - value[0].x) + value[0].x) + 1;
			for (int j = xl; j < xr; j++) {
				float z, Dz, Dzm;
				float dzdx, dzdy;
				z = -(A * j + B * i + D) / C;
				Dz = -(A * (j + 1) + B * (i + 1) + D) / C;
				Dzm = -(A * (j - 1) + B * (i - 1) + D) / C;
				dzdx = -(A * (j + 1) + B * i + D) / C;
				dzdy = -(A * j + B * (i + 1) + D) / C;

				float red, green, blue;
				GzTextureIndex UV;
				GzTextureIndex DUV, DUVm;
				UV[0] = -(AU * j + BU * i + DU) / C;
				UV[1] = -(AV * j + BV * i + DV) / C;
				DUV[0] = -(AU * (j + 1) + BU * (i + 1) + DU) / C;
				DUV[1] = -(AV * (j + 1) + BV * (i + 1) + DV) / C;
				DUVm[0] = -(AU * (j - 1) + BU * (i - 1) + DU) / C;
				DUVm[1] = -(AV * (j - 1) + BV * (i - 1) + DV) / C;
				float Vz = z / (INT_MAX - z);
				UV[0] = UV[0] * (Vz + 1);
				UV[1] = UV[1] * (Vz + 1);

				float DVz = Dz / (INT_MAX - Dz);
				DUV[0] = DUV[0] * (DVz + 1);
				DUV[1] = DUV[1] * (DVz + 1);

				float DVzm = Dzm / (INT_MAX - Dzm);
				DUVm[0] = DUVm[0] * (DVzm + 1);
				DUVm[1] = DUVm[1] * (DVzm + 1);

				int scaleSizeU = round((DUV[0] - UV[0]) * 50);
				int scaleSizeV = round((DUV[1] - UV[1]) * 50);
				int scaleSizeUm = round((DUVm[0] - UV[0]) * 50);
				int scaleSizeVm = round((DUVm[1] - UV[1]) * 50);
				int scaleSize = scaleSizeU * scaleSizeV;



				GzColor texK;
				if (tex_fun) {
					tex_fun(UV[0], UV[1], texK, scaleSizeU, scaleSizeV, scaleSizeUm, scaleSizeVm);
				}


				if (interp_mode == GZ_FLAT) {
					red = min(1.0, value[0].red);
					green = min(1.0, value[0].green);
					blue = min(1.0, value[0].blue);
				}
				else if (interp_mode == GZ_COLOR) {
					if (tex_fun) {
						red = min(1.0, -(AR * j + BR * i + DR) * texK[0] / C);
						green = min(1.0, -(AG * j + BG * i + DG) * texK[1] / C);
						blue = min(1.0, -(AB * j + BB * i + DB) * texK[2] / C);
					}
					else {
						red = min(1.0, -(AR * j + BR * i + DR) / C);
						green = min(1.0, -(AG * j + BG * i + DG) / C);
						blue = min(1.0, -(AB * j + BB * i + DB) / C);
					}
				}
				else if (interp_mode == GZ_NORMALS) {
					point p;
					p.x = j;
					p.y = i;
					p.z = z;
					p.xn = -(AX * j + BX * i + DX) / C;
					p.yn = -(AY * j + BY * i + DY) / C;
					p.zn = -(AZ * j + BZ * i + DZ) / C;
					GzShade(p, texK);
					red = min(1.0, p.red);
					green = min(1.0, p.green);
					blue = min(1.0, p.blue);
				}

				GzPut(j, i, ctoi(red), ctoi(green), ctoi(blue), 1, (int)z, a);
			}
		}

		for (int i = ((int)(value[1].y)) + 1; i < ((int)(value[2].y)) + 1; i++) {
			int xl = (int)(((float)i - value[1].y) / (value[2].y - value[1].y) * (value[2].x - x1) + x1) + 1;
			int xr = (int)(((float)i - value[1].y) / (value[2].y - value[1].y) * (value[2].x - x2) + x2) + 1;
			for (int j = xl; j < xr; j++) {
				float z, Dz, Dzm;
				float dzdx, dzdy;
				z = -(A * j + B * i + D) / C;
				Dz = -(A * (j + 1) + B * (i + 1) + D) / C;
				Dzm = -(A * (j - 1) + B * (i - 1) + D) / C;
				dzdx = -(A * (j + 1) + B * i + D) / C;
				dzdy = -(A * j + B * (i + 1) + D) / C;

				float red, green, blue;
				GzTextureIndex UV;
				GzTextureIndex DUV, DUVm;
				UV[0] = -(AU * j + BU * i + DU) / C;
				UV[1] = -(AV * j + BV * i + DV) / C;
				DUV[0] = -(AU * (j + 1) + BU * (i + 1) + DU) / C;
				DUV[1] = -(AV * (j + 1) + BV * (i + 1) + DV) / C;
				DUVm[0] = -(AU * (j - 1) + BU * (i - 1) + DU) / C;
				DUVm[1] = -(AV * (j - 1) + BV * (i - 1) + DV) / C;
				float Vz = z / (INT_MAX - z);
				UV[0] = UV[0] * (Vz + 1);
				UV[1] = UV[1] * (Vz + 1);

				float DVz = Dz / (INT_MAX - Dz);
				DUV[0] = DUV[0] * (DVz + 1);
				DUV[1] = DUV[1] * (DVz + 1);

				float DVzm = Dzm / (INT_MAX - Dzm);
				DUVm[0] = DUVm[0] * (DVzm + 1);
				DUVm[1] = DUVm[1] * (DVzm + 1);

				int scaleSizeU = round((DUV[0] - UV[0]) * 50);
				int scaleSizeV = round((DUV[1] - UV[1]) * 50);
				int scaleSizeUm = round((DUVm[0] - UV[0]) * 50);
				int scaleSizeVm = round((DUVm[1] - UV[1]) * 50);
				int scaleSize = scaleSizeU * scaleSizeV;



				GzColor texK;
				if (tex_fun) {
					tex_fun(UV[0], UV[1], texK, scaleSizeU, scaleSizeV, scaleSizeUm, scaleSizeVm);
				}


				if (interp_mode == GZ_FLAT) {
					red = min(1.0, value[0].red);
					green = min(1.0, value[0].green);
					blue = min(1.0, value[0].blue);
				}
				else if (interp_mode == GZ_COLOR) {
					if (tex_fun) {
						red = min(1.0, -(AR * j + BR * i + DR) * texK[0] / C);
						green = min(1.0, -(AG * j + BG * i + DG) * texK[1] / C);
						blue = min(1.0, -(AB * j + BB * i + DB) * texK[2] / C);
					}
					else {
						red = min(1.0, -(AR * j + BR * i + DR) / C);
						green = min(1.0, -(AG * j + BG * i + DG) / C);
						blue = min(1.0, -(AB * j + BB * i + DB) / C);
					}
				}
				else if (interp_mode == GZ_NORMALS) {
					point p;
					p.x = j;
					p.y = i;
					p.z = z;
					p.xn = -(AX * j + BX * i + DX) / C;
					p.yn = -(AY * j + BY * i + DY) / C;
					p.zn = -(AZ * j + BZ * i + DZ) / C;
					GzShade(p, texK);
					red = min(1.0, p.red);
					green = min(1.0, p.green);
					blue = min(1.0, p.blue);
				}

				GzPut(j, i, ctoi(red), ctoi(green), ctoi(blue), 1, (int)z, a);
			}
		}
	}

	return GZ_SUCCESS;
}

int cmp(point x, point y) {
	if (x.y != y.y)
		return x.y < y.y;
	return x.x < y.x;
}


void GzRender::sort_y(point* list) {
	std::sort(list, list + 3, cmp);
}

