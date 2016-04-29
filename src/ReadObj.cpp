/* 
	simple cpp source to read & render obj file

	authour: 'ling.shi'
	email: ling.shi@kaust.edu.sa; ling.shi@aliyun.com
	created date: 2011.4
	publish date: 2016.4
*/

#include <stdlib.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <gl/glui.h>
#include <math.h>
using namespace std;
#include "ReadBmp.h" 
#pragma comment (lib, "glui32.lib")
#pragma comment(linker,"/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#define EPS 0.000001

class MyObj
{
public:

	float fCenter[3]/*={0,0,0}*/; // the center position of all vertexss
	float m_fRGB[3][3]; //vertex, face, frames, RGB
	int m_bShowStatus[3];
	float m_fAng[3];
	float m_fScale;
	float m_fPos[3]; //x,y,z position
	int m_nVert;
	int m_nFace;//number
	float **m_pVertex;//positions
	float **m_pFace;

	float **m_pNormal;// normals to vertex
	char m_strName[255];
	MyObj(char* strFile)
	{
		strcpy(m_strName, strFile);
		m_fRGB[0][0]= 0.1f; m_fRGB[0][1] = 0.9f; m_fRGB[0][2] = 0.2f;
		m_fRGB[1][0]= 0.2f; m_fRGB[1][1] = 0.5f; m_fRGB[1][2] = 0.6f;
		m_fRGB[2][0]= 1.0f; m_fRGB[2][1] = 0.9f; m_fRGB[2][2] = 0.9f;
		m_bShowStatus[0]=m_bShowStatus[2]=true; m_bShowStatus[1]=0;
		m_fAng[0] = m_fAng[1] = m_fAng[2] = 0;
		m_fScale = 1.0;
		m_fPos[0] = m_fPos[1] = m_fPos[2] = 0;
		ReadObj(strFile);



		fCenter[0] = fCenter[1] = fCenter[2] = 0;
		for (int i=0; i<m_nVert; ++i)
		{
			for(int j=0;j<3;++j)
				fCenter[j] += m_pVertex[i][j];

		}
		fCenter[0]/=m_nVert;
		fCenter[1]/=m_nVert;
		fCenter[2]/=m_nVert;
		//m_fPos[0] -=fCenter[0];
		//m_fPos[1] -=fCenter[1];
		//m_fPos[2] -= fCenter[2];
		if (strcmp( m_strName, "cube.obj" )==0 )
		{ 
			 //m_fPos[0] = 0.7;
		}
		if (strcmp( m_strName,  "pyramid.obj" )==0 )
		{
			// m_fPos[0] = -0.7;
			 m_fAng[0] = -80;
		}
		if (strcmp( m_strName, "coati.obj" )==0 )
		{ 
			//m_fPos[0] = 0.7;
		}
	}

	~MyObj()
	{
		for (int i=0; i<m_nVert; ++i)
		{
			delete []m_pVertex[i];
		}
		delete []m_pVertex;
		m_pVertex = 0;

		for (int i=0; i<m_nFace; ++i)
		{
			delete []m_pFace;
		}
		delete []m_pFace;

		for (int i=0; i<m_nFace; ++i)
		{
			delete []m_pNormal[i];
		}
		delete []m_pNormal;
		m_pFace = 0;
		{
		}
	}
private:
	void ReadObj(char* strFile)
	{

		string s1;
		float s2,s3,s4;
		ifstream infile(strFile);
		string sline;
		vector<string> vec1;
		vector<float> vec2,vec3,vec4;

		int nNum=0;
		m_nVert = 0;
		m_nFace = 0;
		while(getline(infile, sline))
		{

			istringstream sin(sline);
			sin>>s1>>s2>>s3>>s4;
			if(s1 == "v" )  
			{
				vec1.push_back(s1);
				vec2.push_back(s2);
				vec3.push_back(s3);
				vec4.push_back(s4);
				s1= ""; 
				++ m_nVert;
			} 
			else if (s1 == "f")
			{
				vec1.push_back(s1);
				vec2.push_back(s2);
				vec3.push_back(s3);
				vec4.push_back(s4);
				s1= ""; 
				++ m_nFace;
			}
		}

		//allocate memory
		m_pVertex = new float*[m_nVert];
		for (int i=0; i<m_nVert; ++i)
		{
			m_pVertex[i] = new float[3];
		}
		m_pFace = new float*[m_nFace];
		for (int i=0; i<m_nFace; ++i)
		{
			m_pFace[i] = new float[3];
		}

		m_pNormal = new float *[m_nFace];
		for (int i=0; i<m_nFace; ++i)
		{
			m_pNormal[i] = new float[3];
		}

		//

		// give value
		int j=0, k=0;
		for (int i=0; i<m_nFace+m_nVert; ++i)
		{
			if (vec1[i]=="v")
			{
				m_pVertex[j][0]=vec2[i];
				m_pVertex[j][1]=vec3[i];
				m_pVertex[j][2]=vec4[i];
				++j;
			}
			else if(vec1[i]=="f")
			{
				m_pFace[k][0]=vec2[i];
				m_pFace[k][1]=vec3[i];
				m_pFace[k][2]=vec4[i];
				++k;
			}
		} 

		for (int i=0;i<m_nFace;++i)
		{
			float v[3][3];
			for(int j=0; j<3; ++j)
			{
				int idx=m_pFace[i][j]-1;
				v[j][0]=m_pVertex[idx][0];
				v[j][1]=m_pVertex[idx][1];
				v[j][2]=m_pVertex[idx][2];
			}
			float v2[3]={v[2][0]-v[0][0], v[2][1]-v[0][1], v[2][2]-v[0][2]  };
			float v1[3]={v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]  };
			m_pNormal[i][0] = v1[1]*v2[2]-v1[2]*v2[1];
			m_pNormal[i][1] = v1[2]*v2[0]-v1[0]*v2[2];
			m_pNormal[i][2] = v1[0]*v2[1]-v1[1]*v2[0];
			float N=sqrt(m_pNormal[i][0]*m_pNormal[i][0]+m_pNormal[i][1]*m_pNormal[i][1]+m_pNormal[i][2]*m_pNormal[i][2]);
			
			 
			if (abs(N)<EPS)
			{
				m_pNormal[i][0]=m_pNormal[i][1]=0; m_pNormal[i][2]=1;
			}
			else
			for(int j=0; j<3;++j) m_pNormal[i][j] /=N;
		}



	}


};

MyObj *pObj = NULL;
MyObj *pObjArr[2];
float gfRGB[3];
int gBShow[3];
float translate_xy[2] ;
float rotation_matrix[16]	//  Rotation Matrix Live Variable Array
= { 1.0, 0.0, 0.0, 0.0, 
0.0, 1.0, 0.0, 0.0,
0.0, 0.0, 1.0, 0.0, 
0.0, 0.0, 0.0, 1.0 };	

int color=0; //edge, frame, face
int color2=1;

int objidx=1;
bool mouseLeftDown;
bool mouseRightDown;
float mouseX, mouseY;
float cameraAngleX = 0;
float cameraAngleY = 0;
float cameraDistance;

float AngX=0.0f;
float AngY=0.0f;
float AngZ=0.0f;

float AngVx = 0.0f;
float AngVy = 0.0f;
float AngVz = 0.0f;

float pVertex[50][3];
float pface[50][3];
int nVer=0, nFace =0 ;
int gbTexture = 0;
int gbLight = 0;

float gfLightPos[] = {1.0, 1.0, 2.0, 0.0};
float gfLight1Pos[] = {-1.0, 1.0, 2.0, 0.0};
float gfAmbRGB[] = {1.0, 1.0, 1.0, 0.0};
float gfDiffRGB[] = {1.0, 1.0, 1.0, 0.0};
float gfSpecRGB[] = {1.0, 1.0, 1.0, 0.0};
float *gfPRGB = gfAmbRGB;
int gnColor=0;
//  define the window position on screen
int window_x;
int window_y;
GLuint main_window;//  pointer to the GLUI window
GLUI * glui_window;
GLUI *glui_lighting;
void setupGLUI ();
GLUI_Listbox *color_listbox;
GLUI_Spinner *spinner[3];
GLUI_Panel *op_panel;
GLUI_Checkbox  *MyCheckColor[3];
int gBMateral[4] = { 1, 1, 1, 1};
float gBLight0[4] = {1, 1,1,0};
float gBLight1[4] = {1, 0,0,0};
int gBLightOpen1= 0;
float gLightClr[3][4] = 	{
	1.0 ,0.0,0.0, 0.0,
	0.0, 1.0,0.0, 0.0,
	0.0,0.0, 1.0, 0.0};

	int gbTeapot = 0;


int draw = 1;				 

#define    checkImageWidth 256
#define    checkImageHeight 256
GLubyte checkImage[checkImageWidth][checkImageHeight][3];
void makeCheckImage(void)
{
	int i, j,   c;

	unsigned char *pImage = new unsigned char[256*256*3];
	bmp_read(pImage, 256,256, "head.bmp");

	for (i = 0; i < checkImageWidth; i++) {
		for (j = 0; j < checkImageHeight; j++) {
			c = (i+j*checkImageWidth)*3;
			checkImage[i][j][0] = pImage[c];
			checkImage[i][j][1] = pImage[c+1];
			checkImage[i][j][2] = pImage[c+2];
		}
	}
}

int idx=0;
int idx1=0;
int idx2=0;
int idx3=0;
// Initiliaze the OpenGL window
void init(int control_id=0)
{  


	glLineWidth(0.02);
	if(0==idx)
	{
 
		for(int i=0;i<3;++i)
		{
			gfRGB[i]=pObj->m_fRGB[color][i];
			gBShow[i]=pObj->m_bShowStatus[i];
		}
		idx++;

	}
	else
	{ 
		for(int i=0;i<3;++i)
		{pObj->m_fRGB[color][i]=gfRGB[i];
		pObj->m_bShowStatus[i]=gBShow[i];
		}
 
	}

	switch (control_id)
	{
	case 21://lighting:
		if(gbLight)
			glui_lighting->enable();
		else
			glui_lighting->disable();
		break;
	case 98:
		pObj->m_fScale = 1;
		break;
	case 99:


		pObj->m_fPos[0]=-pObj->fCenter[0];
		pObj->m_fPos[1] =-pObj->fCenter[1];
		pObj->m_fPos[2] =- pObj->fCenter[2];
		//pObj->m_fPos[0]=pObj->m_fPos[1]=pObj->m_fPos[2]=0 ;
		break;

	case 2:
		color = color2 -1;
		 
	case 1:
		pObj = pObjArr[objidx-1];  
		for(int i=0;i<3;++i) 
		{
			gfRGB[i]=pObj->m_fRGB[color][i];
			gBShow[i]=pObj->m_bShowStatus[i];
			MyCheckColor[i]->set_int_val(gBShow[i]);

		} 
			idx++;
	 
		for(int i=0;i<3;++i) 
		{
			spinner[i]->set_float_val(pObj->m_fRGB[color][i]);

		}

		break;  
	}

 


	glClearColor (0.0, 0.0, 0.0, 0.0);			// Clear the color 
	glShadeModel (GL_FLAT);						// Set the shading model to GL_FLAT
	glEnable (GL_LINE_SMOOTH);
	glEnable(GL_NORMALIZE);
	//glEnable(GL_RESCALE_NORMAL);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);		// Set Line Antialiasing

	//light
	GLfloat mat_specular[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess[] = { 150.0 };
	//GLfloat light_position[] = { -1.0, 1.0, 1.0, 0.0 };
	GLfloat light_diffuse[] = { 1.0,1.0, 1.0, 1.0 };
	//GLfloat light_ambient[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat mat_zero[] = {0.0, 0.0, 0.0, 1.0 };
  
	if(gBMateral[0])
	{
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
 
	}
	else
	{
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_zero);
 	}

	if(gBMateral[1])
	{
		//glLightfv(GL_LIGHT0, GL_AMBIENT, light_diffuse);
 		glMaterialfv(GL_FRONT, GL_AMBIENT, gLightClr[1]);

	}
	else
	{
		//glLightfv(GL_LIGHT0, GL_AMBIENT, mat_zero);
 		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_zero);
	}
	if(gBMateral[2])
	{
		//glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, gLightClr[0]);

	}
	else
	{
		//glLightfv(GL_LIGHT0, GL_DIFFUSE, mat_zero);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_zero);
	}
	if(gBMateral[3])
	{
		//glLightfv(GL_LIGHT0, GL_SPECULAR, light_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, gLightClr[2]);

	}
	else
	{
		//glLightfv(GL_LIGHT0, GL_SPECULAR, mat_zero);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_zero);
	}
	//glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);

	//glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);

	//glLightfv(GL_LIGHT0, GL_AMBIENT, gfAmbRGB);    
	glLightfv(GL_LIGHT0, GL_AMBIENT, gBLight0);  
	glLightfv(GL_LIGHT0, GL_DIFFUSE, gBLight0); 
	glLightfv(GL_LIGHT0, GL_SPECULAR, gBLight0);  
	glLightfv(GL_LIGHT0, GL_POSITION, gfLightPos);
      
	if(gbLight)
	{ 
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		if(gBLightOpen1)
		{ 
			glLightfv(GL_LIGHT1, GL_AMBIENT, gBLight1);  
			glLightfv(GL_LIGHT1, GL_DIFFUSE, gBLight1); 
			glLightfv(GL_LIGHT1, GL_SPECULAR, gBLight1);  
			glLightfv(GL_LIGHT1, GL_POSITION, gfLight1Pos);
			glEnable(GL_LIGHT1);
		}
		else
			glDisable(GL_LIGHT1); 
	}
	else
	{
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0); 
	}
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);


	//Texture

	makeCheckImage();
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, checkImageWidth, 
		checkImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, 
		&checkImage[0][0][0]);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
		GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
		GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	if(gbTexture)
		glEnable(GL_TEXTURE_2D);
	else
		glDisable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);

}

void displayobj(MyObj* pObj)
{

	glPushMatrix();

	glScalef(pObj->m_fScale, pObj->m_fScale,pObj->m_fScale); 
	glTranslatef(pObj->m_fPos[0], pObj->m_fPos[1], pObj->m_fPos[2]);	

	glRotatef(pObj->m_fAng[0], 1.0, 0.0, 0.0);
	glRotatef(pObj->m_fAng[1], 0.0, 1.0,  0.0);
	glRotatef(pObj->m_fAng[2], 0.0, 0.0, 1.0);
	glMultMatrixf (rotation_matrix);
	if(pObj->m_bShowStatus[0])
	{
		glColor3fv(pObj->m_fRGB[0]); 
		//glColor3f(0.1, 0.9, 0.2);
		glPointSize(2);
		glBegin(GL_POINTS);
		for (int i=0; i<pObj->m_nVert; ++i)
		{
			//glVertex3f(pVertex[0][i],pVertex[1][i], pVertex[2][i]);
			glVertex3fv(pObj->m_pVertex[i]);
		}

		glEnd();

	}

	if(pObj->m_bShowStatus[1])
	{
		glColor3fv(pObj->m_fRGB[1]);
		//glColor3f(0.2, 0.5, 0.6);
		glBegin(GL_TRIANGLES);

		for (int i=0; i<pObj->m_nFace; ++i) 
		{
			//glNormal3f (1, 0, 0);
			glNormal3fv(pObj->m_pNormal[i]);
			for(int j = 0; j<3; ++j)
			{
				int idx=pObj->m_pFace[i][j]-1; 
				glVertex3fv(pObj->m_pVertex[idx]);
			}

		}

		glEnd();

	}
	if(pObj->m_bShowStatus[2])
	{
		glColor3fv(pObj->m_fRGB[2]);
		//glColor3f(1.0, 0.9, 0.9);
		glLineWidth(0.02);

		for (int i=0; i<pObj->m_nFace; ++i) 
		{
			glBegin(GL_LINE_LOOP);
			for(int j = 0; j<3; ++j)
			{
				int idx=pObj->m_pFace[i][j]-1; 
				glVertex3fv(pObj->m_pVertex[idx]);
			}
			glEnd();
		}

	}

	if (gbTexture)
	{
		if (strcmp(pObj->m_strName, "cube.obj" )==0 )
		{ 
			glBegin(GL_QUADS);
			glTexCoord2f(0.0, 0.0); glVertex3fv(pObj->m_pVertex[0]);
			glTexCoord2f(0.0, 1.0); glVertex3fv(pObj->m_pVertex[1]);
			glTexCoord2f(1.0, 1.0); glVertex3fv(pObj->m_pVertex[3]);
			glTexCoord2f(1.0, 0.0); glVertex3fv(pObj->m_pVertex[2]);  
 

			glTexCoord2f(0.0, 0.0); glVertex3fv(pObj->m_pVertex[0]);
			glTexCoord2f(0.0, 1.0); glVertex3fv(pObj->m_pVertex[1]);
			glTexCoord2f(1.0, 1.0); glVertex3fv(pObj->m_pVertex[5]);
			glTexCoord2f(1.0, 0.0); glVertex3fv(pObj->m_pVertex[4]); 

			glTexCoord2f(0.0, 0.0); glVertex3fv(pObj->m_pVertex[0]);
			glTexCoord2f(0.0, 1.0); glVertex3fv(pObj->m_pVertex[2]);
			glTexCoord2f(1.0, 1.0); glVertex3fv(pObj->m_pVertex[6]);
			glTexCoord2f(1.0, 0.0); glVertex3fv(pObj->m_pVertex[4]);  
			glEnd();
		}

		if (strcmp(pObj->m_strName,  "pyramid.obj" )==0 )
		{
 
			glBegin(GL_QUADS);
			glTexCoord2f(0.0, 0.0); glVertex3fv(pObj->m_pVertex[0]);
			glTexCoord2f(0.0, 1.0); glVertex3fv(pObj->m_pVertex[1]);
			glTexCoord2f(1.0, 1.0); glVertex3fv(pObj->m_pVertex[2]);
			glTexCoord2f(1.0, 0.0); glVertex3fv(pObj->m_pVertex[3]);  

			glTexCoord2f(0.0, 0.0); glVertex3fv(pObj->m_pVertex[0]);
			glTexCoord2f(0.0, 1.0); glVertex3fv(pObj->m_pVertex[1]);
			glTexCoord2f(1.0, 1.0); glVertex3fv(pObj->m_pVertex[4]);
			glTexCoord2f(1.0, 0.0); glVertex3fv(pObj->m_pVertex[4]);  

 
			glEnd();      

		}
	}

	glPopMatrix();
}

void display(void)
{	 
	//init();
	//it may be dangerous...
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();								 

	glTranslatef(0, 0, cameraDistance);
	glRotatef(cameraAngleX, 1, 0, 0);   // pitch
	glRotatef(cameraAngleY, 0, 1, 0);   // heading
	glRotatef(AngX,1.0,0.0,0.0);				// Rotate on x
	glRotatef(AngY,0.0,1.0,0.0);				// Rotate on y
	glRotatef(AngZ,0.0,0.0,1.0);				// Rotate on z
	glBegin(GL_LINES);
	glColor3f (0.0, 1.0, 0.0);			// Green for x axis
	glVertex3f(0,0,0);
	glVertex3f(10,0,0);
	glColor3f(1.0,0.0,0.0);				// Red for y axis
	glVertex3f(0,0,0);
	glVertex3f(0,10,0);
	glColor3f(0.0,0.0,1.0);				// Blue for z axis
	glVertex3f(0,0,0);	
	glVertex3f(0,0,10);
	glEnd();

	// Dotted lines for the negative sides of x,y,z
	glEnable(GL_LINE_STIPPLE);				// Enable line stipple to use a dotted pattern for the lines
	glLineStipple(1, 0x0101);				// Dotted stipple pattern for the lines
	glBegin(GL_LINES);				
	glColor3f (0.0, 1.0, 0.0);			// Green for x axis
	glVertex3f(-10,0,0);
	glVertex3f(0,0,0);
	glColor3f(1.0,0.0,0.0);				// Red for y axis
	glVertex3f(0,0,0);
	glVertex3f(0,-10,0);
	glColor3f(0.0,0.0,1.0);				// Blue for z axis
	glVertex3f(0,0,0);
	glVertex3f(0,0,-10);
	glEnd();

	glDisable(GL_LINE_STIPPLE);				// Disable the line stipple



	glScalef(4.0,4.0,4.0); 
	 
	if(gbTeapot)
	{  
		glColor3fv(pObj->m_fRGB[0]);
		//auxSolidSphere(1.0);
		glColor3f(0.0, 1.0, 1.0);
		glutSolidTeapot(1.2);
	}

	GLuint index = glGenLists(1);

	// compile the display list, store a triangle in it
	glNewList(index, GL_COMPILE); 

	displayobj(pObjArr[1]); 
	glEndList(); 
	//displayobj(pObjArr[1]); 

 
	displayobj(pObjArr[0]); 
	glTranslatef(-2.0, 0.5f, 0);
	glCallList(index);
	 glFlush();
	//auxSolidSphere(1.0);
	//glutSolidSphere(1.2, 20, 20);
	//glFlush(); 


	glEnd();
	glPopMatrix();						 
	glutSwapBuffers();
}


void reshape (int w, int h)
{
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);				// Set the viewport
	glMatrixMode (GL_PROJECTION);								// Set the Matrix mode
	glLoadIdentity ();											
	gluPerspective(75, (GLfloat) w /(GLfloat) h , 0.10, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt (AngVx, AngVy, 25.0 + AngVz, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); 
}

void mouse(int button, int state, int x, int y) {

	mouseX = x;
	mouseY = y;

	if(button == GLUT_LEFT_BUTTON)
	{
		if(state == GLUT_DOWN)
		{
			mouseLeftDown = true;
		}
		else if(state == GLUT_UP)
			mouseLeftDown = false;
	}

	else if(button == GLUT_RIGHT_BUTTON)
	{
		if(state == GLUT_DOWN)
		{
			mouseRightDown = true;
		}
		else if(state == GLUT_UP)
			mouseRightDown = false;
	}
}
void mouseMotionCB(int x, int y)
{
	if(mouseLeftDown)
	{
		cameraAngleY += (x - mouseX);
		cameraAngleX += (y - mouseY);
		mouseX = x;
		mouseY = y;
	}
	if(mouseRightDown)
	{
		cameraDistance += (y - mouseY) * 0.2f;
		mouseY = y;
	}

	glutPostRedisplay();
}


void keyboard (unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'n':
		AngVz = -50;
		break;
	case 'N':
		AngVz =  0;
		break;

	case 'l':
		gbLight = !gbLight;
		init();

		
		break;
	case 't':
		gbTexture = !gbTexture;
		init();
		break;
	case 'r':
		if(gbLight)
		{
			gfPRGB[0]+=0.1;
			if(gfPRGB[0]>1.0)
				gfPRGB[0] = 1.0;
			init();
		}
		else
		{
			pObj->m_fRGB[color][0] +=0.1;
			if(pObj->m_fRGB[color][0]>1.0)
				pObj->m_fRGB[color][0] = 1.0;
		}
		break;
	case 'R':
		if(gbLight)
		{
			gfPRGB[0]-=0.1;
			if(gfPRGB[0]<0.0)
				gfPRGB[0] = 0.0;
			init();
		}
		else
		{
		pObj->m_fRGB[color][0] -=0.1;
		if(pObj->m_fRGB[color][0]<0.0)
			pObj->m_fRGB[color][0] = 0.0;
		}
		break;
	case 'g':
		if(gbLight)
		{
			gfPRGB[1]+=0.1;
			if(gfPRGB[1]>1.0)
				gfPRGB[1] = 1.0;
			init();
		}
		else
		{
		pObj->m_fRGB[color][1] +=0.1;
		if(pObj->m_fRGB[color][1]>1.0)
			pObj->m_fRGB[color][1] = 1.0;
		}
		break;
	case 'G':
		if(gbLight)
		{
			gfPRGB[1]-=0.1;
			if(gfPRGB[1]<0.0)
				gfPRGB[1] = 0.0;
			init();
		}
		else
		{
		pObj->m_fRGB[color][1] -=0.1;
		if(pObj->m_fRGB[color][1]<0.0)
			pObj->m_fRGB[color][1] = 0.0;
		}
		break;
	case 'b':
		if(gbLight)
		{
			gfPRGB[2]+=0.1;
			if(gfPRGB[2]>1.0)
				gfPRGB[2] = 1.0;
			init();
		}
		else
		{
		pObj->m_fRGB[color][2] +=0.1;
		if(pObj->m_fRGB[color][2]>1.0)
			pObj->m_fRGB[color][2] = 1.0;
		}
		break;
	case 'B':
		if(gbLight)
		{
			gfPRGB[2]-=0.1;
			if(gfPRGB[2]<0.0)
				gfPRGB[2] = 0.0;
			init();
		}
		else
		{
		pObj->m_fRGB[color][2] -=0.1;
		if(pObj->m_fRGB[color][2]<0.0)
			pObj->m_fRGB[color][2] = 0.0;
		}
		break;
	case '1': 
		pObj = pObjArr[0];
		break;
	case '2': 
		pObj = pObjArr[1];
		break;
	case '3':  
		gnColor = 0;
		gfPRGB = gfAmbRGB;
		init();
		break;
	case '4':  
		gnColor = 1;
		gfPRGB = gfDiffRGB;
		init();
		break;
	case '5':  
		gnColor = 2;
		gfPRGB = gfSpecRGB;
		init();
		break; 
	case 'q':
		pObj->m_fAng[0] += 2.0;
		if (pObj->m_fAng[0]>360.0)
		{
			pObj->m_fAng[0] = 0.0f;
		}
		break;
	case 'Q':
		pObj->m_fAng[0] -= 2.0;
		if (pObj->m_fAng[0]<0.0 )
		{
			pObj->m_fAng[0] = 358.0f;
		}
		break;
	case 'w':
		pObj->m_fAng[1] += 2.0;
		if (pObj->m_fAng[1]>360.0)
		{
			pObj->m_fAng[1] = 0.0f;
		}
		break;
	case 'W':
		pObj->m_fAng[1] -= 2.0;
		if (pObj->m_fAng[1]<0.0 )
		{
			pObj->m_fAng[1] = 358.0f;
		}
		break;
	case 'e':
		pObj->m_fAng[2] += 2.0;
		if (pObj->m_fAng[2]>360.0)
		{
			pObj->m_fAng[2] = 0.0f;
		}
		break;
	case 'E':
		pObj->m_fAng[2] -= 2.0;
		if (pObj->m_fAng[2]<0.0 )
		{
			pObj->m_fAng[2] = 358.0f;
		}
		break;
	case '-':
		pObj->m_fScale *= 0.48;
		break;
	case '=':
		pObj->m_fScale /= 0.48f;
		break;
	case 'O':
		pObj->m_fScale = 1.0f;
		break;
	case 'o':
	//Obj->m_fPos[0]=pObj->m_fPos[1]=pObj->m_fPos[2]=0 ;
		pObj->m_fPos[0]=-pObj->fCenter[0];
		pObj->m_fPos[1] =-pObj->fCenter[1];
		pObj->m_fPos[2] =- pObj->fCenter[2];

		break;

	case 'i':
		pObj->m_fPos[0] -= 0.05f;
		break;
	case 'I':
		pObj->m_fPos[0] += 0.05f;
		break;
	case 'j':
		pObj->m_fPos[1] -= 0.05f;
		break;
	case 'J':
		pObj->m_fPos[1] += 0.05f;
		break;
	case 'k':
		pObj->m_fPos[2] -= 0.05f;
		break;
	case 'K':
		pObj->m_fPos[2] += 0.05f;
		break;

	case '6':
		pObj->m_bShowStatus[0] = !pObj->m_bShowStatus[0];
		color = 0;
		break;
	case '7':
		pObj->m_bShowStatus[1] = !pObj->m_bShowStatus[1];
		color = 1;
		break;
	case '8':
		color = 2;
		pObj->m_bShowStatus[2] = !pObj->m_bShowStatus[2];
		break;

	case '9':
		AngVz -= 1.0f;
		break;
	case '0':
		AngVz += 1.0f;
		break;	
	case 's':
		if (gbLight)
		{
			gfLightPos[0] -=0.05;
			init();
		} 
		break;
	case 'x':		
		AngX -= 0.5f;
		break;
	case 'S':	
		if (gbLight)
		{
			gfLightPos[0] +=0.05;
			init();
		} 	 
		break;
	case 'X':	
		AngX += 0.5f;
		break;
	case 'd':
		if (gbLight)
		{
			gfLightPos[1] -=0.05;
			init();
		} 	 
		break;
	case 'y':		
		AngY -= 0.5f;
		break;
	case 'D':	
		if (gbLight)
		{
			gfLightPos[1] +=0.05;
			init();
		} 	 
		break;
	case 'Y':	
		AngY += 0.5f;	
		break;	
	case 'f':	
		if (gbLight)
		{
			gfLightPos[2] -=0.05;
			init();
		} 	 
		break;
	case 'z':	
		AngZ -= 0.5f;
		break;
	case 'F':	
		if (gbLight)
		{
			gfLightPos[2] +=0.05;
			init();
		} 	 
		break;
	case 'Z':	
		AngZ += 0.5f;
		break;

	}
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(AngVx, AngVy, 25.0 + AngVz, 0,0,0,0,1,0);
	glutPostRedisplay();
} 

void idle ()
{
	glutSetWindow (main_window);
	glutPostRedisplay ();
	//sleep (50);
}
void setupGLUI ()
{
	//  Set idle function
	GLUI_Master.set_glutIdleFunc (idle);

	//  Create GLUI window
	glui_window = GLUI_Master.create_glui ("Options", 0, window_x - 320, window_y);

	 
	//---------------------------------------------------------------------
	// 'Object Properties' Panel
	//---------------------------------------------------------------------

	//  Add the 'Object Properties' Panel to the GLUI window
	op_panel = glui_window->add_panel ("Object Properties");

	//  Add the Draw Check box to the 'Object Properties' Panel
	glui_window->add_checkbox_to_panel (op_panel, "Lighting", &gbLight, 21, init  );
	glui_window->add_column_to_panel (op_panel, false);
	glui_window->add_checkbox_to_panel (op_panel, "Texture", &gbTexture , -1, init  );//gbTeapot
	glui_window->add_checkbox_to_panel (op_panel, "Teapot", &gbTeapot , -1, init  );


	GLUI_Panel *op_panel2 = glui_window->add_panel ("Character");
	GLUI_Listbox *color_listbox = glui_window->add_listbox_to_panel (op_panel2, 
		"Object", &objidx, 1, init);

	//  Add the items to the listbox
	/*
	color_listbox->add_item (1, "Cube");
	color_listbox->add_item (2, "Pyramid");  */

	color_listbox->add_item (1, pObjArr[0]->m_strName);
	color_listbox->add_item (2, pObjArr[1]->m_strName);

	 gBShow[0] = pObj->m_bShowStatus[0];
	 gBShow[1] = pObj->m_bShowStatus[1];
	 gBShow[2] = pObj->m_bShowStatus[2];
	 MyCheckColor[0]=glui_window->add_checkbox_to_panel (op_panel2, "Vertex", gBShow, 6, init  );
	MyCheckColor[1]= glui_window->add_checkbox_to_panel (op_panel2, "Face", gBShow+1 , 7, init  ); 
	 MyCheckColor[2]=glui_window->add_checkbox_to_panel (op_panel2, "Frame", gBShow+2 , 8, init  ); 
 
	GLUI_Listbox *color_listbox2 = glui_window->add_listbox_to_panel (op_panel2, 
		"", &color2, 2, init); 

	color_listbox2->add_item (1, "Vertex");
	color_listbox2->add_item (2, "Face"); 
	color_listbox2->add_item (3, "Frame"); 
 
	gfRGB[0] = pObj->m_fRGB[color][0];
	gfRGB[1] = pObj->m_fRGB[color][1];
	gfRGB[2] = pObj->m_fRGB[color][2];
	spinner[0]=glui_window->add_spinner_to_panel(
		op_panel2, "Red",  GLUI_SPINNER_FLOAT, gfRGB,5,init);
	spinner[0]->set_float_limits(0.0, 1.0);
	 

	spinner[1]=glui_window->add_spinner_to_panel(
		op_panel2, "Green",  GLUI_SPINNER_FLOAT, gfRGB+1,3,init);
	spinner[1]->set_float_limits(0.0, 1.0);

	spinner[2]=glui_window->add_spinner_to_panel(
		op_panel2, "Blue",  GLUI_SPINNER_FLOAT, gfRGB+2,4,init);
	spinner[2]->set_float_limits(0.0, 1.0);

	//transform
 	GLUI_Panel *transformation_panel = glui_window->add_panel ("Transformation XY");

 	GLUI_Panel *transformation_panel1 = glui_window->add_panel_to_panel (transformation_panel, "");
 	GLUI_Translation *translation_xy = glui_window->add_translation_to_panel (
		transformation_panel1, "Translation Cube", GLUI_TRANSLATION_XY, pObjArr[0]->m_fPos, 10, init );

  	translation_xy->set_speed( 0.005 );	
 	glui_window->add_column_to_panel (transformation_panel1, false);
 	GLUI_Translation *translation_xy2 = glui_window->add_translation_to_panel (
		transformation_panel1, "Translation Pyramid", GLUI_TRANSLATION_XY, pObjArr[1]->m_fPos, 11, init );
	translation_xy2->set_speed( 0.005 ); 


	GLUI_Panel *transformation_panel2 = glui_window->add_panel ("Transformation Z");
  
	GLUI_Translation *translation_xy3 = glui_window->add_translation_to_panel (
		transformation_panel2, "Translation Cube", GLUI_TRANSLATION_Z, pObjArr[0]->m_fPos+2, 12, init );
 	GLUI_Panel *transformation_panel3 = glui_window->add_panel_to_panel (transformation_panel2, "");

	//  Add the rotation control
	glui_window->add_rotation_to_panel (transformation_panel3, "Rotation", rotation_matrix, 14 );
	glui_window->add_column_to_panel (transformation_panel2, false);
	GLUI_Translation *translation_xy4 = glui_window->add_translation_to_panel (
		transformation_panel2, "Translation Cube", GLUI_TRANSLATION_Z, pObjArr[1]->m_fPos+2, 13, init );
	translation_xy3->set_speed( 0.005 );	
	translation_xy4->set_speed( 0.005 );	



	//  Add separator
	glui_window->add_separator_to_panel (transformation_panel3);


	GLUI_Spinner *Sp1=glui_window->add_spinner_to_panel(
		transformation_panel2, "Scale Cube",  GLUI_SPINNER_FLOAT, &pObjArr[0]->m_fScale);
	Sp1->set_float_limits(0.0, 10.0);
	GLUI_Spinner *Sp2=glui_window->add_spinner_to_panel(
		transformation_panel2, "Scale Pyra",  GLUI_SPINNER_FLOAT, &pObjArr[1]->m_fScale);
	Sp2->set_float_limits(0.0, 10.0);

	glui_window->add_button_to_panel (transformation_panel2,"Unit", 98, init);
	glui_window->add_button_to_panel (transformation_panel2,"Origin", 99, init);

	//  Let the GLUI window know where its main graphics window is
	glui_window->set_main_gfx_window( main_window );



	glui_lighting= GLUI_Master.create_glui ("Light", 0, window_x+525, window_y); 
	GLUI_Panel *op_pl1 = glui_lighting->add_panel ("Lighting"); 

	glui_lighting->add_checkbox_to_panel (op_pl1, "Shining", &gBMateral[0], 101, init  );
	glui_lighting->add_checkbox_to_panel (op_pl1, "Diffusion", &gBMateral[2], 103, init  );
	 glui_lighting->add_column_to_panel (op_pl1, false);
	glui_lighting->add_checkbox_to_panel (op_pl1, "Ambient", &gBMateral[1] , 102, init  );
	glui_lighting->add_checkbox_to_panel (op_pl1, "Specular", &gBMateral[3], 104, init  );



	GLUI_Panel *op_pl3[3];
	GLUI_Spinner *pSpin;
	char *strNames[3]={"diffusion", "Ambient", "Specular"};

	for(int i=0; i<3; i++)
	{
		op_pl3[i] = glui_lighting->add_panel (strNames[i]); 
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl3[i], "Red",  GLUI_SPINNER_FLOAT, &gLightClr[i][0],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl3[i], "Green",  GLUI_SPINNER_FLOAT, &gLightClr[i][1],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
	glui_lighting->add_column_to_panel (op_pl3[i], false); 
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl3[i], "Blue",  GLUI_SPINNER_FLOAT, &gLightClr[i][2],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
	 
	}



	//glui_lighting->add_checkbox_to_panel (op_pl2, "Diffusion", &gBLight0[0], 105, init  );
	//glui_lighting->add_checkbox_to_panel (op_pl2, "Specular", &gBLight0[2], 107, init  );
	//glui_lighting->add_column_to_panel (op_pl2, false);
	//glui_lighting->add_checkbox_to_panel (op_pl2, "Ambient", &gBLight0[1] , 106, init  ); 

	GLUI_Panel *op_pl2 = glui_lighting->add_panel ("Light0"); 

	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl2, "Red",  GLUI_SPINNER_FLOAT, &gBLight0[0],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl2, "Green",  GLUI_SPINNER_FLOAT, &gBLight0[1],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
	glui_lighting->add_column_to_panel (op_pl2, false); 
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl2, "Blue",  GLUI_SPINNER_FLOAT, &gBLight0[2],-1,init);
	pSpin->set_float_limits(0.0, 1.0);
	pSpin->set_speed(40);
 
	 
	GLUI_Panel *op_pl5 = glui_lighting->add_panel ("Light0 Position"); 
	 
		pSpin=glui_lighting->add_spinner_to_panel(
		op_pl5, "X",  GLUI_SPINNER_FLOAT, &gfLightPos[0],-1,init);
	pSpin->set_float_limits(-20.0, 20.0);
	pSpin->set_speed(4 );
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl5, "Y",  GLUI_SPINNER_FLOAT, &gfLightPos[1],-1,init);
	pSpin->set_float_limits(-20.0, 20.0);
	pSpin->set_speed(4 ); 

	glui_lighting->add_column_to_panel (op_pl5, false); 
	pSpin=glui_lighting->add_spinner_to_panel(
		op_pl5, "Z",  GLUI_SPINNER_FLOAT, &gfLightPos[2],-1,init);
	pSpin->set_float_limits(-20.0, 20.0);
	pSpin->set_speed(4 );

	GLUI_Panel *op_pl7 = glui_lighting->add_panel ("Light1"); 
	glui_lighting->add_checkbox_to_panel (op_pl7, "Enable Light1", &gBLightOpen1 , -1, init  );
	GLUI_Spinner *pSp[6];
	pSp[0]=glui_lighting->add_spinner_to_panel(
		op_pl7, "Red",  GLUI_SPINNER_FLOAT, &gBLight1[0],-1,init);
	pSp[0]->set_float_limits(0.0, 1.0);
	pSp[0]->set_speed(40);
	glui_lighting->add_column_to_panel (op_pl7, false); 
	pSp[1]=glui_lighting->add_spinner_to_panel(
		op_pl7, "Green",  GLUI_SPINNER_FLOAT, &gBLight1[1],-1,init);
	pSp[1]->set_float_limits(0.0, 1.0);
	pSp[1]->set_speed(40);
	pSp[2]=glui_lighting->add_spinner_to_panel(
		op_pl7, "Blue",  GLUI_SPINNER_FLOAT, &gBLight1[2],-1,init);
	pSp[2]->set_float_limits(0.0, 1.0);
	pSp[2]->set_speed(40); 


	GLUI_Panel *op_pl6 = glui_lighting->add_panel ("Light1 Position"); 
	pSp[3]=glui_lighting->add_spinner_to_panel(
		op_pl6, "X",  GLUI_SPINNER_FLOAT, &gfLight1Pos[0],-1,init);
	pSp[3]->set_float_limits(-20.0, 20.0);
	pSp[3]->set_speed(4 ); 
	pSp[4]=glui_lighting->add_spinner_to_panel(
		op_pl6, "Y",  GLUI_SPINNER_FLOAT, &gfLight1Pos[1],-1,init);
	pSp[4]->set_float_limits(-20.0, 20.0);
	pSp[4]->set_speed(4 ); 

	glui_lighting->add_column_to_panel (op_pl6, false); 
	pSp[5]=glui_lighting->add_spinner_to_panel(
		op_pl6, "Z",  GLUI_SPINNER_FLOAT, &gfLight1Pos[2],-1,init);
	pSp[5]->set_float_limits(-20.0, 20.0);
	pSp[5]->set_speed(4 ); 
	for(int i=0; i<6; ++i) pSp[i]->disable();

	glui_lighting->disable();
	glui_lighting->set_main_gfx_window( main_window );
}

 

int main(int argc, char** argv)
{


	window_x = (glutGet (GLUT_SCREEN_WIDTH) - 640)/2;
	window_y = (glutGet (GLUT_SCREEN_HEIGHT) - 640)/2;

	//MyObj obj("coati.obj");
	MyObj obj("bunny.obj");
	MyObj obj2("cube.obj");
	pObjArr[0] = &obj;
	pObjArr[1] = &obj2;
	pObj = pObjArr[0];

	glutInit(&argc, argv); 
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);		// Setup display mode to double buffer and RGB color
	glutInitWindowSize(512,512);
	glutInitWindowPosition (window_x, window_y);						// Set the screen size
	main_window =  glutCreateWindow("Assignment2");	
	init ();
	glutReshapeFunc(reshape);		
	glutDisplayFunc(display);	
	glutKeyboardFunc(keyboard);		 
	glutMouseFunc(mouse);  
	glutMotionFunc(mouseMotionCB);
	setupGLUI(); 
	glutMainLoop();
	return 0;
}