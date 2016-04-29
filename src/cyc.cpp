/* 
	cpp source to compute circles and webs from given cyclides.

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
//#include <gl/GLAux.h> 

#include <iostream>
#include <fstream>

using namespace std;


#include <math.h>
using namespace std;
#pragma comment(linker,"/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#define PI 3.1415926
#define EPS 0.0000001
#define SQRT_2 sqrt(2.0)
#define SQRT_3 sqrt(3.0)

// For GUI
int window_x;
int window_y;
GLuint main_window;//  pointer to the GLUI window
GLUI * glui_window; 
void setupGLUI ();
int gbShowAxis=0;


bool mouseLeftDown;
bool mouseRightDown;
float mouseX, mouseY;
float cameraAngleX = 0;
float cameraAngleY = 0;
float cameraDistance;

float AngX=30.0f;
float AngY=0.0f;
float AngZ=0.0f;

float AngVx = 0.0f;
float AngVy = 0.0f;
float AngVz = 0.0f;
float gfScale = 1.5f;

int gbFlag = 0;
float len = 55.83; //length of the curve  
//10.4537873265 25
//len = 22.48f; great parameter 15?



float gfT = 0.04f;
float h = -3.2f; // parameter
int gbFlags[12] = {0, 0, 0, 0,0,0, 1,0,1,0,1,0}; //six families of circles
float gfOffSet = 0.111f; 
int gbInverse=0;

void CirclesOnMSpheres( float **pfCirData, float *  );
void CirclesOnMSpheres2( float **pfCirData  );
void CirclesOnMSpheres3( float **pfCirData  );
void CirclesOnMSpheres4( float **pfCirData  );
void CirclesOnMSpheres5( float **pfCirData  );
void CirclesOnMSpheres6( float **pfCirData  );

float thetax=3.15;//2.1;
float cirpara=0.07;
void OutF(ofstream* out, float *pos, float *pos2=NULL)
{
	if(pos2)
	{
		*out <<"cylinder { ";
	}
	*out<<"<"<<pos[0]<<", "<<pos[1]<<", "<<-pos[2]<<">";

	if(pos2)
	{
		*out<<", <"<<pos2[0]<<", "<<pos2[1]<<", "<<-pos2[2]<<">";
		*out <<", r_edge texture { edge_tex } } "<<endl;
	}
}

void init()
{

	glClearColor (0.0, 0.0, 0.0, 0.0);			// Clear the color 
	glShadeModel (GL_FLAT);						// Set the shading model to GL_FLAT
	glEnable (GL_LINE_SMOOTH);
	glEnable(GL_NORMALIZE);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);		// Set Line Antialiasing

	glShadeModel(GL_FLAT);

}

float FindSquare(float *p)
{
	return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}
void CrossProduct(float *v1, float *v2, float *v3)
{
	v3[0]=v1[1]*v2[2] - v1[2]*v2[1];
	v3[1]=v1[2]*v2[0] - v1[0]*v2[2];
	v3[2]=v1[0]*v2[1] - v1[1]*v2[0];
	float dis= sqrt(FindSquare(v3));
	for(int i=0; i<3; ++i) v3[i]/=dis;
}


float DotProduct(float *v1, float *v2)
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

float FindDis3D(float *p1, float*p2)
{
	return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
}

void MShperes(float data[])
{
	float h, lamda, a, b ,c ,d  ;

	h = data[0]; //parameter;
	lamda = data[1];
	a = data[2];
	b = data[3];
	c = data[4];
	d = data[5];

	//a/=lamda; b/=lamda; c/=lamda; d/=lamda;
	float r = (a*a+b*b+c*c)/4 - d;
	r = sqrt(r);
	glPushMatrix();
	glTranslatef(-a/2.0, -b/2.0, -c/2.0);
	glutWireSphere(r, 20, 20); 

	glPopMatrix();

}

ofstream out("webs.inc");

void CirRotation(float *pCirPara, int num, float **data = NULL)
{

	float vect[3];
	float cen[3];
	float r;
	for(int i=0;i<3;++i) 
	{
		vect[i]=pCirPara[i];
		cen[i] = pCirPara[3+i];
	}
	r = pCirPara[6];

	if (vect[1]<-EPS) 
	{
		for(int i=0;i<3;++i) vect[i]=-vect[i];
	}
	float q1, q2;
	q1 = acos(vect[2]/sqrt(FindSquare(vect)));
	if (abs(vect[1]) < 0.0001)
	{
		q2 = PI/2;
		if (vect[0]<-EPS)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(vect[0]/vect[1]);

	glBegin(GL_LINE_LOOP);
	for (int i=0; i<num; ++i)
	{
		float t = 2.0*PI*i/num;

		float pos[3]  ={  cen[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ,
			cen[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)),
			  cen[2]+ r*sin(q1)*sin(t)};
 
				glVertex3fv(pos); 

					if(gbFlag)
					{
						t = 2.0*PI*(1+i)/num;

						float pos2[3]  ={  cen[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ,
							cen[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)),
							cen[2]+ r*sin(q1)*sin(t)};

						OutF(&out, pos, pos2);

					}


	}
	glEnd();

    
}



//Circle data:
//vect, center, radius

void CirclesByPoint(float *pCurData, float *pCirData, int idx=0, float *pCirPara = NULL)
{

	float lamda= -h /SQRT_2;
	float data[4]={SQRT_2, -h, 1,h/SQRT_2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = - SQRT_2/h/2.0; 
	float data2[4]={SQRT_2, 1.0f/h, -1.0f, SQRT_2/h/2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;  
	
	if(idx ==2)
	{
		//1
		lamda= -h /SQRT_2;
		data[0]=SQRT_2;
		data[1]=-h;
		data[2]=1;
		data[3]=h/SQRT_2; 
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda = - SQRT_2/h/2.0; 
		data2[0]=SQRT_2;
		data2[1]=1.0f/h;
		data2[2]=-1.0f;
		data2[3]=SQRT_2/h/2.0;
		for(int i=0;i<4;++i) data2[i]/=lamda;  

	}
	else if (idx==1)
	{//2
		float lamda= -h /SQRT_2; 
		data[0]=SQRT_2;
		data[1]=h;
		data[2]=1;
		data[3]=h/SQRT_2;  
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda = - SQRT_2/h/2.0;  
		data2[0]=SQRT_2;
		data2[1]=-1.0f/h;
		data2[2]=-1.0f;
		data2[3]=SQRT_2/h/2.0; 
		for(int i=0;i<4;++i) data2[i]/=lamda;    
	}
	else if(idx==3)
	{
		//3
		float lamda= -(SQRT_3+1)*h / 2.0; 
		data[0]=1;
		data[1]=0;
		data[2]=SQRT_2;
		data[3]=(SQRT_3-1)*h/2;  
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda = -(SQRT_3-1) /h/ 2.0; 
		data2[0]=1;
		data2[1]=0;
		data2[2]=-SQRT_2;
		data2[3]=(SQRT_3+1)/h/2; 
		for(int i=0;i<4;++i) data2[i]/=lamda;  

	}
	else if(idx==4)
	{
		//4
		float lamda=  (1-SQRT_3 )*h / 2.0; 
		data[0]=1;
		data[1]=0;
		data[2]=SQRT_2;
		data[3]=(SQRT_3+1)*h/2;  
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda = -(SQRT_3+1) /h/ 2.0; 
		data2[0]=1;
		data2[1]=0;
		data2[2]=-SQRT_2;
		data2[3]=(SQRT_3-1)/h/2; 
		for(int i=0;i<4;++i) data2[i]/=lamda;      

	}
	else if(idx==6)
	{
		//5
		float lamda=  (1-h ) / 2.0; 
		data[0]=SQRT_3;
		data[1]=-SQRT_2*h;
		data[2]=0;
		data[3]=-(1+h)/2; 
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda =-(1.0/h+1)/2; 
		data2[0]=SQRT_3;
		data2[1]=SQRT_2/h;
		data2[2]=0;
		data2[3]=(  1.0 -1.0 /h)/ 2.0; 
		for(int i=0;i<4;++i) data2[i]/=lamda;    

	}
	else if(idx==5)
	{
		//6
		lamda=  (1+h)/2;
		data[0]=SQRT_3;
		data[1]=-SQRT_2*h;
		data[2]=0;
		data[3]=-(1-h)/2; 
		for(int i=0;i<4;++i) data[i]/=lamda; 
		lamda = -(  1.0 -1.0 /h)/ 2.0; 
		data2[0]=SQRT_3;
		data2[1]=SQRT_2/h;
		data2[2]=0;
		data2[3]=(  1.0 +1.0 /h)/ 2.0; 
		for(int i=0;i<4;++i) data2[i]/=lamda;      
	}



	float m1[3]={-data[0]/2, -data[1]/2, -data[2]/2};
	float m2[3]={-data2[0]/2, -data2[1]/2, -data2[2]/2};//center
	float r1 = sqrt(FindSquare(data)/4 -data[3]);
	float r2 = sqrt(FindSquare(data2)/4-data2[3]); //radius
	float D =  FindDis3D(m1, m2 );
	float cosTheta=(r1*r1+D*D-r2*r2)/r1/D/2;
	float D1=r1*cosTheta;
	float r=sqrt(r1*r1-D1*D1);//radius of intersected circle of the two spheres
	float alph = D1/D;
	float m0[3] = {alph * m2[0] + (1-alph)*m1[0], alph * m2[1] + (1-alph)*m1[1],alph * m2[2] + (1-alph)*m1[2]}; //center of the circle


	float Vect[3]={data2[0]-data[0], data2[1]-data[1], data2[2]-data[2]};
	if (Vect[1]<0.00001) for(int i=0;i<3;++i) Vect[i]=-Vect[i];

	if(pCirPara)
	{
		for(int i=0;i<3;++i) 
		{
			pCirPara[i] = Vect[i];
			pCirPara[i+3] = m0[i];
		}
		pCirPara[6] = r;

	}

	float q1, q2;
	q1 = acos(Vect[2]/sqrt(FindSquare(Vect)));
	if (abs(Vect[1]) < EPS)
	{
		q2 = PI/2;
		if (Vect[0]<-EPS)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(Vect[0]/Vect[1]);


	glBegin(GL_LINE_LOOP);
	for (int i = 0; i<40; ++i)
	{  
		float t = PI*2*i/40;
		pCirData[3*i] =  m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ;
		pCirData[3*i+1] = m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t));
		pCirData[3*i+2] = m0[2]+ r*sin(q1)*sin(t);

			// glVertex3fv(&pCirData[3*i]);
		 
	}

	glEnd(); 



	//////////////////////////////////////////////////////////////////////////
	//finding the right parameter for the point;
	float x,y,z ;
	x = pCurData[1];
	y = pCurData[2]; 
	z = pCurData[3];

	float t0 = PI- asin((z-m0[2])/r/sin(q1));
	float pos[3] = {m0[0]+ r*(cos(q2)*cos(t0)-sin(q2)*cos(q1)*sin(t0) ),
	m0[1]- r*(sin(q2)*cos(t0)+cos(q2)*cos(q1)*sin(t0)),
	m0[2]+ r*sin(q1)*sin(t0)};
	t0 = PI-t0;
	float pos2[3] = {m0[0]+ r*(cos(q2)*cos(t0)-sin(q2)*cos(q1)*sin(t0) ),
		m0[1]- r*(sin(q2)*cos(t0)+cos(q2)*cos(q1)*sin(t0)),
		m0[2]+ r*sin(q1)*sin(t0)};

	if (FindDis3D(&pCurData[1], pos) < FindDis3D(&pCurData[1], pos2))
	{
		t0  = PI - t0;
	}

	// end
	//////////////////////////////////////////////////////////////////////////


	float theta =pCurData[0]/r;
	glColor3f(1.0, 0.9, 0.9);
	glBegin(GL_LINE_STRIP); 



	int num = 10;
	for (int i=0; i<= 10; ++i)
	{
		float t=-theta*i/num + t0;
		pCurData[3*i] = m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
		pCurData[3*i+1] = m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t));
		pCurData[3*i+2] = m0[2]+ r*sin(q1)*sin(t);
		//glVertex3fv(  &pCurData[3*i] 	);
	} 
	float t = theta + t0;
	pCurData[33] = m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
	pCurData[34] = m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t));
	pCurData[35] = m0[2]+ r*sin(q1)*sin(t);

	glEnd();

}

void DrawCurve( float *pos, int idx=1, float *ppCirPara = NULL, int bExtend = 0)
{
	float x, y, z;
	x = pos[0]; y = pos[1]; z = pos[2];
	float a, b, c;
	if(idx==3)
	{

		//family 3
		a = (2*x-2*SQRT_2*z)/(SQRT_3-1);
		b = 2*SQRT_3;
		c = -(2*x+2*SQRT_2*z)/(SQRT_3+1);
	}
	else if(idx==6)
	{

		//family 5
		a = SQRT_3*x+SQRT_2*y+1;
		b = 2*SQRT_2*y - 2*SQRT_3*x;
		c = -(SQRT_3*x+SQRT_2*y-1);
	}
	else if(idx==2)
	{

		//family 1
		a = 2*x-SQRT_2*z;
		b = 2*SQRT_2*y;
		c = -(2*x+SQRT_2*z);
	}
	else if(idx==1)
	{
		//family 2

		a = 2*x-SQRT_2*z;
		b = -2*SQRT_2*y;
		c = -(2*x+SQRT_2*z);
	}
	else if (4==idx)
	{

		//family 4
		a = (2*x-2*SQRT_2*z)/(SQRT_3+1);
		b = -2*SQRT_3;
		c = -(2*x+2*SQRT_2*z)/(SQRT_3-1);
	}
	else if(5==idx)
	{

		//family 6
		a = -SQRT_3*x+SQRT_2*y-1;
		b = -2*SQRT_3*x - 2*SQRT_2*y;
		c = SQRT_3*x-SQRT_2*y-1;

	}


	//looking for next point
	float pos2[3];
	float pos3[3];

	float pCurData[36];
	float pCirData[120];
	float pCurData2[36];
	float pCirData2[120];
	float pCirPara[7];
	float pCirPara2[7];
	pCurData[0] = len/100 ;
	pCurData[1] = x;
	pCurData[2] = y;
	pCurData[3] = z;
	pCurData2[0] = len/100 ;
	pCurData2[1] = x;
	pCurData2[2] = y;
	pCurData2[3] = z;
	h = -b-sqrt(b*b-4*a*c); h/=2*a;

	CirclesByPoint(pCurData, pCirData, idx, pCirPara);
	h = -b+sqrt(b*b-4*a*c); h/=2*a;
	CirclesByPoint(pCurData2, pCirData2, idx, pCirPara2);  

	if (FindDis3D(pCurData, pos) < FindDis3D(pCurData2, pos))
	{
		glColor3f(0.10, 0.9, 0.09);
		glBegin(GL_LINE_LOOP);

		for (int i=0; i<40; ++i)
		{
			//glVertex3fv(&pCirData[3*i]);
		}

		glEnd();
		glColor3f(1.0, 0.9, 0.9);
		glBegin(GL_LINE_STRIP); 
		for (int i=0; i<=10; ++i)
		{
			//glVertex3fv(&pCurData[3*i]);
		}

		glEnd();

		pos2[0] = pCirData[30];
		pos2[1] = pCirData[31];
		pos2[2] = pCirData[32];
		pos3[0] = pCirData[33];
		pos3[1] = pCirData[34];
		pos3[2] = pCirData[35];

		if (ppCirPara)
		{
			memcpy(ppCirPara, pCirPara, sizeof(float)*7);
		}

	}
	else
	{
		glColor3f(0.10, 0.9, 0.09);
		glBegin(GL_LINE_LOOP);

		for (int i=0; i<40; ++i)
		{
			//glVertex3fv(&pCirData2[3*i]);
		}

		glEnd();
		glColor3f(1.0, 0.9, 0.9);
		glBegin(GL_LINE_STRIP); 
		for (int i=0; i<=10; ++i)
		{
			//glVertex3fv(&pCurData2[3*i]);
		}

		glEnd();

		pos2[0] = pCurData2[30];
		pos2[1] = pCurData2[31];
		pos2[2] = pCurData2[32];
		pos3[0] = pCurData2[33];
		pos3[1] = pCurData2[34];
		pos3[2] = pCurData2[35];
		
		if (ppCirPara)
		{
			memcpy(ppCirPara, pCirPara2, sizeof(float)*7);
		}


	} 

	//pos2 - POINT B

	if(bExtend)
	{
		x = pos2[0]; y = pos2[1]; z = pos2[2];

		DrawCurve(pos2, 3, ppCirPara+7);



		DrawCurve(pos2, 5, ppCirPara+7*2);

		//x = pos3[0]; y = pos3[1]; z = pos3[2];
		////family 3
		//a = (2*x-2*SQRT_2*z)/(SQRT_3-1);
		//b = 2*SQRT_3;
		//c = -(2*x+2*SQRT_2*z)/(SQRT_3+1);

		//DrawCurve(a, b, c, pos3, 3);


		////family 5
		//a = SQRT_3*x+SQRT_2*y+1;
		//b = 2*SQRT_2*y - 2*SQRT_3*x;
		//c = -(SQRT_3*x+SQRT_2*y-1);

		//DrawCurve(a, b, c, pos3, 5);
	}



}

void Interpt2Cirs(float *pCirPara1, float *pCirPara2, float *pos)
{
	//two planes: fing intersection lines
	//找到交线
	float vect[3];
	CrossProduct(pCirPara1, pCirPara2, vect);
	float a1,b1,c1,d1,a2,b2,c2,d2;
	a1= pCirPara1[0]; b1 = pCirPara1[1]; c1 = pCirPara1[2]; d1 = -DotProduct(pCirPara1, pCirPara1+3);
	a2= pCirPara2[0]; b2 = pCirPara2[1]; c2 = pCirPara2[2]; d2 = -DotProduct(pCirPara2, pCirPara2+3);
	float x0,  y0,  z0;

	if(abs(b1) > EPS)
	{
		x0 = 0;
		if(abs(b2)<EPS)
		{
			z0 = - d2/c2;
		}
		else
			z0 = (d2 - d1*b2/b1)/(c1*b2/b1  - c2);
		y0 = -c1*z0/b1 - d1/b1;

	}
	else
	{
		x0 = 0;
		z0 = -d1/c1;
		y0 = -(c2*z0+d2)/b2;

		y0 = 0;
		z0 = (a2*d1/a1 - d2)/(c2-a2*c1/a1);
		x0 = -c1*z0/a1 - d1/a1;
	}

	float a = 1.0;
	float b = 2*(vect[0]*(x0-pCirPara1[3]) + vect[1]*(y0-pCirPara1[4]) + vect[2]*(z0-pCirPara1[5]));
	float c = (x0-pCirPara1[3])*(x0-pCirPara1[3]) + (y0-pCirPara1[4])*(y0-pCirPara1[4]) + (z0-pCirPara1[5])*(z0-pCirPara1[5]) - pCirPara1[6]*pCirPara1[6];


	float t1, t2;
	t1 = -b + sqrt(b*b-4*a*c); t1 /= 2*a;
	t2 = -b - sqrt(b*b-4*a*c); t2 /= 2*a;

	float pos1[3] = {x0+t1*vect[0], y0 + t1*vect[1], z0 + t1*vect[2]};
	float pos2[3] = {x0+t2*vect[0], y0 + t2*vect[1], z0 + t2*vect[2]};

	float dis1, dis2;
	dis1 = FindDis3D(pos1, pCirPara2+3);
	dis2 = FindDis3D(pos2, pCirPara2+3);

	if (abs(pCirPara2[6]-dis2) <abs(pCirPara2[6]-dis1))
	{
		memcpy(pos, pos2, sizeof(float)*3);
	}
	else
		memcpy(pos, pos1, sizeof(float)*3);



}



//hexogals
void FindCirclebyPt(float *pos)
{


	float x, y, z;

	//point A
	x = pos[0]; y = pos[1]; z = pos[2];

	float a, b, c;



	float ppCirPara[7*3];
	//family 1
	if(gbFlags[6])
	{
		DrawCurve(pos, 1, ppCirPara, 1);

		//point B
		glColor3f(0.2, 0.8, 0.3);
		CirRotation(ppCirPara, 50);
		CirRotation(ppCirPara+7, 50);
		glColor3f(0.9, 0.9, 0.1);
		CirRotation(ppCirPara+14, 50);
 


	} 
	////family 3
	float ppCirPara3[7];
	float pCirPara1A[7];
	if(gbFlags[8])
	{
 
		DrawCurve(pos, 3, ppCirPara3);
		CirRotation(ppCirPara3, 50);

		//point C
		float pos2[3];
		Interpt2Cirs(ppCirPara+14, ppCirPara3, pos2);

		DrawCurve(pos, 1, pCirPara1A);
		glColor3f(0.2, 0.8, 0.3);
		CirRotation(pCirPara1A, 50);

  
	}



	//family 5
	float ppCirPara5[7];
	float pCirPara2A[7];
	if(gbFlags[10])
	{
		DrawCurve(pos, 5,ppCirPara5);
		glColor3f(0.9, 0.9, 0.1);
		CirRotation(ppCirPara5, 50);


		//point D
		float pos2[3];
		Interpt2Cirs(ppCirPara+7, ppCirPara5, pos2);

		DrawCurve(pos, 1, pCirPara2A);
		glColor3f(0.2, 0.8, 0.3);
		CirRotation(pCirPara2A, 50);
	}


	//point E
	float pos3[3];
	Interpt2Cirs(pCirPara1A, ppCirPara5, pos3);
	float ppCirPara3B[7];
	DrawCurve(pos3, 3, ppCirPara3B);
	CirRotation(ppCirPara3B, 50);



	//point F
	float pos4[3];
	Interpt2Cirs(pCirPara2A, ppCirPara3, pos4);

   
	float ppCirPara5B[7];
	DrawCurve(pos4, 5, ppCirPara5B);
	glColor3f(0.9, 0.9, 0.1);
	CirRotation(ppCirPara5B, 50); 



	//family 2
	if(gbFlags[7])
	{

		DrawCurve(pos, 2);
	}
	//family 4
	if(gbFlags[9])
	{
 
		DrawCurve(pos, 4);
	}



	//family 6
	if(gbFlags[11])
	{
		DrawCurve(pos, 6);
	}


}
int nMode = 25;
int nMode2 = nMode  ;
void FindCirclebyPt3B(float *pos, float *pos2, int mode)
{
	float clr[3][3] = {0.2, 0.8, 0.3, 
		0.9, 0.9, 0.1,
		0.8, 0.9, 0.8};  
	float pCirParaA5[7];
	float pCirParaB1[7];
	float pCirParaB5[7];
	float pCirParaC3[7];

	DrawCurve(pos, 5, pCirParaA5);
	DrawCurve(pos2, 1, pCirParaB1);
	
	//A', A5 v B1 (A'=C)
	float pos3[3];
	Interpt2Cirs(pCirParaA5, pCirParaB1, pos3);
	DrawCurve(pos3, 3, pCirParaC3);
	glColor3fv(clr[1]);
	if(gbFlag)
		out<<"#declare edge_tex=my_color_texture_3;" <<endl;
	CirRotation(pCirParaC3, 50);

	DrawCurve(pos2, 5, pCirParaB5);
	//B', B5 v A'3
	float pos4[3];
	Interpt2Cirs(pCirParaB5, pCirParaC3, pos4);

	if (mode <nMode2 )
	{
		FindCirclebyPt3B(pos3, pos4, 1+mode);
	}


}

void FindCirclebyPt3(float *pos, float *pos2, int mode)
{
	float clr[3][3] = {0.2, 0.8, 0.3, 
		0.9, 0.9, 0.1,
		0.8, 0.9, 0.8};  
	float pCirParaA1[7];
	float pCirParaA3[7];
	float pCirParaA5[7];

	float pCirParaB3[7];
	float pCirParaB5[7];

	float pCirParaC3[7];
	float pCirParaD5[7]; 

	float pos3[3]; //B'
	float pos4[3]; //A'

	//B', A1 v B5
	DrawCurve(pos, 1, pCirParaA1);
	DrawCurve(pos2, 5, pCirParaB5);

	Interpt2Cirs(pCirParaB5, pCirParaA1, pos3);
	
	DrawCurve(pos3, 3, pCirParaC3);
	glColor3fv(clr[1]);
	if(gbFlag)
		out<<"#declare edge_tex=my_color_texture_3;" <<endl;
 
	CirRotation(pCirParaC3, 50);

	//A', A5 v C3 (C=B')
	DrawCurve(pos, 5, pCirParaA5);
	Interpt2Cirs(pCirParaC3, pCirParaA5, pos4); 
	if (mode <nMode )
	{
		FindCirclebyPt3(pos4, pos3, 1+mode);
	}

}

// another direction of pt2
void FindCirclebyPt2B(float *pos, float *pos2, int mode)
{

	float clr[3][3] = {0.2, 0.8, 0.3, 
		0.9, 0.9, 0.1,
		0.8, 0.9, 0.8};  

	float pCirParaA1[7];
	float pCirParaA3[7];
	float pCirParaA5[7];

	float pCirParaB1[7];
	float pCirParaB3[7];

	float pCirParaC1[7];
	float pCirParaD5[7]; 

	//point A,  circle5;
	DrawCurve(pos, 5, pCirParaA5);
	glColor3fv(clr[2]);
	if(gbFlag) 
	{
 
		out<<"#declare edge_tex=my_color_texture_7;" <<endl;
	}
	CirRotation(pCirParaA5, 50); 

	//point B, circle 1
	DrawCurve(pos2, 1, pCirParaB1);
	glColor3fv(clr[1]);
	if(gbFlag)
		out<<"#declare edge_tex=my_color_texture_2;" <<endl;
	CirRotation(pCirParaB1, 50);

	// for next vertex:
	//A', A3 v B1
	float pos3[3];
	DrawCurve(pos, 3, pCirParaA3);
	Interpt2Cirs(pCirParaA3, pCirParaB1, pos3);
	//B', B3 v A'5 (A'=D)
	float pos4[3];
	DrawCurve(pos2, 3, pCirParaB3);
	DrawCurve(pos3, 5, pCirParaD5);
	Interpt2Cirs(pCirParaD5, pCirParaB3, pos4);
 
	if (mode <nMode2)
	{
		FindCirclebyPt2B(pos3, pos4, 1+ mode);
	}


}

//hexogals in after initalization
void FindCirclebyPt2(float *pos, float *pos2, int mode, float **dataA = NULL, float **dataB=NULL)
{

	float clr[3][3] = {0.2, 0.8, 0.3, 
		0.9, 0.9, 0.1,
		0.8, 0.9, 0.8};  


	float pCirParaA1[7];
	float pCirParaA3[7];
	float pCirParaA5[7];

	float pCirParaB3[7];
	float pCirParaB5[7];

	float pCirParaC1[7];
	float pCirParaD5[7]; 
 
		//point A
		DrawCurve(pos, 1, pCirParaA1); 
		glColor3fv(clr[0]);
		if(gbFlag)
			out<<"#declare edge_tex=my_color_texture_2;" <<endl;
		 CirRotation(pCirParaA1, 50); 
		 if(dataA && dataB)
		 memcpy(dataA[mode], pCirParaA1, sizeof(float)*7);



		DrawCurve(pos, 3, pCirParaA3); //*

		//point B
		DrawCurve(pos2, 3, pCirParaB3); //* 

		float pos5[3]; //B'=D
		Interpt2Cirs(pCirParaB3, pCirParaA1, pos5);

		//point D
		DrawCurve(pos5, 5, pCirParaD5);
		glColor3fv(clr[2]);
		if(gbFlag) 
			out<<"#declare edge_tex=my_color_texture_4;" <<endl;
		CirRotation(pCirParaD5, 50); 
		if(dataA && dataB)
		memcpy(dataB[mode], pCirParaD5, sizeof(float)*7);


		float pos4[3]; //A'=C

		Interpt2Cirs(pCirParaD5, pCirParaA3, pos4); 
 


	if (mode <nMode+1)
	{
		               // A',    B'
		FindCirclebyPt2(pos4,   pos5, 1+mode, dataA, dataB);
	}



}


//hexogals in one direction
//initailizing
void FindCirclebyPt(float *pos, int mode)
{

	float clr[3][3] = {0.2, 0.8, 0.3, 
					0.9, 0.9, 0.1,
					0.8, 0.9, 0.8};

	float x, y, z;

	//point A
	x = pos[0]; y = pos[1]; z = pos[2];


	float ppCirPara[7*3];  //1, 3, 5 in point B

	float pCirParaA1[7];
	float pCirParaA3[7];
	float pCirParaA5[7];

	float *pCirParaB1 = pCirParaA1;
	float pCirParaB3[7];
	float pCirParaB5[7];

	float pCirParaC1[7];
	float pCirParaD1[7];

	float pCirParaE3[7];
	float pCirParaF5[7];
 
	//family 1
	if(gbFlags[6])
	{
		// the circle 1,3,5 at point B

		//point B
		DrawCurve(pos, 1, ppCirPara, 1);

		memcpy(pCirParaA1, ppCirPara, sizeof(float)*7);
		memcpy(pCirParaB3, ppCirPara+7, sizeof(float)*7);
		memcpy(pCirParaB5, ppCirPara+14, sizeof(float)*7); 
		glColor3fv(clr[0]);
		if(gbFlag)
			out<<"#declare edge_tex=my_color_texture_2;" <<endl;
		CirRotation(pCirParaA1, 50);
		glColor3fv(clr[1]);
		if(gbFlag)
			out<<"#declare edge_tex=my_color_texture_3;" <<endl;
		CirRotation(pCirParaB3, 50);
		glColor3fv(clr[2]);
		if(gbFlag)
			out<<"#declare edge_tex=my_color_texture_4;" <<endl;
		CirRotation(pCirParaB5, 50);   

	} 

	//family 3
	//point C, the intersection of circle 3 at A and circle 5 At B
	float pos2[3];
	if(gbFlags[8])
	{

		//circle 3 at point A
		DrawCurve(pos, 3, pCirParaA3);
 
		glColor3fv(clr[1]);
		if(gbFlag)
			out<<"#declare edge_tex=my_color_texture_3;" <<endl;
		CirRotation(pCirParaA3, 50); 

		//point C, the intersection of circle 3 at A and circle 5 At B
		//float pos2[3];
		Interpt2Cirs(pCirParaB5, pCirParaA3, pos2);

		//find circle 1 at point C 
		DrawCurve(pos, 1, pCirParaC1);
		glColor3fv(clr[0]);
		//CirRotation(pCirParaC1, 50);  
	} 

	//family 5
	float pos6[3];
	if(gbFlags[10])
	{

		//circle 5 at point A
		DrawCurve(pos, 5, pCirParaA5); 
		glColor3fv(clr[2]);
		//CirRotation(pCirParaA5, 50);  
		//point D, the intersection of circle 5 at A and circle 3 at B
		Interpt2Cirs(pCirParaB3, pCirParaA5, pos6);

		//find circle circle 1 at D
		DrawCurve(pos6, 1, pCirParaD1); 
		glColor3fv(clr[0]);
		//CirRotation(pCirParaD1, 50);
	}


		//point E
		float pos3[3];
		//point F
		float pos4[3];
	if(0 && gbFlags[10] && gbFlags[8])
	{
		Interpt2Cirs(pCirParaC1, pCirParaA5,  pos3);
		DrawCurve(pos3, 3, pCirParaE3);
		glColor3fv(clr[1]);
		CirRotation(pCirParaE3, 50);



		Interpt2Cirs( pCirParaD1, pCirParaA3, pos4);

		DrawCurve(pos4, 5, pCirParaF5);
		glColor3fv(clr[2]);
		CirRotation(pCirParaF5, 50);  

	}


	glBegin(GL_LINES);
	glVertex3fv(pos2);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, 0);
	glVertex3fv(pos6);
	glEnd();
 
 

	float pos5[3];
	Interpt2Cirs(pCirParaB3, pCirParaC1, pos5); 


	if(1)
	{
		// circle 1, circle5
		FindCirclebyPt2(pos2, pos5,  mode ); 
		//-              A,   D
		FindCirclebyPt2B(pos, pos6, mode); 

		// pos6[3]-D
		//circle 3 
		 FindCirclebyPt3(pos6, pos5, mode);
		 FindCirclebyPt3B(pos, pos2, mode);

	} 

}


void DrawCircles(float data[], float data2[], float **pfCirData=NULL, float *TwoPoints=NULL)
{



	float m1[3]={-data[0]/2, -data[1]/2, -data[2]/2};
	float m2[3]={-data2[0]/2, -data2[1]/2, -data2[2]/2};//center
	float r1 = sqrt(FindSquare(data)/4 -data[3]);
	float r2 = sqrt(FindSquare(data2)/4-data2[3]); //radius
	float D =  FindDis3D(m1, m2 );
	float cosTheta=(r1*r1+D*D-r2*r2)/r1/D/2;
	float D1=r1*cosTheta;
	float r=sqrt(r1*r1-D1*D1);//radius of intersected circle of the two spheres
	float alph = D1/D;
	float m0[3] = {alph * m2[0] + (1-alph)*m1[0], alph * m2[1] + (1-alph)*m1[1],alph * m2[2] + (1-alph)*m1[2]}; //center of the circle


	glColor3f(0.2, 0.8, 0.1);
	glBegin(GL_LINE_LOOP);
	int idx = 0; 
	float Vect[3]={data2[0]-data[0], data2[1]-data[1], data2[2]-data[2]};
	if (Vect[1]<0.00001) for(int i=0;i<3;++i) Vect[i]=-Vect[i];
	float q1, q2;
	q1 = acos(Vect[2]/sqrt(FindSquare(Vect)));
	if (abs(Vect[1]) < EPS)
	{
		q2 = PI/2;
		if (Vect[0]<-EPS)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(Vect[0]/Vect[1]);



	for (float t=0.0f; t<2*PI; t+=gfT)
	{  
		if (NULL != pfCirData)
		{
			pfCirData[idx][0] = m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
			pfCirData[idx][1] = m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t));
			pfCirData[idx][2] = m0[2]+ r*sin(q1)*sin(t); 
			++idx;
		}
		else
		{
			
			glVertex3f( m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ,
			m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)), 
			m0[2]+ r*sin(q1)*sin(t));
		}
	}

	glEnd(); 



	if (TwoPoints != NULL)
	{
		float t0 = PI- asin((TwoPoints[1]-m0[2])/r/sin(q1));
		float theta =TwoPoints[0]/r;
		glColor3f(1.0, 0.9, 0.9);
		glBegin(GL_LINE_STRIP); 
 


		int num = 10;
		for (int i=0; i<= 10; ++i)
		{
			float t=-theta*i/num + t0;
			glVertex3f( m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ,
				m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)), 
				m0[2]+ r*sin(q1)*sin(t));
		} 
		glEnd();
		float t = t0-theta;

		float NPoint[3] = {m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ) ,
			m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)), 
			m0[2]+ r*sin(q1)*sin(t)};
	}
}

void RandPoint(float *pos)
{
	if (cirpara < 0.25)
	{
		h = -1.0/4/cirpara;
	}
	else if(cirpara < 0.5)
	{
		h = -(0.5-cirpara)*4;
	}
	else if(cirpara <0.75)
	{
		h = (cirpara-0.5)*4;
	}
	else
	{
		h = 1.0/4/(1.0-cirpara);
	}
	float lamda=  (1-h ) / 2.0;
	float data[4]={ SQRT_3, -SQRT_2*h,  0,-(1+h)/2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda =-(1.0/h+1)/2;
	float data2[4]={SQRT_3, SQRT_2/h, 0, (  1.0 -1.0 /h)/ 2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;   


	float m1[3]={-data[0]/2, -data[1]/2, -data[2]/2};
	float m2[3]={-data2[0]/2, -data2[1]/2, -data2[2]/2};//center
	float r1 = sqrt(FindSquare(data)/4 -data[3]);
	float r2 = sqrt(FindSquare(data2)/4-data2[3]); //radius
	float D =  FindDis3D(m1, m2 );
	float cosTheta=(r1*r1+D*D-r2*r2)/r1/D/2;
	float D1=r1*cosTheta;
	float r=sqrt(r1*r1-D1*D1);//radius of intersected circle of the two spheres
	float alph = D1/D;
	float m0[3] = {alph * m2[0] + (1-alph)*m1[0], alph * m2[1] + (1-alph)*m1[1],alph * m2[2] + (1-alph)*m1[2]}; //center of the circle







	float t = thetax+PI;  
	if(cirpara < 0.396447 || cirpara > 0.896446)
		t-=PI;


	float t0=t;


	float t2 = -atan((data2[1]-data[1])/(data2[0]-data[0]));


	pos[0] = m0[0]+ r*cos(t0)*sin(t2);
	pos[1] = m0[1]+ r*cos(t0)*cos(t2);
	pos[2] = m0[2]+ r*sin(t0); 
			 

 
 

}

void CirclesOnMSpheres6( float **pfCirData  )
{   
	float lamda=  (1+h)/2;
	float data[4]={ SQRT_3, -SQRT_2*h,  0,-(1-h)/2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = -(  1.0 -1.0 /h)/ 2.0;
	float data2[4]={SQRT_3, SQRT_2/h, 0, (  1.0 +1.0 /h)/ 2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;       
	DrawCircles(data, data2, pfCirData);
}
void CirclesOnMSpheres5( float **pfCirData  )
{ 

	float lamda=  (1-h ) / 2.0;
	float data[4]={ SQRT_3, -SQRT_2*h,  0,-(1+h)/2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda =-(1.0/h+1)/2;
	float data2[4]={SQRT_3, SQRT_2/h, 0, (  1.0 -1.0 /h)/ 2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;      
	DrawCircles(data, data2, pfCirData);
}
void CirclesOnMSpheres4( float **pfCirData  )
{   
	float lamda=  (1-SQRT_3 )*h / 2.0;
	float data[4]={1, 0,  SQRT_2,(SQRT_3+1)*h/2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = -(SQRT_3+1) /h/ 2.0;
	float data2[4]={1, 0, -SQRT_2, (SQRT_3-1)/h/2};
	for(int i=0;i<4;++i) data2[i]/=lamda;      
	DrawCircles(data, data2, pfCirData);
}
void CirclesOnMSpheres3( float **pfCirData  )
{   
	float lamda= -(SQRT_3+1)*h / 2.0;
	float data[4]={1, 0,  SQRT_2,(SQRT_3-1)*h/2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = -(SQRT_3-1) /h/ 2.0;
	float data2[4]={1, 0, -SQRT_2, (SQRT_3+1)/h/2};
	for(int i=0;i<4;++i) data2[i]/=lamda;       
	DrawCircles(data, data2, pfCirData);
}
void CirclesOnMSpheres2( float **pfCirData  )
{  

	float lamda= -h /SQRT_2;
	float data[4]={SQRT_2, h, 1,h/SQRT_2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = - SQRT_2/h/2.0; 
	float data2[4]={SQRT_2, -1.0f/h, -1.0f, SQRT_2/h/2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;     
	DrawCircles(data, data2, pfCirData);

} 





void CirclesOnMSpheres( float **pfCirData, float *TwoPoint = NULL  )
{  
	float lamda= -h /SQRT_2;
	float data[4]={SQRT_2, -h, 1,h/SQRT_2} ; 
	for(int i=0;i<4;++i) data[i]/=lamda; 
	lamda = - SQRT_2/h/2.0; 
	float data2[4]={SQRT_2, 1.0f/h, -1.0f, SQRT_2/h/2.0};
	for(int i=0;i<4;++i) data2[i]/=lamda;    
	DrawCircles(data, data2, pfCirData, TwoPoint);
}


void Hyperbloid()
{
	glColor3f(0.7f, 0.0f, 0.7f);
	for (float x=0.0f; x<4.0f; x+=0.12f)
	{
		float a, b; //long, short axises.
		a=sqrt((1+x*x)/2);
		b=sqrt((1+x*x)/3);

		glBegin(GL_LINE_STRIP);
		for (float theta=0; theta < 2*PI ; theta += 0.05 )
		{
			glVertex3f(x, a*cos(theta), b*sin(theta));
		}

		glEnd();
	}

	glColor3f(0.0f, 0.7f, 0.6f);
	//glBegin(GL_LINE_STRIP);  // the smallest circle
	//for (float theta=0; theta < 2*PI ; theta += 0.05 )
	//{
	//	glVertex3f(0, sqrt(1.0/2.0f)*cos(theta), sqrt(1.0f/3.0f)*sin(theta));
	//}
	//glEnd();
	glBegin(GL_LINES); // rulings   
	for (float theta=0.01; theta < PI/2 ; theta += 0.15 )
	{

		//glVertex3f(0, sqrt(1.0/2.0f)*cos(theta), sqrt(1.0f/3.0f)*sin(theta));

		float theta2=theta+PI/3;
		glVertex3f(sqrt(3.0f), 2*sqrt(1.0/2.0f)*cos(theta2), 2*sqrt(1.0f/3.0f)*sin(theta2));
		theta2=theta-PI/3;
		glVertex3f(-sqrt(3.0f), 2*sqrt(1.0/2.0f)*cos(theta2), 2*sqrt(1.0f/3.0f)*sin(theta2));

		glVertex3f(sqrt(3.0f), 2*sqrt(1.0/2.0f)*cos(theta2), 2*sqrt(1.0f/3.0f)*sin(theta2));
		theta2=theta+PI/3;
		glVertex3f(-sqrt(3.0f), 2*sqrt(1.0/2.0f)*cos(theta2), 2*sqrt(1.0f/3.0f)*sin(theta2));

	}
	glEnd();


} 
void Inversion(float *pos, float *pos2)
{
	float fSphPos[3] = {2.0, 2.0, 2.0};
	float r = 1.5;

	//the inversion of xy, yz, xz planes

	//xy plane
	// 1.5*1.5/2.0
	float dis = FindDis3D(pos, fSphPos);
	float alp = r*r/dis/dis; 
	pos2[0] =  (1-alp)*fSphPos[0] + alp * pos[0] ;
	pos2[1] =  (1-alp)*fSphPos[1] + alp * pos[1] ;
	pos2[2] =  (1-alp)*fSphPos[2] + alp * pos[2] ;

}

void Cyclide6Fam()
{

	ofstream out("circles9.inc");
	ofstream out2("circles7.inc");
	//1:
	for (int k=0; k<6; ++k)
	{
 		if(0==gbFlags[k])
		{
			continue;
		} 
		float fClr[6][3] ={.7, 0.2, 0.75,.2, 0.7, 0.75,.7, 0.7, 0.15,
						.1, 0.7, 0.15,0.6, 0.2, 0.15, 0.1, 0.2, 0.65} ;
		glColor3fv(fClr[k]);
		//void* pFun[6]={CirclesOnMSpheres, CirclesOnMSpheres2, CirclesOnMSpheres3, 
		//	CirclesOnMSpheres4, CirclesOnMSpheres5, CirclesOnMSpheres6,};
 
		float fSet[2] = { 0.956f, 0.001f};
		int nCircles = int((fSet[0]-fSet[1])/gfOffSet)+1 ; // number of circles
 
		float ***pfCircleData = new float**[4*nCircles];
		int idx=0;

		int nCirData = int(2*PI/gfT) +1;  // data numbers of one circle
		for (float f=fSet[0];f>fSet[1];  f -= gfOffSet)
		{ 
			//++nCircles;
			h = f;

 			pfCircleData[3*nCircles-1-idx] = new float*[nCirData];//theta
			for (int i=0; i<nCirData; ++i)
			{
				pfCircleData[3*nCircles-1-idx][i] = new float[3]; //x,y,z
			} 

			if(int(nCircles*2.5) == 3*nCircles-1-idx)
			{
				float w = h;
			}
			//if(k==5 && idx >nCircles/2+3 )	gbInverse=1; //0.044
			if(k==5 && idx >nCircles/2  )	gbInverse=1; //0.084
 			if(0==k)
				CirclesOnMSpheres(pfCircleData[3*nCircles-1-idx]); 
			if(1==k)
				CirclesOnMSpheres2(pfCircleData[3*nCircles-1-idx]); 
			if(2==k)
				CirclesOnMSpheres3(pfCircleData[3*nCircles-1-idx]); 
			if(3==k)
				CirclesOnMSpheres4(pfCircleData[3*nCircles-1-idx]); 
			if(4==k)
				CirclesOnMSpheres5(pfCircleData[3*nCircles-1-idx]); 
			if(5==k)
				CirclesOnMSpheres6(pfCircleData[3*nCircles-1-idx]); 
			gbInverse=0;
 

			h = 1.0f/f;

			pfCircleData[3*nCircles+idx] = new float*[nCirData];//theta
			for (int i=0; i<nCirData; ++i)
			{
				pfCircleData[3*nCircles+idx][i] = new float[3]; //x,y,z
			}     
			//if(k==4 && idx >nCircles/2+3 )	gbInverse=1; //0.044
			if(k==4 && idx >nCircles/2  )	gbInverse=1; //0.084
			//if(k==5 && idx >nCircles/2+15 )	gbInverse=1; 
			if(0==k)
				CirclesOnMSpheres(pfCircleData[3*nCircles+idx]); 
			if(1==k)
				CirclesOnMSpheres2(pfCircleData[3*nCircles+idx]); 
			if(2==k)
				CirclesOnMSpheres3(pfCircleData[3*nCircles+idx]); 
			if(3==k)
				CirclesOnMSpheres4(pfCircleData[3*nCircles+idx]); 
			if(4==k)
				CirclesOnMSpheres5(pfCircleData[3*nCircles+idx]); 
			if(5==k)
				CirclesOnMSpheres6(pfCircleData[3*nCircles+idx]);   
			 gbInverse=0;

			h = -f;

			pfCircleData[nCircles+idx] = new float*[nCirData];//theta
			for (int i=0; i<nCirData; ++i)
			{
				pfCircleData[nCircles+idx][i] = new float[3]; //x,y,z
			} 

			//if(k==4 && idx <=nCircles/2+3    )	gbInverse=1;//0.084
			if(k==4 && idx <=nCircles/2    )	gbInverse=1;//0.044
			if(k==5  )	gbInverse=1;
			  
 			if(0==k)
				CirclesOnMSpheres(pfCircleData[nCircles+idx]); 
			if(1==k)
				CirclesOnMSpheres2(pfCircleData[nCircles+idx]); 
			if(2==k)
				CirclesOnMSpheres3(pfCircleData[nCircles+idx]); 
			if(3==k)
				CirclesOnMSpheres4(pfCircleData[nCircles+idx]); 
			if(4==k)
				CirclesOnMSpheres5(pfCircleData[nCircles+idx]); 
			if(5==k)
				CirclesOnMSpheres6(pfCircleData[nCircles+idx]);  
			gbInverse=0;


			h = -1.0f/f;

			pfCircleData[nCircles-1-idx] = new float*[nCirData];//theta
			for (int i=0; i<nCirData; ++i)
			{
				pfCircleData[nCircles-1-idx][i] = new float[3]; //x,y,z
			}
 

			if(k==4   )	gbInverse=1;
			//if(k==5 && idx <nCircles/2+4 )	gbInverse=1; // 0.044
			if(k==5 && idx <nCircles/2+1 )	gbInverse=1; // 0.084
			if(0==k)
				CirclesOnMSpheres(pfCircleData[nCircles-1-idx]); 
			if(1==k)
				CirclesOnMSpheres2(pfCircleData[nCircles-1-idx]); 
			if(2==k)
				CirclesOnMSpheres3(pfCircleData[nCircles-1-idx]); 
			if(3==k)
				CirclesOnMSpheres4(pfCircleData[nCircles-1-idx]); 
			if(4==k)
				CirclesOnMSpheres5(pfCircleData[nCircles-1-idx]); 
			if(5==k)
				CirclesOnMSpheres6(pfCircleData[nCircles-1-idx]); 
			gbInverse=0;
	 

			++idx;
		}


		glColor3fv(fClr[k]);
		for ( int i=0; i<4*nCircles; i+=1)
		{
			glBegin(GL_LINE_STRIP);

			for (int j=0; j<nCirData; ++j)
			{
				  glVertex3fv(pfCircleData[i][j]); 

				  if(gbFlag)
				  {
					  float pos2[3], pos3[3];
					  Inversion(pfCircleData[i][j], pos2);
					  Inversion(pfCircleData[i][(1+j)%nCirData], pos3);
					  OutF(&out2, pos2, pos3);
					  //OutF(&out,pfCircleData[i][j], pfCircleData[i][(1+j)%nCirData] );

				  }

			}
			glEnd();
		}

		if( k==0   )
		{ 
			
			//cylinder { <0.798145, 0.311514, 0.933333>, <0.811154, 0.234879, 0.933333>, r_edge texture { edge_tex } }
			ofstream txtOut("Plane.txt");
			txtOut <<"CData1=[ ";
			ofstream txt2("Cicle.inc");
			for (int j=0; j<nCirData; ++j)
			{
				txt2 <<"cylinder { <" << pfCircleData[ int( nCircles*2.5)][j][0] <<", " <<pfCircleData[ int( nCircles*2.5)][j][1]
				<<", "<<pfCircleData[ int( nCircles*2.5)][j][2] <<">, <" 
					<< pfCircleData[ int( nCircles*2.5)][(j+1)%nCirData][0] <<", " <<pfCircleData[ int( nCircles*2.5)][(j+1)%nCirData][1]
				<<", "<<pfCircleData[ int( nCircles*2.5)][(j+1)%nCirData][2] <<" >, r_edge texture { edge_tex } }"<<endl;
			}
			for (int i=0; i<3; ++i)
			{
				for (int j=0; j<nCirData; ++j)
					txtOut << pfCircleData[ int( nCircles*2.5)][j][i] << " ";
				txtOut <<";"<<endl;


			}
			txtOut <<"];"<<endl;
			glBegin(GL_LINE_LOOP);
			for (int j=0; j<nCirData; ++j)
			{

			//pfCircleData[2.5 * nCircles][i]
			glVertex3fv(pfCircleData[ int( nCircles*2.5)][j]); 

			}
			glEnd();
		}
		if( k==1 )
		{
			ofstream txtOut("Plane2.txt");
			txtOut <<"CData2=[ ";
			for (int i=0; i<3; ++i)
			{
				for (int j=0; j<nCirData; ++j)
					txtOut << pfCircleData[int( nCircles*1.5)-1][j][i] << " ";
				txtOut <<";"<<endl;

			}
			ofstream txt2("Cicle2.inc");
			for (int j=0; j<nCirData; ++j)
			{
				txt2 <<"cylinder { <" << pfCircleData[ int( nCircles*1.5)-1][j][0] <<", " <<pfCircleData[ int( nCircles*1.5)-1][j][1]
				<<", "<<pfCircleData[int( nCircles*1.5)-1][j][2] <<">, <" 
					<< pfCircleData[int( nCircles*1.5)-1][(j+1)%nCirData][0] <<", " <<pfCircleData[ int( nCircles*1.5)-1][(j+1)%nCirData][1]
				<<", "<<pfCircleData[int( nCircles*1.5)-1][(j+1)%nCirData][2] <<" >, r_edge texture { edge_tex } }"<<endl;
			}
			txtOut <<"];"<<endl;
			glBegin(GL_LINE_LOOP);
			for (int j=0; j<nCirData; ++j)
			{

				//pfCircleData[2.5 * nCircles][i]
				glVertex3fv(pfCircleData[int( nCircles*1.5)-1][j]); 

			}
			glEnd();

		}

		static int nflag[6]={1,1,1,1,1,1};
		if (0 &&nflag[k])
		{
			nflag[k]=0;
			char strName[255];
			sprintf(strName, "%dplus.obj", 1+k);
			char strEdges[255];
			sprintf(strEdges, "Edges%d.inc", 1+k);
			ofstream EdgeOut(strEdges);
			ofstream txtOut(strName);
			for ( int i=0; i<4*nCircles; ++i) //Vertex
			{ 

				for (int j=0; j<nCirData; ++j)
				{ 
					txtOut<<"v "<<pfCircleData[i][j][0] <<" "
						<<pfCircleData[i][j][1] <<" "
						<<pfCircleData[i][j][2] <<endl;

				} 
			}
			
			for ( int i=0; i<4*nCircles; ++i)  //faces
			{
				//if (i==1.5*nCircles || i==3.5*nCircles) continue;
				for ( int j=0; j<nCirData; ++j)
				{
					int data[4];
					data[0]=i*nCirData+j+1;
					data[1]=i*nCirData+(j+1)%nCirData+1;
					data[2]=((i+1)%(4*nCircles))*nCirData+j+1;
					data[3]=((i+1)%(4*nCircles))*nCirData+(j+1)%nCirData+1;
			 
						txtOut<<"f "<<data[1]<<" "
							<<data[2]<<" "
							<<data[0]<<endl<<"f "
							<<data[3]<<" "
							<<data[2]<<" "
							<<data[1]<<endl; 

						EdgeOut<<"cylinder { <"
							<<pfCircleData[i][j][0]<<", "
							<<pfCircleData[i][j][1]<<", "
							<<pfCircleData[i][j][2]<<">, <"
							<<pfCircleData[i][(j+1)%nCirData][0]<<", "
							<<pfCircleData[i][(j+1)%nCirData][0]<<", "
							<<pfCircleData[i][(j+1)%nCirData][0]
						<<">, r_edge texture { edge_tex } }"<<endl;


				} 
			} 
		}


		for (int i=0; i<4*nCircles; ++i)
		{
			for (int j=0; j<nCirData; ++j)
				delete []pfCircleData[i][j];
			delete []pfCircleData[i];
		}
		delete []pfCircleData; 



	}


	//2:
	//if(gbFlags[1])
	//{
	//	glColor3f(.2, 0.7, 0.75);
	//	for (float f=1.00f;f>0.001f;  f -= gfOffSet)
	//	{ 
	//		h = f;
	//		CirclesOnMSpheres2(); 
	//		h = 1.0f/f;
	//		CirclesOnMSpheres2();
	//		h = -f;
	//		CirclesOnMSpheres2(); 
	//		h = -1.0f/f;
	//		CirclesOnMSpheres2();
	//	}
	//}
	//3:

	//if(gbFlags[2])
	//{
	//	glColor3f(.7, 0.7, 0.15);
	//	for (float f=1.00f;f>0.001f;  f -= gfOffSet)
	//	{ 
	//		h = f;
	//		CirclesOnMSpheres3(); 
	//		h = 1.0f/f;
	//		CirclesOnMSpheres3();
	//		h = -f;
	//		CirclesOnMSpheres3(); 
	//		h = -1.0f/f;
	//		CirclesOnMSpheres3();
	//	}
	//}
	////4:


	//if(gbFlags[3])
	//{
	//	glColor3f(.1, 0.7, 0.15);
	//	for (float f=1.00f;f>0.001f;  f -= gfOffSet)
	//	{ 
	//		h = f;
	//		CirclesOnMSpheres4(); 
	//		h = 1.0f/f;
	//		CirclesOnMSpheres4();
	//		h = -f;
	//		CirclesOnMSpheres4(); 
	//		h = -1.0f/f;
	//		CirclesOnMSpheres4();
	//	}
	//}

	//if(gbFlags[4])
	//{
	//	glColor3f(0.6, 0.2, 0.15);
	//	for (float f=1.00f;f>0.001f;  f -= gfOffSet)
	//	{ 
	//		h = f;
	//		CirclesOnMSpheres5(); 
	//		h = 1.0f/f;
	//		CirclesOnMSpheres5();
	//		h = -f;
	//		CirclesOnMSpheres5(); 
	//		h = -1.0f/f;
	//		CirclesOnMSpheres5();
	//	}
	//}

	//if(gbFlags[5])
	//{
	//	glColor3f(0.1, 0.2, 0.65);
	//	for (float f=1.00f;f>0.001f;  f -= gfOffSet)
	//	{ 
	//		h = f;
	//		CirclesOnMSpheres6(); 
	//		h = 1.0f/f;
	//		CirclesOnMSpheres6();
	//		h = -f;
	//		CirclesOnMSpheres6(); 
	//		h = -1.0f/f;
	//		CirclesOnMSpheres6();
	//	}
	//}


	//glColor3f(205.0/255, 133.0/255, 63.0/255);
	//int num = 200;
 //

	//glBegin(GL_LINE_STRIP);
	//for (float h = 1+SQRT_3 ; h<2+2*SQRT_2; h+=0.02)
	//{ 
	//	float k = h; 
	//	float x = sqrt(k*k-2*k-2)/SQRT_2;
	//	float y = sqrt( 4+4*k-k*k)/SQRT_2;
	//	glVertex3f(x, y, 0);
	//}
	//glEnd();
	//glBegin(GL_LINE_STRIP);
	//for (float h = 1+SQRT_3 ; h<2+2*SQRT_2; h+=0.02)
	//{ 
	//	float k = h; 
	//	float x = sqrt(k*k-2*k-2)/SQRT_2;
	//	float y = -sqrt( 4+4*k-k*k)/SQRT_2;
	//	glVertex3f(x, y, 0);
	//}
	//glEnd();
	//glBegin(GL_LINE_STRIP);
	//for (float h = 1+SQRT_3 ; h<2+2*SQRT_2; h+=0.02)
	//{ 
	//	float k = h; 
	//	float x = -sqrt(k*k-2*k-2)/SQRT_2;
	//	float y = sqrt( 4+4*k-k*k)/SQRT_2;
	//	glVertex3f(x, y, 0);
	//}
	//glEnd();
	//glBegin(GL_LINE_STRIP);
	//for (float h = 1+SQRT_3 ; h<2+2*SQRT_2; h+=0.02)
	//{ 
	//	float k = h; 
	//	float x = -sqrt(k*k-2*k-2)/SQRT_2;
	//	float y = -sqrt( 4+4*k-k*k)/SQRT_2;
	//	glVertex3f(x, y, 0);
	//}
	//glEnd();


	 
 

	if(1)
	{
		float pos[3] ;
		RandPoint(pos);
		glColor3f(0.8, 0.1, 0.7);
		glBegin(GL_LINES);
		glVertex3fv(pos);
		glVertex3f(0.0, 0.0, 0.0);
		glEnd();


		glColor3f(0.0, 0.78, 0.1);
		FindCirclebyPt(pos, 0);

	}

}



void display(void)
{	 
	//init();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLineWidth(0.01);
	glPushMatrix();								 

	glScalef(5.0,5.0,5.0);
	glTranslatef(0, 0, cameraDistance);
	glRotatef(-90, 1, 0, 0);   // pitch
	glRotatef(cameraAngleX, 1, 0, 0);   // pitch
	glRotatef(cameraAngleY, 0, 1, 0);   // heading
	glRotatef(AngX,1.0,0.0,0.0);				// Rotate on x
	glRotatef(AngY,0.0,1.0,0.0);				// Rotate on y
	glRotatef(AngZ,0.0,0.0,1.0);				// Rotate on z

	if(gbShowAxis)
	{
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

	}

	glColor3f(1.0,1.0, 1.0);

	//glutWireSphere(3.2, 20, 20);
	glFlush();
	//auxSolidSphere(1.0);
	//
	//glFlush(); 

	glScalef(gfScale, gfScale, gfScale);

	//Hyperbloid();

	Cyclide6Fam();


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
	gluLookAt (AngVx, AngVy, 70.0 + AngVz, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); 
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
	case '1':

		h*=1.5f; 


		break;
	case '2':

		h /= 1.5f;


		break;
	case '-':
		gfScale *= 0.91;
		break;
	case '=':
		gfScale /= 0.91f;
		break;
	case '9':
		AngVz -= 1.0f;
		break;
	case '0':
		AngVz += 1.0f;
		break;	

	case 'x':		
		AngX -= 0.5f;
		break;  
	case 'X':	
		AngX += 0.5f;
		break; 
	case 'y':		
		AngY -= 0.5f;
		break;  
	case 'Y':	
		AngY += 0.5f;	
		break;	 
	case 'z':	
		AngZ -= 0.5f;
		break; 
	case 'Z':	
		AngZ += 0.5f;
		break;
	case '27':
		glutDestroyWindow(main_window);
		break;
	case 'n':
		AngVz = -50;
		break;
	case 'N':
		AngVz =  0;
		break;


	}

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(AngVx, AngVy, 70.0 + AngVz, 0,0,0,0,1,0);
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
	glui_window = GLUI_Master.create_glui ("Options", 0, window_x - 200, window_y);


	//---------------------------------------------------------------------
	// 'Object Properties' Panel
	//--------------------------------------------------------------------- 
	GLUI_Panel *op_panel = glui_window->add_panel (" ==Families of circles== ");

	//

	//  Add the Draw Check box to the 'Object Properties' Panel
	glui_window->add_checkbox_to_panel(op_panel, "Show Axis", &gbShowAxis );
	glui_window->add_checkbox_to_panel(op_panel, "Saving", &gbFlag );


	glui_window->add_checkbox_to_panel(op_panel, "Family 1", &gbFlags[0] );
	glui_window->add_checkbox_to_panel(op_panel, "Family 2", &gbFlags[1] );
	glui_window->add_checkbox_to_panel(op_panel, "Family 3", &gbFlags[2] );
	glui_window->add_checkbox_to_panel(op_panel, "Family 4", &gbFlags[3] );
	glui_window->add_checkbox_to_panel(op_panel, "Family 5", &gbFlags[4] );
	glui_window->add_checkbox_to_panel(op_panel, "Family 6", &gbFlags[5] );

	GLUI_Spinner *spinner =
		glui_window->add_spinner_to_panel(op_panel, "Circle Offset",  GLUI_SPINNER_FLOAT, &gfOffSet);

	spinner->set_float_limits(0.01, 0.5);


	glui_window->add_checkbox_to_panel(op_panel, "Circle 1", &gbFlags[0+6] );
	glui_window->add_checkbox_to_panel(op_panel, "Circle 2", &gbFlags[1+6] );
	glui_window->add_checkbox_to_panel(op_panel, "Circle 3", &gbFlags[2+6] );
	glui_window->add_checkbox_to_panel(op_panel, "Circle 4", &gbFlags[3+6] );
	glui_window->add_checkbox_to_panel(op_panel, "Circle 5", &gbFlags[4+6] );
	glui_window->add_checkbox_to_panel(op_panel, "Circle 6", &gbFlags[5+6] );


	//
	//float thetax=0;
	//float cirpara=0.1;
	GLUI_Spinner *spinner2 =
		glui_window->add_spinner_to_panel(op_panel, "vect y/x",  GLUI_SPINNER_FLOAT, &thetax); 
	spinner2->set_float_limits(0.0, PI*2);
	GLUI_Spinner *spinner3 =
		glui_window->add_spinner_to_panel(op_panel, "vect z/x",  GLUI_SPINNER_FLOAT, &cirpara); 
	spinner3->set_float_limits(0.00000001, 0.9999999); 

	//len
	GLUI_Spinner *spinner4 =
	glui_window->add_spinner_to_panel(op_panel, "Curve len",  GLUI_SPINNER_FLOAT, &len); 
	spinner4->set_float_limits(10, 99.99); 


	glui_window->set_main_gfx_window( main_window );
}

int main(int argc, char** argv)
{


	window_x = (glutGet (GLUT_SCREEN_WIDTH) - 400)/2;
	window_y = (glutGet (GLUT_SCREEN_HEIGHT) - 600)/2;

	glutInit(&argc, argv); 
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);		 
	glutInitWindowSize(1024,768);
	glutInitWindowPosition (window_x, window_y);						// Set the screen size
	main_window =  glutCreateWindow("Cyclides");	
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
