/* 
	cpp source to design cyclides(2 families of circles)

	authour: 'ling.shi'
	email: ling.shi@kaust.edu.sa; ling.shi@aliyun.com
	created date: 2011.4
	publish date: 2016.4
*/

#include <iostream>
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
#pragma comment(linker,"/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#define PI 3.1415926

// For GUI
int window_x;
int window_y;
GLuint main_window;//  pointer to the GLUI window
GLUI * glui_window; 
void setupGLUI ();

// For mouse contrl
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
float gfScale = 1.0f;

// For cyclides

// displaying colors
float gfCirPos1[3] = {0.664, 1.324, 0.0f}; 
float gfCirPos2[3] = {0.322, 0.054, 0.0f}; 
float gfCirPos3[3] =  {-0.3, -0.8, 0};
float gfTheta[3] = {-0.018, 1.89, 1.64};
float gfRad[3] = {0.83 , 0.8, 0.6};


//Sphere
float gfSphR = 17.0f;
float gfSphPos[3] = { -2.5,  -3.7, 0.7};
int gbShowAxis=0;
int gbShowSph=0;

//quadratics
float gfDis = 4.10f; //height
float gfAngStep = PI/20;//Angle around the circle in quadratic
float gfAngOffSet = PI/1.7;//Angle offset
int gbCrossLine = 1;//cross lines in quadratic
int gbCrossCyc=0; //tangent lines to sphere
int gbCirCyc=1; //tangent lines to sphere


//For Inverstion
int gbShowSp2 = 0;

//flags for Inversion
int gbFlags[10]={1, 0, 0, 0,0,0,  0,    0,0,  1};
static int nFlag=0; //saving data to file
float gfSphdata[2] = {2.0, 1.5f};

float gfr2 = 2.0f;
float fAlpha = 0.23f;  //define where this point lie along the Line 0~1
float t0 = 1.5*PI+0.44; //define which line to select
float gLightClr[3][4] = 	{
	0.9 ,0.2,0.2, 0.0,
	0.2, 0.9,0.2, 0.0,
	0.0,0.0, 1.0, 0.0};


float P1[3] = {  -0.9, 1.0, 0.8};
float P2[3] = {0.96, 0.98, 0.7};
float P3[3] = {1.56, -0.3, 1.38};
float P4[3] = {0.02, -0.8, 0.92};

// the further point
float P[3] = {0.53, 0.05, 1.43};


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


	//light
	GLfloat mat_specular[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess[] = { 150.0 };
	//GLfloat light_position[] = { -1.0, 1.0, 1.0, 0.0 };
	GLfloat light_diffuse[] = { 1.0,1.0, 1.0, 1.0 };
	//GLfloat light_ambient[] = { 0.0, 1.0, 1.0, 1.0 };
	GLfloat mat_zero[] = {0.0, 0.0, 0.0, 1.0 };
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glMaterialfv(GL_FRONT, GL_AMBIENT, gLightClr[1]);

	glMaterialfv(GL_FRONT, GL_DIFFUSE, gLightClr[0]);


	float gBLight0[4] = {1, 1,1,0};
	float gfLightPos[] = {1.0, 1.0, 1.0, 0.0};

	//glLightfv(GL_LIGHT0, GL_AMBIENT, gBLight0);  
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, gBLight0); 
	//glLightfv(GL_LIGHT0, GL_SPECULAR, gBLight0);  
	glLightfv(GL_LIGHT0, GL_POSITION, gfLightPos);

	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
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

void Normalizer(float *v) 
{
	float dis = sqrt(FindSquare(v));
	for(int i=0; i<3; ++i) v[i] /= dis;

}
 
float FindDis3D(float *p1, float*p2)
{
	return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]));
}

void CirRotationXY(float *vect1, float *cen, float r, int num, float **data)
{

	float vect[3];
	for(int i=0;i<3;++i) vect[i]=vect1[i];
	if (vect[0] < 0.0)
	{
		for(int i=0;i<3;++i) vect[i]=-vect[i];
	}

	float theta;
	if (abs(vect[0]) < 0.000001)
	{
		theta = PI/2;
	}
	else theta = atan(vect[1]/vect[0]);


	glBegin(GL_LINE_LOOP);
	for (int i=0; i<num; ++i)
	{
		float t = 2*PI*i/num; 
		float pos[3] = {cos(t)*sin(theta), cos(t)*cos(theta), sin(t)};

		if(data)
			for(int j=0;j<3;++j)
			{
				data[i][j] = cen[j]+r*pos[j];  
			}
		else
		{
			for(int j=0;j<3;++j)
			{
				pos[j] = cen[j]+r*pos[j];
				glVertex3fv(pos); 
			}

		}

	}
	glEnd();


}

void CirRotation(float *vect1, float *cen, float r, int num, float **data)
{

	float vect[3];
	for(int i=0;i<3;++i) vect[i]=vect1[i];

	int cirflag=0;
	if (vect[1]<-0.000001) 
	{
		cirflag = 1;
		for(int i=0;i<3;++i) vect[i]=-vect[i];
	}
	float q1, q2;
	q1 = acos(vect[2]/sqrt(FindSquare(vect)));
	if (abs(vect[1]) < 0.0001)
	{
		q2 = PI/2;
		if (vect[0]<-0.000001)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(vect[0]/vect[1]);

	glBegin(GL_LINE_LOOP);
	for (int i=0; i<num; ++i)
	{
		float t = 2*PI*i/num;
		//if(!cirflag)
		//	t = PI-t;
		float pos[3] = {(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) ),
			-(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t)),
			sin(q1)*sin(t)};

		if(data)
			for(int j=0;j<3;++j)
			{
				data[i][j] = cen[j]+r*pos[j];  
			}
		else
		{
			for(int j=0;j<3;++j)
			{
				pos[j] = cen[j]+r*pos[j];
				glVertex3fv(pos); 
			}

		}

	}
	glEnd();
}


void CirRotation2(float *vect1, float *cen, float r, float t, float *pos)
{

	float vect[3];
	for(int i=0;i<3;++i) vect[i]=vect1[i];

	if (vect[1]<-0.000001) for(int i=0;i<3;++i) vect[i]=-vect[i];
	float q1, q2;
	q1 = acos(vect[2]/sqrt(FindSquare(vect)));
	if (abs(vect[1]) < 0.0001)
	{
		q2 = PI/2;
		if (vect[0]<-0.000001)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(vect[0]/vect[1]);


	pos[0] = (cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
	pos[1] = 	-(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t));
	pos[2] = 	sin(q1)*sin(t);



	for(int j=0;j<3;++j)
	{
		pos[j] = cen[j]+r*pos[j]; 
	} 

}


void CirInversion(float *CirVect, float *cen, float r, float *sphPos, float SphRad,
	float *invVect, float *invCen, float *r2)
{
	float vect0[3] = {cen[0]-sphPos[0], cen[1]-sphPos[1], cen[2]-sphPos[2]}; 


	float v[3];
	CrossProduct(vect0, CirVect, v);

	 
	float v3[3]={1.0, 0, 0};
	float theta = PI/2 +acos(DotProduct(v3, v));

	float t1, t2; 
	t1 = theta+PI/2; t2 = theta-PI/2;


	for(int i=0;i<3;++i) CirVect[i]=CirVect[i];

	if (CirVect[1]<-0.000001) for(int i=0;i<3;++i) CirVect[i]=-CirVect[i];
	float q1, q2;
	q1 = acos(CirVect[2]/sqrt(FindSquare(CirVect)));
	if (abs(CirVect[1]) < 0.0001)
	{
		q2 = PI/2;
		if (CirVect[0]<-0.000001)
		{
			q2=-q2;
		}

	}
	else
		q2 = atan(CirVect[0]/CirVect[1]);


	float pos[3] ={(cos(q2)*cos(t1)-sin(q2)*cos(q1)*sin(t1) ),
		-(sin(q2)*cos(t1)+cos(q2)*cos(q1)*sin(t1)),
		sin(q1)*sin(t1)};
	for (int i=0; i<3; ++i)
	{
		pos[i] = cen[i] + r*pos[i];
	}
	
	 


	float pos2[3] ={(cos(q2)*cos(t2)-sin(q2)*cos(q1)*sin(t2) ),
		-(sin(q2)*cos(t2)+cos(q2)*cos(q1)*sin(t2)),
		sin(q1)*sin(t2)};
	for (int i=0; i<3; ++i)
	{
		pos2[i] = cen[i] + r*pos2[i];
	}


	float pos3[3] ={(cos(q2)*cos(t2+PI/2)-sin(q2)*cos(q1)*sin(t2+PI/2) ),
		-(sin(q2)*cos(t2+PI/2)+cos(q2)*cos(q1)*sin(t2+PI/2)),
		sin(q1)*sin(t2+PI/2)};
	for (int i=0; i<3; ++i)
	{
		pos3[i] = cen[i] + r*pos3[i];
	}


	float dis = FindDis3D(pos, sphPos);
	float alp = SphRad*SphRad/dis/dis;
	for(int i=0; i<3; ++i)
		pos[i] = (1-alp)*sphPos[i] + alp*pos[i];

	dis = FindDis3D(pos2, sphPos);
	alp = SphRad*SphRad/dis/dis;
	for(int i=0; i<3; ++i)
		pos2[i] = (1-alp)*sphPos[i] + alp*pos2[i];

	(*r2) = FindDis3D(pos, pos2)/2;

	glColor3f(0.1, 0.7, 0.2);
	glBegin(GL_LINES);
	glVertex3fv(pos);
	glVertex3fv(pos2);
	glEnd();

	//Finding center
	for (int i=0; i<3; ++i)
	{
		invCen[i] = (pos[i]+pos2[i])/2;
	}



	//Finding normal vector
	float v1[3] = {pos[0] - invCen[0], pos[1] - invCen[1], pos[2] - invCen[2]};
	dis = FindDis3D(pos3, sphPos);
	alp = SphRad*SphRad/dis/dis;
	for(int i=0; i<3; ++i)
		pos3[i] = (1-alp)*sphPos[i] + alp*pos3[i];

	float v2[3] = {pos3[0] - invCen[0], pos3[1] - invCen[1], pos3[2] - invCen[2]};

	CrossProduct(v1, v2, invVect);






}


void Cyclides2()
{

	// two cirles of Ruled Quadrics
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	if(gbShowSph)
		for (float theta=0; theta < 2*PI ; theta += 0.05 )
		{
			glVertex3f(cos(theta)+sin(theta), sin(theta), 0);

		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		if(gbShowSph)
			for (float theta=0; theta < 2*PI ; theta += 0.05 )
				glVertex3f(cos(theta)+ sin(theta)+ gfDis, sin(theta)+gfDis, gfDis);
		glEnd();

		// Side lines of Quadrics
		if(gbCrossLine && gbShowSph)
		{
			glBegin(GL_LINES);
			glColor3f(0.7, 0.7, 0.4);
			for (float theta=0; theta <2*PI; theta += gfAngStep)
			{
				glVertex3f(cos(theta)+sin(theta),sin(theta), 0);
				glVertex3f(cos(theta +gfAngOffSet)+ sin(theta +gfAngOffSet)+ gfDis,sin(theta +gfAngOffSet)+ gfDis ,gfDis);
			}
			//glColor3f(0.2, 0.7, 0.1);
			//for (float theta=0; theta <2*PI; theta += gfAngStep)
			//{
			//	glVertex3f(cos(theta), sin(theta), 0);
			//	glVertex3f(cos(theta -gfAngOffSet), sin(theta -gfAngOffSet), gfDis);
			//}
			glEnd();
		}

		//Sphere
		if(gbShowSph)
		{
			glPushMatrix();
			glColor3f(0.30, .0,0.2f);
			glTranslatef(gfSphPos[0], gfSphPos[1], gfSphPos[2]);
			glutSolidSphere(gfSphR/10,24, 24);
			glPopMatrix();
		}

		//Connet one point to the center of the sphere
		//cyclide data
		float fData[2][100][100][3];
		int nNums;

		//float fData2[500][100][3];
		int nCirVerNum=0;
		for (int AngFlag=0; AngFlag<2; ++AngFlag) 
		{
			int idx1=0;
			float t0 = 0;//1.5*PI+0.44;
			for (; t0 < 2*PI;  t0 += gfAngStep )
			{

				float fDraLinePair[2][3] = {cos(t0), sin(t0), 0, cos(t0+gfAngOffSet), sin(t0+gfAngOffSet), gfDis};
				if(1==AngFlag) {fDraLinePair[1][0]=cos(t0-gfAngOffSet); fDraLinePair[1][1]=sin(t0-gfAngOffSet);}

				fDraLinePair[0][1] += fDraLinePair[0][2];
				fDraLinePair[0][0] += fDraLinePair[0][1];
				fDraLinePair[1][1] += fDraLinePair[1][2];
				fDraLinePair[1][0] += fDraLinePair[1][1];






				float fAlpha = 0.23f;   //define where this point lie along the Line 0~1
				float tmp=FindDis3D(fDraLinePair[0], fDraLinePair[1]); tmp*=tmp;
				fAlpha = -((fDraLinePair[1][0]-fDraLinePair[0][0])*(fDraLinePair[0][0]-gfSphPos[0])
					+(fDraLinePair[1][1]-fDraLinePair[0][1])*(fDraLinePair[0][1]-gfSphPos[1])
					+(fDraLinePair[1][2]-fDraLinePair[0][2])*(fDraLinePair[0][2]-gfSphPos[2]))/tmp;
				float fDraLinePos[3] = { fAlpha*fDraLinePair[1][0] +(1-fAlpha)*fDraLinePair[0][0],
					fAlpha*fDraLinePair[1][1] +(1-fAlpha)*fDraLinePair[0][1],
					fAlpha*fDraLinePair[1][2] +(1-fAlpha)*fDraLinePair[0][2]};

				if(gbShowSph)
				{

					glColor3f(1.0, 0.0, .7);
					glBegin(GL_LINES);
					glVertex3fv(fDraLinePos);
					glVertex3fv(gfSphPos);
					glEnd();
				}
				float cr = (fDraLinePair[0][0]-fDraLinePair[1][0])*(fDraLinePos[0]-gfSphPos[0])
					+ (fDraLinePair[0][1]-fDraLinePair[1][1])*(fDraLinePos[1]-gfSphPos[1])
					+ (fDraLinePair[0][2]-fDraLinePair[1][2])*(fDraLinePos[2]-gfSphPos[2]);
				// cr= 0 means it is right

				////distance and vectors
				float fDisC2P= FindDis3D(fDraLinePos, gfSphPos);
				float fVectX0P[3] = {fDraLinePair[1][0]-fDraLinePair[0][0], fDraLinePair[1][1]-fDraLinePair[0][1],fDraLinePair[1][2]-fDraLinePair[0][2]};
				float fDisV = FindDis3D(fDraLinePair[0], fDraLinePair[1]);

				float fUnitV[3] =   {fVectX0P[0]/fDisV, fVectX0P[1]/fDisV, fVectX0P[2]/fDisV};
				//so:
				// theta: angel between n(z axis) and vector x0p
				float theta =   acos(fUnitV[2]);
				float fSmallR = sqrt(fDisC2P*fDisC2P-gfSphR*gfSphR/100);  

				//glColor3f(1.0, 1.0, 1.0);
				//glBegin(GL_LINES);

				//glVertex3fv(fDraLinePos);
				//glVertex3f(fDraLinePos[0]+fUnitV[0],fDraLinePos[1]+fUnitV[1],fDraLinePos[2]+fUnitV[2] );
				//glEnd();
				nCirVerNum=0;
				//  float d= FindDis3D(fDraLinePos, gfSphPos);
				float theta1= atan(fVectX0P[1]/fVectX0P[2]);
				float theta2 = asin(fVectX0P[0]/fDisV);
				float fCirVerPos[100][3];
				for (float t=0.0; t < 2*PI; t += PI/50)
				{
					float fV[3] = {cos(t), sin(t), 0};
					fV[2] = -(cos(t)*fVectX0P[0]+sin(t)*fVectX0P[1])/fVectX0P[2];
					float dis = sqrt(1+fV[2]*fV[2]); 



					//float k3=(fUnitV[0]*sin(t)-fUnitV[1]*cos(t))*(1-cos(theta));
					//fV[0] = cos(t)*cos(theta)  - k3*fUnitV[1];
					//fV[1] = sin(t)*cos(theta) + k3*fUnitV[0];
					//fV[2] = -(fUnitV[1]*sin(t)+fUnitV[0]*cos(t))*sin(theta);


					//float origin[3]={cos(t), sin(t), 0};

					//fV[0] = cos(t)*cos(theta2)+sin(theta2)*(sin(t)*sin(theta1));
					//fV[1] = sin(t)*cos(theta1);
					//fV[2] = -cos(t)*sin(theta2)+cos(theta2)*sin(t)*sin(theta1);

					//fV[0] = cos(t); 
					//fV[1] = sin(t);
					//fV[2] = 0;
					//1: scale

					for(int i=0;i<3;++i) fV[i]*=fSmallR/dis;

					//2: translate
					for(int i=0;i<3;++i) fV[i]+=fDraLinePos[i];


					for(int i=0;i<3;++i) fData[AngFlag][idx1][nCirVerNum][i]=fV[i];
					++nCirVerNum;
				}

				//glPushMatrix();
				//glTranslatef(fDraLinePos[0], fDraLinePos[1], fDraLinePos[2]);
				//glutSolidSphere(fSmallR,20,20);
				//glPopMatrix();
				++idx1;
			}
			nNums = idx1;
			//glColor3f(1.0,0.0,0.4);
			//if(gbCrossCyc) //cross lines
			//	for (int j=0; j<nCirVerNum; ++j)
			//	{
			//		glBegin(GL_LINE_LOOP);
			//		for (int i=0; i<idx1;++i)
			//		{
			//			glVertex3fv(fData[AngFlag][i][j]);
			//		}
			//		glEnd();
			//	}
				if(1==AngFlag) 
					glColor3f(0.1,0.70,8.1);
				else
					glColor3f(0.9,0.70,0.17);

				if(gbCirCyc && gbFlags[0] && 0==AngFlag ||gbCirCyc && gbFlags[1] && 1==AngFlag )
					for (int i=0; i<idx1;++i)
					{
						glBegin(GL_LINE_LOOP);
						for (int j=0; j<nCirVerNum; ++j)
						{
							//glVertex3fv(fDraLinePos);
							glVertex3fv(fData[AngFlag][i][j]);
						}
						glEnd();

					}

					int bShowSide = 0;
					if(bShowSide && gbCirCyc)
					{


						glEnable(GL_LIGHTING);
						glEnable(GL_LIGHT0);

						glBegin(GL_TRIANGLES);


						if(1==nFlag)
						{
							ofstream outfile2("Circles.inc");
							for (int i=0; i<idx1; i+=2)
							{
								for (int j=0; j<nCirVerNum; ++j)
								{ 

									outfile2 << "cylinder { <"
										<< fData[1][i][j][0] <<", " 
										<< fData[1][i][j][1] <<", "
										<< -fData[1][i][j][2] <<">, <" 
										<< fData[1][i][(j+1)%nCirVerNum][0] <<", " 
										<< fData[1][i][(j+1)%nCirVerNum][1] <<", "
										<< -fData[1][i][(j+1)%nCirVerNum][2] <<">, r_edge texture { edge_tex } }"
										<<endl;

								}

							}





							ofstream outfile("Cyclide1.obj"); 
							outfile <<"# A cyclide generated by a ruled quadric and a sphere" <<endl;
							outfile <<"# number of vertex: " << idx1*nCirVerNum <<endl;
							outfile <<"# number of faces: "<< (idx1-1)*nCirVerNum*2<<endl;


							for (int i=0; i<idx1; ++i)
							{
								for (int j=0; j<nCirVerNum; ++j)
								{ 
									outfile <<"v " <<fData[AngFlag][i][j][0]<<" "
										<<fData[AngFlag][i][j][1]<<" "
										<<fData[AngFlag][i][j][2]<<endl;
								}
							}

							for (int i=0; i<idx1; ++i)
							{
								for (int j=0; j<nCirVerNum; ++j)
								{ 
									outfile <<"f "<< i*nCirVerNum+j+1 <<" "
										<<i*nCirVerNum+(j+1)%nCirVerNum+1<<" "
										<<((i+1)%idx1) * nCirVerNum +j+1<<endl;



									outfile <<"f "<<i*nCirVerNum+(j+1)%nCirVerNum+1<<" "
										<<((i+1)%idx1) * nCirVerNum+(j+1)%nCirVerNum+1<<" "
										<<((i+1)%idx1) * nCirVerNum +j+1<<endl;

								}
							}

						}
						for (int i=0; i<idx1; ++i)
						{
							for (int j=0; j<nCirVerNum; ++j)
							{
								//i+1 = (i+1)%idx1
								float v2[3]={fData[AngFlag][(i+1)%idx1][j][0]-fData[AngFlag][i][j][0], fData[AngFlag][(i+1)%idx1][j][1]-fData[AngFlag][i][j][1], fData[AngFlag][(i+1)%idx1][j][2]-fData[AngFlag][i][j][2]  };
								float v1[3]={fData[AngFlag][i][(j+1)%nCirVerNum][0]-fData[AngFlag][i][j][0], fData[AngFlag][i][(j+1)%nCirVerNum][1]-fData[AngFlag][i][j][1], fData[AngFlag][i][(j+1)%nCirVerNum][2]-fData[AngFlag][i][j][2]  };

								float fNormal[3];
								fNormal[0] = v1[1]*v2[2]-v1[2]*v2[1];
								fNormal[1] = v1[2]*v2[0]-v1[0]*v2[2];
								fNormal[2] = v1[0]*v2[1]-v1[1]*v2[0];
								float N=sqrt(fNormal[0]*fNormal[0]+fNormal[1]*fNormal[1]+fNormal[2]*fNormal[2]);


								if (0==N)
								{
									fNormal[0]=fNormal[1]=0; fNormal[2]=1;
								}
								else
									for(int j=0; j<3;++j) fNormal[j] /= N;
								glNormal3fv(fNormal);
								glVertex3fv(fData[AngFlag][i][j]);
								glVertex3fv(fData[AngFlag][i][(j+1)%nCirVerNum]);
								glVertex3fv(fData[AngFlag][(i+1)%idx1][j]); 

								v2[0] = fData[AngFlag][(i+1)%idx1][(j+1)%nCirVerNum][0]-fData[AngFlag][i][(j+1)%nCirVerNum][0];
								v2[1] = fData[AngFlag][(i+1)%idx1][(j+1)%nCirVerNum][1]-fData[AngFlag][i][(j+1)%nCirVerNum][1];
								v2[2] = fData[AngFlag][(i+1)%idx1][(j+1)%nCirVerNum][2]-fData[AngFlag][i][(j+1)%nCirVerNum][2];
								v1[0] = fData[AngFlag][(i+1)%idx1][j][0]-fData[AngFlag][i][(j+1)%nCirVerNum][0];
								v1[1] = fData[AngFlag][(i+1)%idx1][j][1]-fData[AngFlag][i][(j+1)%nCirVerNum][1];
								v1[2] = fData[AngFlag][(i+1)%idx1][j][2]-fData[AngFlag][i][(j+1)%nCirVerNum][2];


								fNormal[0] = v1[1]*v2[2]-v1[2]*v2[1];
								fNormal[1] = v1[2]*v2[0]-v1[0]*v2[2];
								fNormal[2] = v1[0]*v2[1]-v1[1]*v2[0];
								N=sqrt(fNormal[0]*fNormal[0]+fNormal[1]*fNormal[1]+fNormal[2]*fNormal[2]);


								if (0==N)
								{
									fNormal[0]=fNormal[1]=0; fNormal[2]=1;
								}
								else
									for(int j=0; j<3;++j) fNormal[j] /=N;
								glNormal3fv(fNormal);
								glVertex3fv(fData[AngFlag][i][(j+1)%nCirVerNum]);
								glVertex3fv(fData[AngFlag][(i+1)%idx1][(j+1)%nCirVerNum]);
								glVertex3fv(fData[AngFlag][(i+1)%idx1][j]);
							}
						}
						glEnd();

						glDisable(GL_LIGHT0); 
						glDisable(GL_LIGHTING);
					}



		}

		// Inversion and something else
		float fPosSp2[3] = {-7.0f, 0.0f, 0.0f};

		//float gfSphR = 1.7f;
		//float gfSphPos[3] = { -2.5,  -3.7, 0.7};
		fPosSp2[2] = gfSphPos[2]-gfSphR/10;
		fPosSp2[1] = gfSphPos[1]; fPosSp2[0] = gfSphPos[0];
		glPushMatrix();
		glColor3f(0.4f, 0.4f, 0.4f);
		glTranslatef(fPosSp2[0], fPosSp2[1], fPosSp2[2]);

		if(gbShowSp2)
			glutSolidSphere(gfr2, 20, 20);
		//glutWireSphere(r2, 20, 20);

		glPopMatrix();

		//Sphere
		//float gfSphR = 0.36f;
		//float gfSphPos[3] = { -2.5,  -3.7, 0.7};
		if(gbShowSp2)
			for (int i=0; i< nNums;++i)
			{
				glColor3f(0.2f, 0.1f, 0.9f);  // blue lines, for
				glBegin(GL_LINE_LOOP);
				for (int j=0; j<nCirVerNum; ++j)
				{ 
					//float fPos[3]= { gfSphR*cos(u)*cos(v)+gfSphPos[0], gfSphR*cos(u)*sin(v)+gfSphPos[1], gfSphR*sin(u)+gfSphPos[2]};
					float *fPos = fData[0][i][j];
					float dis = FindDis3D(fPos, fPosSp2);

					float alp = gfr2*gfr2/dis/dis; 
					glVertex3f((1-alp)*fPosSp2[0] + alp * fPos[0] , (1-alp)*fPosSp2[1] + alp * fPos[1] ,(1-alp)*fPosSp2[2] + alp * fPos[2]  );

				}
				glEnd();


				glColor3f(1.0f, 1.0f, 1.0f); //white lines
				glBegin(GL_LINE_LOOP);
				float u = (1.0*i/nNums - 0.5)*PI;
				for (int j=0; j<2*nCirVerNum; ++j)
				{
					float v=2.0*PI*j/2/nCirVerNum;
					float fPos[3];
					fPos[0] = gfSphR*cos(u)*cos(v)+gfSphPos[0];
					fPos[1] = gfSphR*cos(u)*sin(v)+gfSphPos[1];
					fPos[2] = gfSphR*sin(u)+gfSphPos[2];
					float dis = FindDis3D(fPos, fPosSp2);
					float alp = gfr2*gfr2/dis/dis;
					glVertex3f((1-alp)*fPosSp2[0] + alp * fPos[0] , (1-alp)*fPosSp2[1] + alp * fPos[1] ,(1-alp)*fPosSp2[2] + alp * fPos[2]  );

				}

				glEnd();
				//break;
			}




}


void Cyclides()
{

	// two cirles of Ruled Quadrics
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (float theta=0; theta < 2*PI ; theta += 0.05 )
	{
		glVertex3f(cos(theta), sin(theta), 0);

	}
	glEnd();
	glBegin(GL_LINE_STRIP);
	for (float theta=0; theta < 2*PI ; theta += 0.05 )
		glVertex3f(cos(theta), sin(theta), gfDis);
	glEnd();

	// Side lines of Quadrics
	glColor3f(0.1, 0.7, 0.4);
	glBegin(GL_LINES);
	for (float theta=1.5*PI; theta <2*PI; theta += 0.11)
	{
		glVertex3f(cos(theta), sin(theta), 0);
		glVertex3f(cos(theta - PI/6), sin(theta-PI/6), gfDis);
	}
	glEnd();

	//Sphere
	glPushMatrix();
	glColor3f(1.0, 1.0,1.0f);
	glTranslatef(gfSphPos[0], gfSphPos[1], gfSphPos[2]);
	glutSolidSphere(gfSphR/10,24, 24);
	glPopMatrix();

	//Connet one point to the center of the sphere
	//float t0 = 1.5*PI+0.44;
	float fDraLinePair[2][3] = {cos(t0), sin(t0), 0, cos(t0-PI/6), sin(t0-PI/6), gfDis};
	//float fAlpha = 0.23f;   //define where this point lie along the Line 0~1
	float fDraLinePos[3] = { fAlpha*fDraLinePair[1][0] +(1-fAlpha)*fDraLinePair[0][0],
		fAlpha*fDraLinePair[1][1] +(1-fAlpha)*fDraLinePair[0][1],
		fAlpha*fDraLinePair[1][2] +(1-fAlpha)*fDraLinePair[0][2]};

	glColor3f(0.0, 0.0, .7);
	glBegin(GL_LINES);
	glVertex3fv(fDraLinePos);
	glVertex3fv(gfSphPos);
	glEnd();

	//distance and vectors
	float fDisC2P= FindDis3D(fDraLinePos, gfSphPos);
	float fVectX0P[3] = {gfSphPos[0]-fDraLinePos[0], gfSphPos[1]-fDraLinePos[1],gfSphPos[2]-fDraLinePos[2]};
	float fUnitV[3] = {fVectX0P[0]/fDisC2P, fVectX0P[1]/fDisC2P, fVectX0P[2]/fDisC2P};
	//so:
	float fVectV[3] = {-fUnitV[2], fUnitV[1], 0}; 
	// theta: angel between n(z axis) and vector x0p
	float theta =  acos(fUnitV[2]);
	float fSmallR = sqrt(fDisC2P*fDisC2P-gfSphR*gfSphR/100);//新球的半径
	float fAlpha0 = fSmallR/fDisC2P; fAlpha0 *= fAlpha0;
	float fCPos[3] = {(1-fAlpha0)*fDraLinePos[0]+fAlpha0*gfSphPos[0],
		(1-fAlpha0)*fDraLinePos[1]+fAlpha0*gfSphPos[1],
		(1-fAlpha0)*fDraLinePos[2]+fAlpha0*gfSphPos[2]};

	float fCRadius = fSmallR*gfSphR/10/fDisC2P;
	glColor3f(1.0,0.0,0.4);
	glBegin(GL_LINE_STRIP);
	int nCirVerNum=0;
	float fCirVerPos[100][3];
	for (float t=0.0; t < 2*PI; t += 0.06)
	{
		float fV[3];
		float k3=(fUnitV[0]*sin(t)-fUnitV[1]*cos(t))*(1-cos(theta));
		fV[0] = cos(t)*cos(theta)  - k3*fUnitV[1];
		fV[1] = sin(t)*cos(theta) + k3*fUnitV[0];
		fV[2] = -(fUnitV[1]*sin(t)+fUnitV[0]*cos(t))*sin(theta);
		//1: scale
		for(int i=0;i<3;++i) fV[i]*=fCRadius;

		//2: translate
		for(int i=0;i<3;++i) fV[i]+=fCPos[i];

		glVertex3fv(fV); 
		for(int i=0;i<3;++i) fCirVerPos[nCirVerNum][i]=fV[i];
		++nCirVerNum;
	}
	glEnd();
	glColor3f(0.4,0.4,0.4);
	glBegin(GL_LINES);
	for (int i=0; i<nCirVerNum; ++i)
	{
		glVertex3fv(fDraLinePos);
		glVertex3fv(fCirVerPos[i]);

	}
	glEnd();

}
void InvSphInversion(float *fCirPos, float fRad, float theta, float *fSphPos, float fSphRad)
{

}


ofstream outfile2("lines15.inc");
ofstream outfile4("circles14.inc");
void SphInversion(float *fCirPos, float fRad, float theta, float *fSphPos, float fSphRad, float *fAxis = NULL)
{ 
	//glColor3f(0.5, 0.8, 0.2);
	glColor3f(1.0f, 0.0, 0.0);
	// circle 1
	glBegin(GL_LINE_LOOP);
	if(gbFlags[0]) //show circle
		for (float t=0.0; t < 2*PI; t += 0.06)
		{
			float Pos[3];
			Pos[0] = fRad*cos(t)*cos(theta) + fCirPos[0];
			Pos[1] = fRad*cos(t)*sin(theta) + fCirPos[1];
			Pos[2] = fRad*sin(t) + fCirPos[2];
			glVertex3fv(Pos);

			float Pos2[3];
			float t2 = t+0.06f;
			Pos2[0] = fRad*cos(t2)*cos(theta) + fCirPos[0];
			Pos2[1] = fRad*cos(t2)*sin(theta) + fCirPos[1];
			Pos2[2] = fRad*sin(t2) + fCirPos[2];
			if(nFlag)
			{ 
			OutF(&outfile4, Pos, Pos2);
			}

		}
		glEnd();

		if(1)
		{
			if (nFlag  )
			{


				//outfile4<<"#declare edge_tex=  my_color_texture_1; " ;
			}
			float t=0;
			float Pos[3];
			Pos[0] = fRad*cos(t)*cos(theta) + fCirPos[0];
			Pos[1] = fRad*cos(t)*sin(theta) + fCirPos[1];
			Pos[2] = fRad*sin(t) + fCirPos[2];
			float Pos2[3];
			float t2 = PI;
			Pos2[0] = fRad*cos(t2)*cos(theta) + fCirPos[0];
			Pos2[1] = fRad*cos(t2)*sin(theta) + fCirPos[1];
			Pos2[2] = fRad*sin(t2) + fCirPos[2];
		//	if(nFlag)
			//	OutF(&outfile4, Pos, Pos2);
		}

		//find the shortest and longest distances from the circle to the sphere center
		float A, B;
		A = (fCirPos[0]-fSphPos[0])*cos(theta) + (fCirPos[1]-fSphPos[1])*sin(theta);
		B = fCirPos[2] - fSphPos[2];
		glBegin(GL_LINES);
		if(gbFlags[0])
			for (float t= PI/2+asin(A/sqrt(A*A+B*B)); t < PI*2  ; t += PI)
			{  
				float Pos[3];
				Pos[0] = fRad*cos(t)*cos(theta) + fCirPos[0];
				Pos[1] = fRad*cos(t)*sin(theta) + fCirPos[1];
				Pos[2] = fRad*sin(t) + fCirPos[2];
				//glVertex3fv(Pos);
			}
			glEnd();

			//Inversion
			glColor3f(0.5, 0.2, 0.8);
			glBegin(GL_LINE_LOOP);


			float Pos2[3];
			//Inversion of the circle
			for (float t=0.0; t < PI*2  ; t += 0.06)
			{
				float Pos[3];
				Pos[0] =fRad* cos(t)*cos(theta) + fCirPos[0];
				Pos[1] =fRad* cos(t)*sin(theta) + fCirPos[1];
				Pos[2] =fRad* sin(t) + fCirPos[2];
				float dis = FindDis3D(Pos, fSphPos);

				float alp = fSphRad*fSphRad/dis/dis; 
				Pos2[0] =  (1-alp)*fSphPos[0] + alp * Pos[0] ;
				Pos2[1] =  (1-alp)*fSphPos[1] + alp * Pos[1] ;
				Pos2[2] =  (1-alp)*fSphPos[2] + alp * Pos[2] ;
				if(gbFlags[1]) //show inv circle
					glVertex3fv( Pos2 );

			}
			glEnd();


			glBegin(GL_LINES);
			float fChord[2][3]; int i=0;

			for (float t=   PI/2+asin(A/sqrt(A*A+B*B)); t < PI*2  ; t += PI)
			{

				float Pos[3];
				Pos[0] = fRad*cos(t)*cos(theta) + fCirPos[0];
				Pos[1] = fRad*cos(t)*sin(theta) + fCirPos[1];
				Pos[2] = fRad*sin(t) + fCirPos[2];
				float dis = FindDis3D(Pos, fSphPos);

				float alp = fSphRad*fSphRad/dis/dis; 
				float Pos2[3] = {(1-alp)*fSphPos[0] + alp * Pos[0] , (1-alp)*fSphPos[1] + alp * Pos[1] ,(1-alp)*fSphPos[2] + alp * Pos[2] };
				if(gbFlags[1])
					glVertex3fv( Pos2 );
				//glVertex3f(0.0f, 0.0f, 0.0f);
				memcpy(fChord[i++], Pos2, 3*sizeof(float));

			}
			float fCenter[3] = {(fChord[0][0]+fChord[1][0])/2,
				(fChord[0][1]+fChord[1][1])/2,
				(fChord[0][2]+fChord[1][2])/2};


			//glVertex3fv(fChord[0]);
			for(int i=0;i<3;++i)
			{

				fChord[0][i]=fChord[0][i]-fCenter[i];
				fChord[1][i]=Pos2[i]-fCenter[i];
			}
			//glVertex3fv(fCenter); 
			//glVertex3fv(Pos2); 
			//glVertex3fv(fCenter); 
			CrossProduct(fChord[0], fChord[1], Pos2);

			//float dis=sqrt(FindSquare(Pos2)); for(int i=0; i<3; ++i) Pos2[i]/=dis; //normalization
			fAxis[0]=fCenter[0]-Pos2[0]; fAxis[1]=fCenter[1]-Pos2[1];fAxis[2]=fCenter[2]-Pos2[2];
			fAxis[3]=Pos2[0]+fCenter[0]; fAxis[4]=Pos2[1]+fCenter[1];fAxis[5]=Pos2[2]+fCenter[2];
			//glVertex3fv(fAxis);
			//glVertex3fv(fAxis+3);

			glEnd();


}


ofstream  outx("circles15.inc");
void CircleBy2Point(float *p1, float *p2, float *SphPos, float *pR, 
	float *fCirPos = NULL, float *fTheta=NULL, float *radius=NULL)
{
	float v[3];
	int num = 50;

	float qq2[3]={p2[0], p2[1], -p2[2]};
	float a, b,c;
	a=FindDis3D(p1, p2); b=FindDis3D(p2, qq2); c=FindDis3D(p1, qq2);
	float r=a*b*c/b/2/sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));

	a = sqrt(r*r-p1[2]*p1[2]);
	b = sqrt(r*r-p2[2]*p2[2]);
	float cen[3] ={ (b*p1[0]+a*p2[0])/(a+b), (b*p1[1]+a*qq2[1])/(a+b), 0};


	float pp1[3];
	float pp2[3];
	for (int i=0; i<3; ++i)
	{
		pp1[i] = p1[i]-cen[i];
		pp2[i] = p2[i]-cen[i];
	}
	CrossProduct(pp1, pp2, v);

	float q1, q2;
	q1 = acos(v[2]);
	if (abs(v[1]) < 0.0001)
		q2 = PI/2;
	else
		q2 = atan(v[0]/v[1]);

	float q0 = atan(v[1]/v[0]);
	if (NULL != fTheta)
	{
		(*fTheta) = q0 + PI/2;
	}


	glBegin(GL_LINE_LOOP);
	glColor3f(192.0/255, 192.0/255, 192.0/255);
	for (int i=0; i<num; ++i)
	{
		float t = 2.0*PI*i/num;
		float Pos[3];
		Pos[0] = cen[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
		Pos[1] = cen[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t) );
		Pos[2] = cen[2]+ r*sin(q1)*sin(t); 
		 
		glVertex3fv(Pos);
		if(nFlag)
		outfile2 <<"cylinder { <" << Pos[0] <<", " << Pos[1]<<", "<< Pos[2] <<">, <";

		t=2.0*PI*(i+1)/num;

		Pos[0] = cen[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
		Pos[1] = cen[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t) );
		Pos[2] = cen[2]+ r*sin(q1)*sin(t); 
		if(nFlag)
		outfile2<< Pos[0] <<", " << Pos[1]<<", "<< Pos[2] <<">, r_edge texture { my_color_texture_0 } } "<<endl;



	}

	if (NULL != radius)
	{
		(*radius) = r;
	}
	if (NULL != fCirPos)
	{
		fCirPos[0] = cen[0]; fCirPos[1] = cen[1]; fCirPos[2] = cen[2];
	}
	glEnd();


	//////////////////////////////////////////////////////////////////////////
	//calculating sphere by this circle and a further point

	float Px[3] = {P[0], P[1], 0.0f}; //projection to the xy plane
	
	glColor3f(0.5f, 1.0f, 170.0/255);
	glBegin(GL_LINES);
	glVertex3fv(P);
	glVertex3fv(Px);
	glEnd();



	float D = FindDis3D(cen, Px); 
	a = (D*D + P[2]*P[2]  - r*r)/2/D;
	a = (r*r-(cen[0]-P[0])*(cen[0]-P[0]) -(cen[1]-P[1])*(cen[1]-P[1]) - P[2]*P[2] );
	a /=  2*((cen[0]-P[0])*v[0] + (cen[1]-P[1])*v[1]);

	float R = sqrt(r*r+a*a); //radius of sphere
	(*pR) = R;
	//v[0] = cen[0]-Px[0]; v[1] = cen[1]-Px[1];
	//float dis=sqrt(v[0]*v[0]+v[1]*v[1]);
	//v[0]/=dis; v[1]/=dis;
	float SpCen[3] = {cen[0]+a*v[0], cen[1]+a*v[1], 0.0f};

	float d = FindDis3D(SpCen, P);
	float d2 = FindDis3D(SpCen,p1);
	float d3 = FindDis3D(SpCen, p2);

	memcpy(SphPos, SpCen, 3*sizeof(float)); 


}





void CircleInversion(float *cen, float r, float *CirVect, float *fCen2, float *R2, float *vect2)
{


	float fSphPos[3] = {0.0f, 0.0f,  2.0f};
	float fSphRad = 1.5f;
	//float theta = atan((fSphPos[1]-cen[1])/(fSphPos[0]-cen[0]));


	// the vector of the circle center and the sphere center
	float vect0[3] = {cen[0]-fSphPos[0], cen[1]-fSphPos[1], cen[2]-fSphPos[2]}; 

	float x1, y1, z1; //projection point
	x1 = CirVect[0] * (DotProduct(cen, CirVect)-2.0f*CirVect[2]);
	y1 = x1*CirVect[1]/CirVect[0];
	z1 = x1*CirVect[2]/CirVect[0] +2.0f;
	float v0[3] = {x1-cen[0], y1-cen[1], z1-cen[2]};
	Normalizer(v0);

	//float postion[2][3]; 

	float theta = atan(vect0[1]/vect0[0]);

	if (theta < 0.0)
	{
		theta += PI;
	}



	// find the paired points having longest and shortest distance with the circle center.
	float pos[2][3];
	int idx =0;
	for (float t = theta; t < 2*PI; t+=PI)
	{
		CirRotation2(CirVect, cen, r, t, pos[idx]);

		pos[idx][0] = cos(t);
		pos[idx][1] = sin(t);

		pos[idx][2] = -(cos(t)*CirVect[0]+sin(t)*CirVect[1])/CirVect[2];

		float dis1 = sqrt(1+pos[idx][2]*pos[idx][2]);  

		for(int i=0;i<3;++i)
		{
			pos[idx][i]*= r/dis1 ;
			pos[idx][i] +=cen[i];
		} 
		++idx;
	}

	for (int i=0; i<3; ++i)
	{
		pos[0][i] = cen[i] + r*v0[i];
		pos[1][i] = cen[i] - r*v0[i];

	} 
	//glColor3f(0.99, 0.1, 0.1);
	//glBegin(GL_LINES);
	//glVertex3fv(pos[0]);
	//glVertex3fv(pos[1]);
	//glEnd();



	// a third point
	float pos2[3];
	float t = theta +PI/2;
	pos2[0] = cos(t);
	pos2[1] = sin(t);
	pos2[2] = -(cos(t)*CirVect[0]+sin(t)*CirVect[1])/CirVect[2];
	float dis1 = sqrt(1+pos2[2]*pos2[2]);   
	for(int i=0;i<3;++i)
	{
		pos2[i]*= r/dis1 ;
		pos2[i] +=cen[i];
	}


	//inversion
	for (int i=0; i<2; ++i)
	{
		float dis = FindDis3D(pos[i], fSphPos); 
		float alp = fSphRad*fSphRad/dis/dis;  
		pos[i][0] =  (1-alp)*fSphPos[0] + alp * pos[i][0] ;
		pos[i][1] =  (1-alp)*fSphPos[1] + alp * pos[i][1] ;
		pos[i][2] =  (1-alp)*fSphPos[2] + alp * pos[i][2] ; 
	}
	//glVertex3fv(pos[0]); 
	//glVertex3fv(pos[1]); 
	float dis = FindDis3D(pos2, fSphPos); 
	float alp = fSphRad*fSphRad/dis/dis;  
	pos2[0] =  (1-alp)*fSphPos[0] + alp * pos2[0] ;
	pos2[1] =  (1-alp)*fSphPos[1] + alp * pos2[1] ;
	pos2[2] =  (1-alp)*fSphPos[2] + alp * pos2[2] ; 
	 

	for (int i=0; i<3; ++i)
	{
		fCen2[i] = (pos[0][i]+pos[1][i])/2;
	}

	(*R2) = FindDis3D(pos[0], pos[1])/2;

	float v1[3], v2[3];
	for (int i=0; i<3; ++i)
	{
		v1[i] = pos[0][i]-fCen2[i];
		v2[i] = pos2[i]-fCen2[i];
	}

	CrossProduct(v1, v2, vect2);





}


void InterCircle2Spheres(float m1[], float r1, float m2[], float r2,
	float *fCirPos = NULL, float *fTheta=NULL, float *radius=NULL)
{
	float D =  FindDis3D(m1, m2 );
	float cosTheta=(r1*r1+D*D-r2*r2)/r1/D/2;
	float D1=r1*cosTheta;
	float r=sqrt(r1*r1-D1*D1);//radius of intersected circle of the two spheres
	float alph = D1/D;
	float m0[3] = {alph * m2[0] + (1-alph)*m1[0], alph * m2[1] + (1-alph)*m1[1],alph * m2[2] + (1-alph)*m1[2]}; //center of the circle


	int idx = 0; 
	float Vect[3]={m2[0]-m1[0], m2[1]-m1[1], m2[2]-m1[2]};
	float q1, q2;
	q1 = acos(Vect[2]/sqrt(FindSquare(Vect)));
	if (abs(Vect[1]) < 0.0001)
		q2 = PI/2;
	else
		q2 = atan(Vect[0]/Vect[1]);

	float q0 = atan(Vect[1]/Vect[0]);
	if(NULL != fTheta)
		(*fTheta) = q0 +PI/2;
	//-1.5338290 + PI/2

	// -0.0012483435 0.19625540 0.0000

	glColor3f(0.9, 0.2, 0.1);

	glBegin(GL_LINE_LOOP);

	float fAngStep=0.032f;
	for (float t=0.0; t < 2*PI; t += fAngStep)
	{
		float fV[3] = {cos(t), sin(t), 0};

		//}
		fV[0] = m0[0]+ r*(cos(q2)*cos(t)-sin(q2)*cos(q1)*sin(t) );
		fV[1] = m0[1]- r*(sin(q2)*cos(t)+cos(q2)*cos(q1)*sin(t) );
		fV[2] = m0[2]+ r*sin(q1)*sin(t);  

		glVertex3fv(fV);

	}
	glEnd();
	if (NULL != fCirPos)
	{
		fCirPos[0] = m0[0]; fCirPos[1] = m0[1]; fCirPos[2] = m0[2];
	}

	if (NULL != radius)
	{
		(*radius) = r;
	}
}

void RuledSurfaces2()
{
	glPushMatrix();
	glScalef(3.0f, 3.0f, 3.0f);

	int num=70;

	//plane: - xy


	if (gbFlags[7])
	{
		glColor3f(1.0f, 1.0f, 1.0f);
		glBegin(GL_QUADS);
		glVertex3f(2.0f, 2.0f, 0.0f);
		glVertex3f(-2.0f, 2.0f, 0.0f);
		glVertex3f(-2.0f, -2.0f, 0.0f);
		glVertex3f(2.0f,  -2.0f, 0.0f);
		glEnd(); 
	}

	//float P1[3] = {-cos(PI/3), 1.0, sin(PI/3)};
	//float P2[3] = {cos(PI/4), 1.0, sin(PI/4)};
	//float P3[3] = {0.9, -0.3, 1.2};
	//float P4[3] = {-1.1,-0.8, 0.97};



	float Q1[3] = {P1[0], P1[1], -P1[2]};
	float Q2[3] = {P2[0], P2[1], -P2[2]};
	float Q3[3] = {P3[0], P3[1], -P3[2]};
	float Q4[3] = {P4[0], P4[1], -P4[2]};
	
	glColor3f(0.0/255, 255.0/255, 127.0/255);
	glBegin(GL_LINE_LOOP);
	glVertex3fv(P1);
	glVertex3fv(P2);
	glVertex3fv(P3);
	glVertex3fv(P4);
	glEnd();
	glBegin(GL_LINE_LOOP);
	glVertex3fv(Q1);
	glVertex3fv(Q2);
	glVertex3fv(Q3);
	glVertex3fv(Q4);
	glEnd();

	
	float SpPos1[3], R1;
	float SpPos2[3], R2;

	CircleBy2Point(P2, P3, SpPos1, &R1, gfCirPos1, &gfTheta[0], &gfRad[0]);
	CircleBy2Point(P4, P1, SpPos1, &R1, gfCirPos3, &gfTheta[2], &gfRad[2]);


	CircleBy2Point(P1, P2, SpPos1, &R1); 

	CircleBy2Point(P3, P4, SpPos2, &R2);


	//glPushMatrix(); 
	//glColor3f(30.0/255, 144.0/255, 1.0f);
	//glTranslatef(SpPos2[0], SpPos2[1], 0.0);
	//glutWireSphere(R2, 30, 30);
	//glPopMatrix();

	InterCircle2Spheres(SpPos1, R1, SpPos2, R2, gfCirPos2, &gfTheta[1], &gfRad[1]);


	float fSphPos[3] = {0.0f, 0.0f,  2.0f};
	fSphPos[2] = gfSphdata[0];
	float fSphRad = gfSphdata[1];


	float fVectors[6][3] = 
	{ 0, 0, 0, 2, 2, -0.1,
	0, 0, 1, 3, 1, 1,
	-1, 0, 1.9, 1, -0.8, 2};


	glPushMatrix();
	glTranslatef(fSphPos[0], fSphPos[1], fSphPos[2]);
	//glutSolidSphere(fSphRad, 20, 20);

	//    glutWireSphere(fSphRad, 10, 10);
	glPopMatrix();
	//Inversion of plane
	float fSphCenter = fSphRad*fSphRad/sqrt(FindSquare(fSphPos));
	float fSphRad2 = fSphCenter/2;
	fSphCenter = fSphPos[2] - fSphRad2;
	float fSphPos2[3]={0, 0, fSphCenter};
	glPushMatrix();
	glTranslatef(fSphPos2[0], fSphPos2[1], fSphPos2[2]); 

	if(gbShowAxis)    glutSolidSphere(fSphRad2, 20, 20);
	glPopMatrix();

	// circle 1
	//float fCirPos1[3] = {1.0f, 0.0f, 0.0f};
	//float theta=PI/3;
	//float rad = 0.6;
	float fAxis[18];
	SphInversion(gfCirPos1, gfRad[0], gfTheta[0], fSphPos, fSphRad, fAxis);
	//memcpy(fVectors , fAxis, 18*sizeof(float));

	//circle 2
	//float fCirPos2[3] = {-0.11, 0.6f, 0.0f};
	//theta= PI/5;
	//rad = 0.6f;
	SphInversion(gfCirPos2, gfRad[1], gfTheta[1], fSphPos, fSphRad, fAxis+6);

	//circle 3
	//float fCirPos3[3] = { -0.71, 0.8f, 0};
	//theta= PI/3;
	//rad = 0.6;
	SphInversion(gfCirPos3, gfRad[2], gfTheta[2], fSphPos, fSphRad, fAxis+12); 
	for(int i=0; i<18; ++i)		fVectors[i/3][i%3]=fAxis[i];




	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_LINES);
	if(gbFlags[2]) //show axis
		for (int i=0; i<6; ++i)
		{
			glVertex3fv(fVectors[i]);
		}
		glEnd();

		float fCroVectors[6][3];
		glColor3f(0.1f, 0.1f, 0.7f); // rulings, group 1
		
		
		
		//////////////////////////////////////////////////////////////////////////
		// for saving obj file...
		//////////////////////////////////////////////////////////////////////////
		
		int nCircles=0;
		num=0;
		ofstream outfile("cyclide15.obj");
		ofstream out15A("lines15A.inc"); 
		ofstream out15B("lines15B.inc"); 


		int tmp = 0; 
		float data2[20] = { -20, -3, -1.2,-0.4, 0.0, 0.2, 0.9, 1.1, 1.4,
			1.9, 2.7, 5.4,  0.932, 0.945, 0.958, 0.975, 0.99, 1.02,   0.915, 0.922};
		//for (int kk=0; kk <20; ++kk)  
		for (float h = 1.1; h>-0.9; h-= 0.0015)
		{ 
			  

			float k = h;//  data2[kk];
			 
		  
			 
			float Point[3] = {k*fVectors[0][0]+(1-k)*fVectors[1][0], k*fVectors[0][1]+(1-k)*fVectors[1][1], k*fVectors[0][2]+(1-k)*fVectors[1][2]};
			float fNorms[2][3];
			float fVect[2][3] = {Point[0]-fVectors[2][0], Point[1]-fVectors[2][1], Point[2]-fVectors[2][2],
				Point[0]-fVectors[3][0], Point[1]-fVectors[3][1], Point[2]-fVectors[3][2]};
			CrossProduct(fVect[0], fVect[1], fNorms[0]); // to get the normal of plane 1
			float fVect2[2][3] = {Point[0]-fVectors[4][0], Point[1]-fVectors[4][1], Point[2]-fVectors[4][2],
				Point[0]-fVectors[5][0], Point[1]-fVectors[5][1], Point[2]-fVectors[5][2]};
			CrossProduct(fVect2[0], fVect2[1], fNorms[1]); // get the normal of plane 2

			float fLineVect[3];
			CrossProduct(fNorms[0], fNorms[1], fLineVect);

			float coe = 0.5f;
			glBegin(GL_LINES);
			if (gbFlags[3]) //show ruling1
			{
				glVertex3f(Point[0]-coe*fLineVect[0], Point[1]-coe*fLineVect[1], Point[2]-coe*fLineVect[2]);

				glVertex3f(Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]);
			}
			glEnd();

			//////////////////////////////////////////////////////////////////////////
			// Circles - For cyclide
			coe = (fSphPos2[0]-Point[0])*fLineVect[0]+
				(fSphPos2[1]-Point[1])*fLineVect[1]+
				(fSphPos2[2]-Point[2])*fLineVect[2];
			float fCCener[3]  // center of new circle
			={Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]};

			float dis=FindDis3D(fCCener, fSphPos2);

			float fCRad = sqrt(  dis*dis-fSphRad2*fSphRad2);

			if(dis>fSphRad2)
			{
				float fAngStep=0.032;
				//int num=2*PI/fAngStep; num++;
				float fCirData[300][3];
				int idx=0;

				float q1, q2;
				q1 = acos(fLineVect[2]/sqrt(FindSquare(fLineVect)));
				if (abs(fLineVect[1]) < 0.0001)
					q2 = PI/2;
				else
					q2 = atan(fLineVect[0]/fLineVect[1]); 


				glBegin(GL_LINE_LOOP);
				glColor3f(0.1, 0.1, 0.8); 
				for (float t=0.0; t <2* PI; t += fAngStep)
				{
					float fV[3] = {cos(t), sin(t), 0};  

					//glVertex3fv(fV);
					CirRotation2(fLineVect, fCCener, fCRad, t, fV);

					fCirData[idx][0]=fV[0]; fCirData[idx][1]=fV[1]; fCirData[idx][2]=fV[2];
					idx++;
				}
				if(gbFlags[8]) // the circle families
				{
					glColor3f(0.7f, 0.7f, 0.1f); // circles
					glBegin(GL_LINE_LOOP);

					for (int i=0;i<idx;++i)
					{
						glVertex3fv(fCirData[i]);
					}
					glEnd();

				} 
					 

				if(gbFlags[9])  // after inversion, the cyclide
				{

					float invCen[3], invVect[3], r2;
					CircleInversion(fCCener, fCRad, fLineVect, invCen, &r2, invVect);
					glColor3f(0.7f, 0.7f, 0.1f);
					float **posfull = new float*[idx];
					for (int i=0; i<idx; ++i)
					{
						posfull[i] = new float[3];
					} 
					CirRotation(invVect, invCen, r2, idx, posfull);

					glBegin(GL_LINE_LOOP);
					for (int i=0; i<idx; ++i)
					{
						if(nFlag)
						{

							out15B <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"; 
							out15B <<  posfull[(1+i)%idx][0] <<", " << posfull[(1+i)%idx][1]<<", "<< -posfull[(1+i)%idx][2] <<">, r_edge texture { l_clr_org } } "<<endl;
							outfile <<"v "<<" "<<posfull[i][0]<<" "<<posfull[i][1]<<" "<<posfull[i][2]<<endl;  
						}
						glVertex3fv(posfull[i]);
					}
					glEnd();


					for (int i=0; i<idx; ++i)
					{
						delete []posfull[i];
					}
					delete []posfull;

					if(nFlag)
					for(int i=0; i<idx; ++i)
					{ 

						//outfile2 <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"
						//	<<posfull[(i+1)%idx][0] <<", " << posfull[(i+1)%idx][1]<<", "<< -posfull[(i+1)%idx][2] <<">, r_edge texture { l_clr_org } } "<<endl;

					}
				}

				++nCircles; //number of circles
				num =idx;   //number of vertex along each circle


			}


			//////////////////////////////////////////////////////////////////////////

			// for ruling group2
			if (tmp == 0)
			{ 
				fCroVectors[0][0] = Point[0]-coe*fLineVect[0]; fCroVectors[0][1] = Point[1]-coe*fLineVect[1]; fCroVectors[0][2] = Point[2]-coe*fLineVect[2];
				fCroVectors[1][0] = Point[0]+coe*fLineVect[0]; fCroVectors[1][1] = Point[1]+coe*fLineVect[1]; fCroVectors[1][2] = Point[2]+coe*fLineVect[2];
			}
			if (tmp == 40)
			{			
				fCroVectors[2][0] = Point[0]-coe*fLineVect[0]; fCroVectors[2][1] = Point[1]-coe*fLineVect[1]; fCroVectors[2][2] = Point[2]-coe*fLineVect[2];

				fCroVectors[3][0] = Point[0]+coe*fLineVect[0]; fCroVectors[3][1] = Point[1]+coe*fLineVect[1]; fCroVectors[3][2] = Point[2]+coe*fLineVect[2];
			}

			if (tmp == 80)
			{			
				fCroVectors[4][0] = Point[0]-coe*fLineVect[0]; fCroVectors[4][1] = Point[1]-coe*fLineVect[1]; fCroVectors[4][2] = Point[2]-coe*fLineVect[2];

				fCroVectors[5][0] = Point[0]+coe*fLineVect[0]; fCroVectors[5][1] = Point[1]+coe*fLineVect[1]; fCroVectors[5][2] = Point[2]+coe*fLineVect[2];
			}
			 tmp ++;


		}  

		//saving faces
		if(nFlag && gbFlags[9])
		{
			for (int i=0; i<nCircles-2 ; ++i)
			{ 

				for (int j=0; j<num; ++j)
				{
					outfile <<"f "<<j+i*num+1<<" "<< j+((i+2)%nCircles)*num+1 <<" "<<i*num+(j+1)%num +1 <<endl;
					outfile <<"f "<<j+((i+2)%nCircles)*num+1<<" "<<(j+1)%num+((i+2)%nCircles)*num+1<<" "<<i*num+(j+1)%num +1<<endl;
				} 
			}

		}


		memcpy(fVectors, fCroVectors, 18*sizeof(float));
		
		if (gbFlags[3]) //show ruling1
		{
			glBegin(GL_LINES);
			for (int i=0; i<6; ++i)
			{
				glVertex3fv(fCroVectors[i]);
			}
			glEnd();
		}		



		nCircles =0;
		float data[21]= {0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.25, 0.31,0.42, 0.8, 4, 290000,
			-0.21,-0.22, -0.227, -0.235, -0.242, -0.25, -0.27,  -0.3, -0.5};
		for (int kk=0; kk<21; ++kk)  
		{
				//float h = dab[jj];
				//if(h>0.4)
				//	h+=0.059;
				//if(h>0.6)
				//	h+=0.159;
				
			float k =   data[kk];
			float Point[3] = {k*fVectors[0][0]+(1-k)*fVectors[1][0], k*fVectors[0][1]+(1-k)*fVectors[1][1], k*fVectors[0][2]+(1-k)*fVectors[1][2]};
			float fNorms[2][3];
			float fVect[2][3] = {Point[0]-fVectors[2][0], Point[1]-fVectors[2][1], Point[2]-fVectors[2][2],
				Point[0]-fVectors[3][0], Point[1]-fVectors[3][1], Point[2]-fVectors[3][2]};
			CrossProduct(fVect[0], fVect[1], fNorms[0]); // to get the normal of plane 1
			float fVect2[2][3] = {Point[0]-fVectors[4][0], Point[1]-fVectors[4][1], Point[2]-fVectors[4][2],
				Point[0]-fVectors[5][0], Point[1]-fVectors[5][1], Point[2]-fVectors[5][2]};
			CrossProduct(fVect2[0], fVect2[1], fNorms[1]); // get the normal of plane 2

			float fLineVect[3];
			CrossProduct(fNorms[0], fNorms[1], fLineVect);

			float coe =  0.5f;
			if (gbFlags[4]) //show ruling2
			{
				glColor3f(0.8f, 0.7f, 0.1f); // yellow lines for ruling2
				glBegin(GL_LINES);
				glVertex3f(Point[0]-coe*fLineVect[0], Point[1]-coe*fLineVect[1], Point[2]-coe*fLineVect[2]);
				glVertex3f(Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]);
				glEnd();
			}


			//float fSphRad2 = fSphCenter/2;
			//fSphCenter = fSphPos[2] - fSphRad2;
			//float fSphPos2[3]={0, 0, fSphCenter};
			coe = (fSphPos2[0]-Point[0])*fLineVect[0]+
				(fSphPos2[1]-Point[1])*fLineVect[1]+
				(fSphPos2[2]-Point[2])*fLineVect[2];
			float fCCener[3]  // center of new circle
			={Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]};

			float dis=FindDis3D(fCCener, fSphPos2);

			float fCRad = sqrt(  dis*dis-fSphRad2*fSphRad2);


			if(dis>fSphRad2)
			{
				float fAngStep=0.032;
				int num=2*PI/fAngStep; num++; 
				float **fCirData = new float *[num];
				for (int i=0; i<num; ++i)
				{
					fCirData[i] = new float[3];
				}
				int idx=num;  
				CirRotation(fLineVect, fCCener, fCRad, num, fCirData);
				if(gbFlags[5])
				{
					glColor3f(0.1f, 0.7f, 0.7f); // rulings, group 2, axises of circles
					glBegin(GL_LINE_LOOP);

					for (int i=0;i<idx;++i)
					{
						glVertex3fv(fCirData[i]);
					}
					glEnd();

				}
				if(gbFlags[6])
				{

					glColor3f(0.1f, 0.7f, 0.7f);
					float **posfull= new float *[num];
					for (int i=0; i<num; ++i)
					{
						posfull[i] = new float[3];
					}
					glBegin(GL_LINE_LOOP);
					for (int i=0;i<idx;++i)
					{
						float dis = FindDis3D(fCirData[i], fSphPos);

						float alp = fSphRad*fSphRad/dis/dis; 
						float Pos[3];
						Pos[0] =  (1-alp)*fSphPos[0] + alp * fCirData[i][0] ;
						Pos[1] =  (1-alp)*fSphPos[1] + alp * fCirData[i][1] ;
						Pos[2] =  (1-alp)*fSphPos[2] + alp * fCirData[i][2] ;
						//glVertex3fv(Pos);
					}
					glEnd();


					float invCen[3], invVect[3], r2; 
					CircleInversion(fCCener, fCRad, fLineVect, invCen, &r2, invVect);
					glColor3f(0.9f, 0.7f, 0.1f);
					CirRotation(invVect, invCen, r2, num, posfull);
					glBegin(GL_LINE_LOOP);
					for (int i=0; i<num; ++i)
					{
						glVertex3fv(posfull[i]);
						//if(nFlag)
							//OutF(&outfile3, posfull[i], posfull[(i+1)%idx]);
					}
					glEnd();


					if(nFlag)
					for(int i=0; i<idx; ++i)
					{ 


						out15A <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"; 
						out15A <<  posfull[(1+i)%idx][0] <<", " << posfull[(1+i)%idx][1]<<", "<< -posfull[(1+i)%idx][2] <<">, r_edge texture { l_clr_blue } } "<<endl;
						//outfile2 <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"
						//	<<posfull[(i+1)%idx][0] <<", " << posfull[(i+1)%idx][1]<<", "<< -posfull[(i+1)%idx][2] <<">, r_edge texture { l_clr_blue } } "<<endl;

					}
					for (int i = 0; i<num; ++i)
					{
						delete []posfull[i];
					}
					delete []posfull;
				}

				++nCircles;

				
			}


		}    



	glPopMatrix();

}

void CircleBy2Point2(float *p1, float *p2, float *fCirPos = NULL, 	float *v=NULL, float *radius=NULL)
{
	if (p1[2] <0.0)
	{
		p1[0] = -p1[0]; p1[1] = -p1[1]; p1[2] = -p1[2];
	}
	if (p2[2] <0.0)
	{
		p2[0] = -p2[0]; p2[1] = -p2[1]; p2[2] = -p2[2];
	}

	float qq2[3]={p2[0], p2[1], -p2[2]};



	float a, b,c;
	a=FindDis3D(p1, p2); b=FindDis3D(p2, qq2); c=FindDis3D(p1, qq2);
	float r=a*b*c/b/2/sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));

	a = sqrt(r*r-p1[2]*p1[2]);
	b = sqrt(r*r-p2[2]*p2[2]);
	float pw[3] = {p2[0], p2[1], 0};
	if (FindDis3D(p1, pw) <r)
	{
		a= -a;
	}
	float py[3] = {p1[0], p1[1], 0}; 


	float cen[3] ={ (b*p1[0]+a*p2[0])/(a+b), (b*p1[1]+a*qq2[1])/(a+b), 0};


	float pp1[3];
	float pp2[3];
	for (int i=0; i<3; ++i)
	{
		pp1[i] = p1[i]-cen[i];
		pp2[i] = p2[i]-cen[i];
	}
	CrossProduct(pp1, pp2, v);


	float q1, q2;
	q1 = acos(v[2]);
	if (abs(v[1]) < 0.0001)
		q2 = PI/2;
	else
		q2 = atan(v[0]/v[1]);

	float q0 = atan(v[1]/v[0]); 




	if (NULL != radius)
	{
		(*radius) = r;
	}
	if (NULL != fCirPos)
	{
		fCirPos[0] = cen[0]; fCirPos[1] = cen[1]; fCirPos[2] = cen[2];
	} 

}

void RuledSurfaces()
{
	glPushMatrix();
	glScalef(3.0f, 3.0f, 3.0f);


	//plane: - xy


	if (gbFlags[7])
	{
		glColor3f(1.0f, 1.0f, 1.0f);
		glBegin(GL_QUADS);
		glVertex3f(3.0f, 3.0f, 0.0f);
		glVertex3f(-3.0f, 3.0f, 0.0f);
		glVertex3f(-3.0f, -3.0f, 0.0f);
		glVertex3f(3.0f,  -3.0f, 0.0f);
		glEnd(); 
	}
	float fSphPos[3] = {0.0f, 0.0f,  2.0f};
	fSphPos[2] = gfSphdata[0];
	float fSphRad = gfSphdata[1];
	glPushMatrix();
	glTranslatef(fSphPos[0], fSphPos[1], fSphPos[2]);
	//glutSolidSphere(fSphRad, 20, 20);

	if(gbShowAxis) glutWireSphere(fSphRad, 10, 10);
	glPopMatrix();

	float fVectors[6][3] = 
	{ 0, 0, 0, 2, 2, -0.1,
	0, 0, 1, 3, 1, 1,
	-1, 0, 1.9, 1, -0.8, 2};
	//Inversion of plane
	float fSphCenter = fSphRad*fSphRad/sqrt(FindSquare(fSphPos));
	float fSphRad2 = fSphCenter/2;
	fSphCenter = fSphPos[2] - fSphRad2;
	float fSphPos2[3]={0, 0, fSphCenter};
	glPushMatrix();
	glTranslatef(fSphPos2[0], fSphPos2[1], fSphPos2[2]); 

	//if(gbShowAxis)    glutSolidSphere(fSphRad2, 20, 20);
	glPopMatrix();

	// circle 1
	//float fCirPos1[3] = {1.0f, 0.0f, 0.0f};
	//float theta=PI/3;
	//float rad = 0.6;
	float fAxis[18];
	SphInversion(gfCirPos1, gfRad[0], gfTheta[0], fSphPos, fSphRad, fAxis);
	//memcpy(fVectors , fAxis, 18*sizeof(float));

	//circle 2
	//float fCirPos2[3] = {-0.11, 0.6f, 0.0f};
	//theta= PI/5;
	//rad = 0.6f;
	SphInversion(gfCirPos2, gfRad[1], gfTheta[1], fSphPos, fSphRad, fAxis+6);

	//circle 3
	//float fCirPos3[3] = { -0.71, 0.8f, 0};
	//theta= PI/3;
	//rad = 0.6;
	SphInversion(gfCirPos3, gfRad[2], gfTheta[2], fSphPos, fSphRad, fAxis+12); 
	for(int i=0; i<18; ++i)		fVectors[i/3][i%3]=fAxis[i];




	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_LINES);
	if(gbFlags[2]) //show axis
		for (int i=0; i<6; ++i)
		{
			glVertex3fv(fVectors[i]);
		}
		glEnd();

		float fCroVectors[6][3];
		glColor3f(0.1f, 0.1f, 0.7f); // rulings, group 1
		
		
		
		//////////////////////////////////////////////////////////////////////////
		// for saving obj file...
		//////////////////////////////////////////////////////////////////////////
		
		int nCircles=0;
		int num=0;
		ofstream outfile("cyclide14.obj");
		//outfile2<<"#declare r_edge =  0.002249;"<<endl; 
		ofstream outfile3("lines14A.inc");
		ofstream outfile5("lines14B.inc");



		int tmp = 0; 
		if(nFlag)
		{
			outfile5 <<"#declare r_edge =  0.004149;  "<<endl;
			outfile5<<"#declare edge_tex=l_clr_3"<<endl;

		}

		//outfile3
		for (int kk=1; kk <2; ++kk) 
		for (float h=-0.0; h<0.5f; h+= 0.056)
		{
			 
			 //if(h<0.07)
				//h-= 0.06;
			 //if(h<-0.2)
				// h-=0.11; 

			float k =   1.0/h;
			if(kk==1)
				k = h;
			 
			float Point[3] = {k*fVectors[0][0]+(1-k)*fVectors[1][0], k*fVectors[0][1]+(1-k)*fVectors[1][1], k*fVectors[0][2]+(1-k)*fVectors[1][2]};
			float fNorms[2][3];
			float fVect[2][3] = {Point[0]-fVectors[2][0], Point[1]-fVectors[2][1], Point[2]-fVectors[2][2],
				Point[0]-fVectors[3][0], Point[1]-fVectors[3][1], Point[2]-fVectors[3][2]};
			CrossProduct(fVect[0], fVect[1], fNorms[0]); // to get the normal of plane 1
			float fVect2[2][3] = {Point[0]-fVectors[4][0], Point[1]-fVectors[4][1], Point[2]-fVectors[4][2],
				Point[0]-fVectors[5][0], Point[1]-fVectors[5][1], Point[2]-fVectors[5][2]};
			CrossProduct(fVect2[0], fVect2[1], fNorms[1]); // get the normal of plane 2

			float fLineVect[3];
			CrossProduct(fNorms[0], fNorms[1], fLineVect);

			float coe = 0.5f;
			glBegin(GL_LINES);
			if (gbFlags[3]) //show ruling1
			{
				glVertex3f(Point[0]-coe*fLineVect[0], Point[1]-coe*fLineVect[1], Point[2]-coe*fLineVect[2]);

				glVertex3f(Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]);
			}
			glEnd();

			//////////////////////////////////////////////////////////////////////////
			// Circles - For cyclide
			coe = (fSphPos2[0]-Point[0])*fLineVect[0]+
				(fSphPos2[1]-Point[1])*fLineVect[1]+
				(fSphPos2[2]-Point[2])*fLineVect[2];
			float fCCener[3]  // center of new circle
			={Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]};

			float dis=FindDis3D(fCCener, fSphPos2);

			float fCRad = sqrt(  dis*dis-fSphRad2*fSphRad2);

			if(dis>fSphRad2)
			{
				float fAngStep=0.032;
				//int num=2*PI/fAngStep; num++;
				float fCirData[200][3];
				int idx=0;

				float q1, q2;
				q1 = acos(fLineVect[2]/sqrt(FindSquare(fLineVect)));
				if (abs(fLineVect[1]) < 0.0001)
					q2 = PI/2;
				else
					q2 = atan(fLineVect[0]/fLineVect[1]);

				glBegin(GL_LINE_LOOP);
				glColor3f(0.1, 0.1, 0.8); 
				for (float t=0.0; t < 2.0*PI; t += fAngStep)
				{
					float fV[3] = {cos(t), sin(t), 0};

					//fV[2] = -(cos(t)*fLineVect[0]+sin(t)*fLineVect[1])/fLineVect[2];

					//float dis1 = sqrt(1+fV[2]*fV[2]);  

					//for(int i=0;i<3;++i)
					//{
					//	fV[i]*= fCRad/dis1 ;
					//	fV[i] +=fCCener[i];
					//}
					//fV[0] = fCCener[0]+ fCRad*(cos(q2)*cos(t)+sin(q2)*cos(q1)*sin(t) );
					//fV[1] = fCCener[1]- fCRad*(sin(q2)*cos(t)-cos(q2)*cos(q1)*sin(t) );
					//fV[2] = fCCener[2]+ fCRad*sin(q1)*sin(t); 

					//glVertex3fv(fV);
					CirRotation2(fLineVect, fCCener, fCRad, t, fV);

					fCirData[idx][0]=fV[0]; fCirData[idx][1]=fV[1]; fCirData[idx][2]=fV[2];
					idx++;
				}
				glEnd();
				if(gbFlags[8]) // the circle families
				{
					glColor3f(0.7f, 0.7f, 0.1f); // circles
					glBegin(GL_LINE_LOOP);

					for (int i=0;i<idx;++i)
					{
						glVertex3fv(fCirData[i]);
					}
					glEnd();

				}
				if(gbFlags[9])  // after inversion, the cyclide
				{


					float invCen[3], invVect[3], r2;
					// CircleBy2Point2(pA[0], pA[1], invCen, invVect, &r2);
					// 
					//CirInversion(fLineVect, fCCener, fCRad, fSphPos, fSphRad, invVect, invCen, &r2);

					//CircleInversion(fCCener, fCRad, fLineVect, invCen, &r2, invVect);


					//CirRotation(invVect, invCen, r2, 100, NULL);


					glColor3f(0.7f, 0.7f, 0.1f);
					float pA[2][3];
					glBegin(GL_LINE_LOOP);
					for (int i=0;i<idx;i+=1)
					{

						float dis = FindDis3D(fCirData[i], fSphPos);

						float alp = fSphRad*fSphRad/dis/dis; 
						float Pos[3];
						Pos[0] =  (1-alp)*fSphPos[0] + alp * fCirData[i][0] ;
						Pos[1] =  (1-alp)*fSphPos[1] + alp * fCirData[i][1] ;
						Pos[2] =  (1-alp)*fSphPos[2] + alp * fCirData[i][2] ; 

						//glVertex3fv(Pos);  
						if(i==0)
						{
							pA[0][0] = Pos[0];
							pA[0][1] = Pos[1];
							pA[0][2] = Pos[2];

						}
						if (i== int(idx*0.76))
						{
							pA[1][0] = Pos[0];
							pA[1][1] = Pos[1];
							pA[1][2] = Pos[2];
						}


					}
					glEnd();
					glColor3f(0.1, 0.7, 0.2);
					//CircleBy2Point2(pA[0], pA[1], invCen, invVect, &r2);
					CircleInversion(fCCener, fCRad, fLineVect, invCen, &r2, invVect);

					float **posfull = new float*[idx];
					for (int i=0; i<idx; ++i)
					{
						posfull[i] = new float[3];
					}
					CirRotation(invVect, invCen, r2, idx, posfull);
					glColor3f(0.2, 0.8, 0.3);
					glBegin(GL_LINE_LOOP);
					for (int i=0; i<idx; ++i)
					{
						glVertex3fv(posfull[i]);
						if(nFlag)
						OutF(&outfile5, posfull[i], posfull[(i+1)%idx]);
					}
					glEnd();


					for (int i=0; i<idx; ++i)
					{ 
						delete []posfull[i];// = new float[3];
					}

					delete []posfull;

			

					if(nFlag)
					for(int i=0; i<idx; ++i)
					{ 

						//outfile2 <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"
						//	<<posfull[(i+1)%idx][0] <<", " << posfull[(i+1)%idx][1]<<", "<< -posfull[(i+1)%idx][2] <<">, r_edge texture { l_clr_org } } "<<endl;

					}
				}

				++nCircles; //number of circles
				num =idx;   //number of vertex along each circle


			}


			//////////////////////////////////////////////////////////////////////////

			// for ruling group2
			if (tmp == 0)
			{ 
				fCroVectors[0][0] = Point[0]-coe*fLineVect[0]; fCroVectors[0][1] = Point[1]-coe*fLineVect[1]; fCroVectors[0][2] = Point[2]-coe*fLineVect[2];
				fCroVectors[1][0] = Point[0]+coe*fLineVect[0]; fCroVectors[1][1] = Point[1]+coe*fLineVect[1]; fCroVectors[1][2] = Point[2]+coe*fLineVect[2];
			}
			if (tmp == 4)
			{			
				fCroVectors[2][0] = Point[0]-coe*fLineVect[0]; fCroVectors[2][1] = Point[1]-coe*fLineVect[1]; fCroVectors[2][2] = Point[2]-coe*fLineVect[2];

				fCroVectors[3][0] = Point[0]+coe*fLineVect[0]; fCroVectors[3][1] = Point[1]+coe*fLineVect[1]; fCroVectors[3][2] = Point[2]+coe*fLineVect[2];
			}

			if (tmp == 9)
			{			
				fCroVectors[4][0] = Point[0]-coe*fLineVect[0]; fCroVectors[4][1] = Point[1]-coe*fLineVect[1]; fCroVectors[4][2] = Point[2]-coe*fLineVect[2];

				fCroVectors[5][0] = Point[0]+coe*fLineVect[0]; fCroVectors[5][1] = Point[1]+coe*fLineVect[1]; fCroVectors[5][2] = Point[2]+coe*fLineVect[2];
			}
			 tmp ++;


		}
		//


		if(gbFlags[9])
		{
			if(0)
			for (int i=0; i<nCircles-1; ++i)
			{
				for (int j=0; j<num; ++j)
				{
					outfile <<"f "<<j+i*num+1 <<" "<< j+((i+1)%nCircles)*num+1<<" "<<i*num+(j+1)%num +1 <<endl;
					outfile <<"f "<<j+((i+1)%nCircles)*num+1<<" "<<(j+1)%num+((i+1)%nCircles)*num+1<<" "<<i*num+(j+1)%num +1<<endl;
				} 
			}



		}


		memcpy(fVectors, fCroVectors, 18*sizeof(float));
		
		if (gbFlags[3]) //show ruling1
		{
			glBegin(GL_LINES);
			for (int i=0; i<6; ++i)
			{
				glVertex3fv(fCroVectors[i]);
			}
			glEnd();
		}		







		if(nFlag)
		{
			outfile3<<"// family 2"<<endl;
			outfile3 <<"#declare r_edge =  0.004149;  "<<endl;
		outfile3<<"#declare edge_tex=l_clr_2"<<endl;

		}
		nCircles =0;
		int num0 = 0;
		float dab[10] = {   -1,  -0.4, -0.1, 0.0, 0.05,0.15,0.28,0.38, 0.48,0.6};
		for (int kk=1; kk<2; ++kk) 
		for (float h = -0.46; h<0.60; h += 0.112)
		{
			//float h = dab[jj];
			//if(h>0.4)
			//	h+=0.059;
			//if(h>0.6)
			//	h+=0.159;
			float k=-1.0/h;
			if (kk==1)
			{
				k = h;
			}
			float Point[3] = {k*fVectors[0][0]+(1-k)*fVectors[1][0], k*fVectors[0][1]+(1-k)*fVectors[1][1], k*fVectors[0][2]+(1-k)*fVectors[1][2]};
			float fNorms[2][3];
			float fVect[2][3] = {Point[0]-fVectors[2][0], Point[1]-fVectors[2][1], Point[2]-fVectors[2][2],
				Point[0]-fVectors[3][0], Point[1]-fVectors[3][1], Point[2]-fVectors[3][2]};
			CrossProduct(fVect[0], fVect[1], fNorms[0]); // to get the normal of plane 1
			float fVect2[2][3] = {Point[0]-fVectors[4][0], Point[1]-fVectors[4][1], Point[2]-fVectors[4][2],
				Point[0]-fVectors[5][0], Point[1]-fVectors[5][1], Point[2]-fVectors[5][2]};
			CrossProduct(fVect2[0], fVect2[1], fNorms[1]); // get the normal of plane 2

			float fLineVect[3];
			CrossProduct(fNorms[0], fNorms[1], fLineVect);

			float coe =  0.5f;
			if (gbFlags[4]) //show ruling2
			{
				glColor3f(0.8f, 0.7f, 0.1f); // yellow lines for ruling2
				glBegin(GL_LINES);
				glVertex3f(Point[0]-coe*fLineVect[0], Point[1]-coe*fLineVect[1], Point[2]-coe*fLineVect[2]);
				glVertex3f(Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]);
				glEnd();
			}


			//float fSphRad2 = fSphCenter/2;
			//fSphCenter = fSphPos[2] - fSphRad2;
			//float fSphPos2[3]={0, 0, fSphCenter};
			coe = (fSphPos2[0]-Point[0])*fLineVect[0]+
				(fSphPos2[1]-Point[1])*fLineVect[1]+
				(fSphPos2[2]-Point[2])*fLineVect[2];
			float fCCener[3]  // center of new circle
			={Point[0]+coe*fLineVect[0], Point[1]+coe*fLineVect[1], Point[2]+coe*fLineVect[2]};

			float dis=FindDis3D(fCCener, fSphPos2);

			float fCRad = sqrt(  dis*dis-fSphRad2*fSphRad2);


			if(dis>fSphRad2)
			{
				float fAngStep=0.04;
				int num=2*PI/fAngStep; num++;
				num0 = num;
				float **fCirData = new float *[num];
				for (int i=0; i<num; ++i)
				{
					fCirData[i] = new float[3];
				}
				int idx=num;  
				CirRotation(fLineVect, fCCener, fCRad, num, fCirData);
				if(gbFlags[5])
				{
					glColor3f(0.1f, 0.7f, 0.7f); // rulings, group 2, axises of circles
					glBegin(GL_LINE_LOOP);

					for (int i=0;i<idx;++i)
					{
						glVertex3fv(fCirData[i]);
					}
					glEnd();

				}

				if(gbFlags[6])
				{

					glColor3f(0.1f, 0.7f, 0.7f);
					float **posfull= new float *[num];
					for (int i=0; i<num; ++i)
					{
						posfull[i] = new float[3];
					}


					float pA[2][3];
					glBegin(GL_LINE_LOOP);
					for (int i=0;i<idx;++i)
					{
						float dis = FindDis3D(fCirData[i], fSphPos);

						float alp = fSphRad*fSphRad/dis/dis; 
						float Pos[3];
						Pos[0] =  (1-alp)*fSphPos[0] + alp * fCirData[i][0] ;
						Pos[1] =  (1-alp)*fSphPos[1] + alp * fCirData[i][1] ;
						Pos[2] =  (1-alp)*fSphPos[2] + alp * fCirData[i][2] ;
						//glVertex3fv(Pos); 


						if(i==0)
						{
							pA[0][0] = Pos[0];
							pA[0][1] = Pos[1];
							pA[0][2] = Pos[2];

						}
						if (i== int(idx*0.81))
						{
							pA[1][0] = Pos[0];
							pA[1][1] = Pos[1];
							pA[1][2] = Pos[2];
						}
					}

					glEnd();
					float invCen[3], invVect[3], r2;
					//CircleBy2Point2(pA[0], pA[1], invCen, invVect, &r2);
					CircleInversion(fCCener, fCRad, fLineVect, invCen, &r2, invVect);
					glColor3f(0.9f, 0.7f, 0.1f);
					CirRotation(invVect, invCen, r2, num, posfull);
					glBegin(GL_LINE_LOOP);
					for (int i=0; i<num; ++i)
					{
						glVertex3fv(posfull[i]);
						if(nFlag)
						OutF(&outfile3, posfull[i], posfull[(i+1)%idx]);
					}
					glEnd();

					if(nFlag)
					{   
						 
						//if(nCircles>18)
						//{
						//	for(int i=0; i<idx; ++i)
						//	{ 
						//		int t = (3*idx/2 -i)%idx;

						//		outfile <<"v "<<posfull[t][0]<<" "<<posfull[t][1]<<" "<<posfull[t][2]<<endl;

						//		outfile2 <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"
						//			<<posfull[(i+1)%idx][0] <<", " << posfull[(i+1)%idx][1]<<", "<< -posfull[(i+1)%idx][2] <<">, r_edge texture { l_clr_blue } } "<<endl;

						//	}

						//}
						//else
						for(int i=0; i<idx; ++i)
						{ 

							outfile <<"v "<<posfull[i][0]<<" "<<posfull[i][1]<<" "<<posfull[i][2]<<endl;

							outfile2 <<"cylinder { <" << posfull[i][0] <<", " << posfull[i][1]<<", "<< -posfull[i][2] <<">, <"
								<<posfull[(i+1)%idx][0] <<", " << posfull[(i+1)%idx][1]<<", "<< -posfull[(i+1)%idx][2] <<">, r_edge texture { l_clr_blue } } "<<endl;

						}

					}
					for (int i = 0; i<num; ++i)
					{
						delete []posfull[i];
					}
					delete []posfull;
				}

				++nCircles;

				for (int i=0; i<num; ++i)
				{
					delete []fCirData[i];
				}
				delete []fCirData;

				

				
			}


		}   

		//faces

		if(nFlag)
		{
			for (int i=0; i<nCircles-2 ; ++i)
			{ 

				for (int j=0; j<num0; ++j)
				{
					outfile <<"f "<<j+i*num0+1 <<" "<<i*num0+(j+1)%num0 +1 <<" "<< j+((i+2)%nCircles)*num0+1<<endl;
					outfile <<"f "<<j+((i+2)%nCircles)*num0+1<<" "<<i*num0+(j+1)%num0 +1<<" "<<(j+1)%num0+((i+2)%nCircles)*num0+1<<endl;
				} 
			}
 
		}

		glPopMatrix();
}


void display(void)
{	 


	//init();
	//it may be dangerous...
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

	Cyclides2();

	//RuledSurfaces2();




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

	case 'u':
		fAlpha +=0.05f;
		if (fAlpha>1.0f)
		{
			fAlpha = 1.0f;
		}
		break;
	case 'd':
		fAlpha -=0.05f;
		if (fAlpha<0.0f)
		{
			fAlpha = 0.0f;
		}
		break;
	case 'l':
		t0 -=0.01f;
		if (t0<0.0f)
		{
			t0 += 2*PI;
		}
		break;
	case 'r':
		t0 +=0.05f;
		if (t0>2*PI)
		{
			t0 -= 2*PI;
		}
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
	GLUI_Panel *op_panel = glui_window->add_panel (" ==Sphere== ");

	//  Add the Draw Check box to the 'Object Properties' Panel
	glui_window->add_checkbox_to_panel(op_panel, "Show Axis", &gbShowAxis );
	glui_window->add_checkbox_to_panel(op_panel, "Saving", &nFlag );
	glui_window->add_checkbox_to_panel(op_panel, "Show Sphere", &gbShowSph  );
	glui_window->add_spinner_to_panel(op_panel, "Sphere Pos-x",  GLUI_SPINNER_FLOAT, gfSphPos);
	glui_window->add_spinner_to_panel(op_panel, "Sphere Pos-y",  GLUI_SPINNER_FLOAT, gfSphPos+1);
	glui_window->add_spinner_to_panel(op_panel, "Sphere Pos-z",  GLUI_SPINNER_FLOAT, gfSphPos+2);
	glui_window->add_spinner_to_panel(op_panel, "Sphere Radius",  GLUI_SPINNER_FLOAT, &gfSphR);


	GLUI_Panel *op_panel2 = glui_window->add_panel (" ==Quadric== ");
	//gfAngStep
	glui_window->add_checkbox_to_panel(op_panel2, "Show Cyclide", &gbCirCyc  );

	//glui_window->add_checkbox_to_panel(op_panel2, "Show Cross", &gbCrossLine  );
	glui_window->add_checkbox_to_panel(op_panel2, "Cross Cyclide", &gbCrossCyc  );

	GLUI_Spinner *spinner = 
		glui_window->add_spinner_to_panel(op_panel2, "Angle Step",  GLUI_SPINNER_FLOAT, &gfAngStep);
	spinner->set_float_limits(0.08, PI*3);
	spinner = glui_window->add_spinner_to_panel(op_panel2, "Angle Offset",  GLUI_SPINNER_FLOAT, &gfAngOffSet);
	spinner->set_float_limits(-PI , PI);
	spinner = glui_window->add_spinner_to_panel(op_panel2, "Height",  GLUI_SPINNER_FLOAT, &gfDis);
	spinner->set_float_limits(0.1 , 99);



	GLUI_Panel *op_panel3 = glui_window->add_panel (" ==Inversion== ");
	glui_window->add_checkbox_to_panel(op_panel3, "Show Inversion", &gbShowSp2  );

	glui_window->add_spinner_to_panel(op_panel3, "Radius",  GLUI_SPINNER_FLOAT, &gfr2);

	GLUI_Panel *op_panel4 = glui_window->add_panel (" ==Circle Inversion== ");
 	glui_window->add_checkbox_to_panel(op_panel4, "Show Cir", gbFlags  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Inv Cir", gbFlags+1  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Axis", gbFlags+2  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Rul 1", gbFlags+3  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Rul 2", gbFlags+4  );


	glui_window->add_checkbox_to_panel(op_panel4, "Show orth circles", gbFlags+5  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Cyclide", gbFlags+6  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Plane", gbFlags+7  );

	glui_window->add_checkbox_to_panel(op_panel4, "Show Cir 2", gbFlags+8  );
	glui_window->add_checkbox_to_panel(op_panel4, "Show Cyclide 2", gbFlags+9  );

	spinner = glui_window->add_spinner_to_panel(op_panel4, "Rad",  GLUI_SPINNER_FLOAT, &gfSphdata[0]);
	spinner->set_float_limits(0.2 , 5.0);
	spinner = glui_window->add_spinner_to_panel(op_panel4, "Z",  GLUI_SPINNER_FLOAT, &gfSphdata[1]);
	spinner->set_float_limits(-2.0 , 5.0);



	int kk = 1;
	if (kk==1)
	{
		glui_window = GLUI_Master.create_glui ("4 Points", 0, window_x - 2400, window_y); 
		GLUI_Panel *op_panel5 = glui_window->add_panel (" ==Points Positions== ");

		spinner = glui_window->add_spinner_to_panel(op_panel5, "P1 x", GLUI_SPINNER_FLOAT, &P1[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P1 y", GLUI_SPINNER_FLOAT, &P1[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f); 
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P1 z", GLUI_SPINNER_FLOAT, &P1[2]) ;
		spinner->set_float_limits(0.10f, 2.0f);

		spinner = glui_window->add_spinner_to_panel(op_panel5, "P2 x", GLUI_SPINNER_FLOAT, &P2[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P2 y", GLUI_SPINNER_FLOAT, &P2[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f); 
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P2 z", GLUI_SPINNER_FLOAT, &P2[2]) ;
		spinner->set_float_limits(0.10f, 2.0f);


		spinner = glui_window->add_spinner_to_panel(op_panel5, "P3 x", GLUI_SPINNER_FLOAT, &P3[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P3 y", GLUI_SPINNER_FLOAT, &P3[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f); 
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P3 z", GLUI_SPINNER_FLOAT, &P3[2]) ;
		spinner->set_float_limits(0.10f, 2.0f);


		spinner = glui_window->add_spinner_to_panel(op_panel5, "P4 x", GLUI_SPINNER_FLOAT, &P4[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P4 y", GLUI_SPINNER_FLOAT, &P4[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f); 
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P4 z", GLUI_SPINNER_FLOAT, &P4[2]) ;
		spinner->set_float_limits(0.10f, 2.0f);



		spinner = glui_window->add_spinner_to_panel(op_panel5, "P x", GLUI_SPINNER_FLOAT, &P[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P y", GLUI_SPINNER_FLOAT, &P[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f); 
		spinner = glui_window->add_spinner_to_panel(op_panel5, "P z", GLUI_SPINNER_FLOAT, &P[2]) ;
		spinner->set_float_limits(0.10f, 2.0f);
	}
	else
	{

		//another panel|
		glui_window = GLUI_Master.create_glui ("Cyclides", 0, window_x - 400, window_y+300); 
		GLUI_Panel *op_panel5 = glui_window->add_panel (" ==Circle Data== ");
		//circle 1
		spinner = glui_window->add_spinner_to_panel(op_panel5, "1-pos x", GLUI_SPINNER_FLOAT, &gfCirPos1[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "1-pos y", GLUI_SPINNER_FLOAT, &gfCirPos1[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "1-Angle", GLUI_SPINNER_FLOAT, &gfTheta[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "1-Radius", GLUI_SPINNER_FLOAT, &gfRad[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);

		//circle 2
		spinner = glui_window->add_spinner_to_panel(op_panel5, "2-pos x", GLUI_SPINNER_FLOAT, &gfCirPos2[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "2-pos y", GLUI_SPINNER_FLOAT, &gfCirPos2[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "2-Angle", GLUI_SPINNER_FLOAT, &gfTheta[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "2-Radius", GLUI_SPINNER_FLOAT, &gfRad[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f);

		//circle 3
		spinner = glui_window->add_spinner_to_panel(op_panel5, "3-pos x", GLUI_SPINNER_FLOAT, &gfCirPos3[0]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "3-pos y", GLUI_SPINNER_FLOAT, &gfCirPos3[1]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "3-Angle", GLUI_SPINNER_FLOAT, &gfTheta[2]) ;
		spinner->set_float_limits(-2.0f, 2.0f);
		spinner = glui_window->add_spinner_to_panel(op_panel5, "3-Radius", GLUI_SPINNER_FLOAT, &gfRad[2]) ;
		spinner->set_float_limits(-2.0f, 2.0f);

	}


	glui_window->set_main_gfx_window( main_window );
}


int main(int argc, char** argv)
{

	window_x = (glutGet (GLUT_SCREEN_WIDTH) - 640)/2;
	window_y = (glutGet (GLUT_SCREEN_HEIGHT) - 640)/2;

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
