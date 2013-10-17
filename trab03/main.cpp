/*
 *	Created by
 *	Rafael Umino Nakanishi	nusp 6793198
 *	ICMC USP
 *
 */

#include <vector>
#include <map>
#include "Meshadj.hpp"
#include "marching.hpp"
#include <GL/glut.h>

bool trackingMouse = false;
int W, H;
int oldX, oldY;
Mesh mesh;
MarchingTetra tetra;

int indface;

double pos = 10;
/**/
struct Vec{
	GLfloat x, y, z;
	Vec(){x=y=z=0;}
	Vec(GLfloat x, GLfloat y, GLfloat z){
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void normalize(){
		GLfloat n = sqrt(x*x + y*y + z*z);
		x /= n; y /= n; z /= n;
	}
	
	Vec operator*(const GLfloat s){
		Vec v;
		v.x = x*s;
		v.y = y*s;
		v.z = z*s;
		return v;
	}

	GLfloat operator*(const Vec& v){
		return x*v.x + y*v.y + z*v.z;
	}

	Vec operator^(const Vec& v){
		Vec r;
		r.x = y*v.z - z*v.y;
		r.y = z*v.x - x*v.z;
		r.z = x*v.y - y*v.x;
		return r;
	}

	Vec operator+(const Vec& v){
		Vec r;
		r.x = x + v.x;
		r.y = y + v.y;
		r.z = z + v.z;
		return r;
	}

	void operator=(const Vec& v){
		x = v.x; y = v.y; z = v.z;
	}

	void print(){
		printf("%f %f %f\n", x, y, z);
	}
};

struct Quaternion{
	GLfloat w;
	Vec p;
	Quaternion(){w = p.x = p.y = p.z = 0;}
	Quaternion(GLfloat angle, GLfloat x, GLfloat y, GLfloat z){
		Vec axis; axis.x = x; axis.y = y; axis.z = z;
		set(angle, axis);
	}

	void operator=(Quaternion q){
		w = q.w;
		p = q.p;
	}

	Quaternion operator!(){
		Quaternion q;
		GLfloat m = magnetude();
		q.w = w/m;
		q.p.x = -p.x/m;
		q.p.y = -p.y/m;
		q.p.z = -p.z/m;
		return q;
	}

	Quaternion operator*(Quaternion& q){
		Quaternion r;
		r.w = w*q.w - p*q.p;
		r.p = p*q.w + q.p*w + (p^q.p);
		return r;
	}

	GLfloat magnetude(){
		return w + p*p;
	}

	void set(GLfloat angle, Vec axis){
		angle *= 0.5;
		axis.normalize();

		GLfloat PI = 2.0*acos(0.0);
		GLfloat sin_theta = sin(angle*PI/180.0);
		GLfloat cos_theta = cos(angle*PI/180.0);

		w = cos_theta;
		p.x = sin_theta*axis.x;
		p.y = sin_theta*axis.y;
		p.z = sin_theta*axis.z;

		GLfloat n = 1.0f/sqrt(p.x*p.x+p.y*p.y+p.z*p.z+w*w);
		p.x *= n;
		p.y *= n;
		p.z *= n;
		w *= n;
	}
	void matrix(GLfloat *R){
		R[0] = 1.0f - 2.0f*p.y*p.y - 2.0f*p.z*p.z;
		R[1] = 2.0f*p.x*p.y - 2.0f*p.z*w;
		R[2] = 2.0f*p.x*p.z + 2.0f*p.y*w;
		R[3] = 0.0f;

		R[4] = 2.0f*p.x*p.y + 2.0f*p.z*w;
		R[5] = 1.0f - 2.0f*p.x*p.x - 2.0f*p.z*p.z;
		R[6] = 2.0f*p.y*p.z - 2.0f*p.x*w;
		R[7] = 0.0f;

		R[8] = 2.0f*p.x*p.z - 2.0f*p.y*w;
		R[9] = 2.0f*p.y*p.z + 2.0f*p.x*w;
		R[10] = 1.0f - 2.0f*p.x*p.x - 2.0f*p.y*p.y;
		R[11] = 0.0f;

		R[12] = 0.0f;
		R[13] = 0.0f;
		R[14] = 0.0f;
		R[15] = 1.0f;
	}
};
int degree;
Quaternion rotation(0,1,0,0);
Quaternion keyboardXRotation(5,1,0,0);
Quaternion keyboardYRotation(5,0,1,0);
Quaternion keyboardZRotation(5,0,0,1);
Quaternion keyboardXiRotation(-5,1,0,0);
Quaternion keyboardYiRotation(-5,0,1,0);
Quaternion keyboardZiRotation(-5,0,0,1);

void motion(int x, int y )
{
	if(trackingMouse){
		GLfloat r = 1.0;

		GLfloat dX = (x - oldX)/r;
		GLfloat dY = (y - oldY)/r;

		GLfloat X = r*dX/W;
		GLfloat Y = -r*dY/H;
		Vec start(X, Y, 1.0 - sqrt(X*X + Y*Y));
		start.normalize();
		X = -r*dX/W;
		Y = r*dY/H;
		Vec end(X, Y, 1.0 - sqrt(X*X + Y*Y));
		end.normalize();
		Vec axis = start^end;
		axis.normalize();
		GLfloat angle = (end*start);
		angle*=2;
		Quaternion q(angle,axis.x,axis.y,axis.z);
		rotation = rotation*q;
		glutPostRedisplay(); 

		oldX = x; oldY = y;
	}
}

void render(void)
{
	glClearColor( 1.0, 1.0, 1.0, 0.0 );
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // glEnable( GL_DEPTH_TEST );
	glLoadIdentity();
	
	
	GLfloat ma[16];
	rotation.matrix(ma);
	glMultMatrixf(ma);
	tetra.display(indface);
	
	glFlush();
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	W = w;
	H = h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, w/h, 0, 1000); 
	gluLookAt( 0,0,pos,
		   0,0,0,
		   0,1,0);
	glMatrixMode(GL_MODELVIEW);
}

void zoomin(){

}

void zoomout(){

}

void mouse(int button, int button_state, int x, int y )
{
	if(button == GLUT_LEFT_BUTTON && button_state == GLUT_UP){
		trackingMouse = false;
	}
	else if(button == GLUT_LEFT_BUTTON){
		trackingMouse = true;
		oldX = x; oldY = y;
	}
	else if(button == GLUT_RIGHT_BUTTON && button_state == GLUT_UP){
		trackingMouse = false;
	}
	else if(button == GLUT_RIGHT_BUTTON){
		trackingMouse = true;
	}
	else if(button == GLUT_MIDDLE_BUTTON && button_state == GLUT_UP){
		trackingMouse = false;
	}
	else if(button == GLUT_MIDDLE_BUTTON){
		trackingMouse = true;
	}
}


void keyboard(unsigned char key, int x, int y){
	switch(key){
		case 'a': rotation = keyboardXRotation*rotation; break;
		case 's': rotation = keyboardYRotation*rotation; break;
		case 'd': rotation = keyboardZRotation*rotation; break;
		case 'A': rotation = keyboardXiRotation*rotation; break;
		case 'S': rotation = keyboardYiRotation*rotation; break;
		case 'D': rotation = keyboardZiRotation*rotation; break;
		case 'z': pos+=0.1;reshape(800,800); break;
		case 'x': pos-=0.1;reshape(800,800); break;
		case 'k': indface++; break;
		case 'K': indface--; break;
		case 't': indface = tetra.getfcount(); break;
		case 'h': indface += tetra.getfcount()/2; break;
		case 'H': indface -= tetra.getfcount()/2; break;
		case 'j': indface += tetra.getfcount()/4; break;
		case 'J': indface -= tetra.getfcount()/4; break;
		case 'p': tetra.printTriangles(); break;
		case 'o': tetra.suavize(); break;
		case 'b': tetra.printHist(); break;
		case 'q': exit(0);
	}
	glutPostRedisplay();
}
/**/
int main(int argc, char** argv){
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH |GLUT_RGBA );
	glutInitWindowSize(800, 800);
	glutCreateWindow("Marching Tetrahedra");
	glutReshapeFunc(reshape);
	glutDisplayFunc(render);
	glutKeyboardFunc(keyboard);
	glutMouseFunc( mouse );
	glutMotionFunc( motion );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );
    glEnable( GL_CULL_FACE );
	
	indface = -1;
	degree = 5;

	tetra.setSlices(30);
	tetra.setBB(new Vertex(-2,-2,-2), new Vertex(2,2,2));
	tetra.setIterMax(10);
	tetra.start();

	// tetra.display();

 
    // Vertex   v1 = Vertex(1., 0., 0.),
    //      v2 = Vertex(0., 1., 0.),
    //      v3 = Vertex(0., 0., 1.),
    //      v4 = Vertex(0.5, 0., 0.);
 
    // // std::cout << "plo: "; v1.printvert();
 
    // mesh.addFace(&v2, &v4, &v3);
    // mesh.addFace(&v4, &v2, &v1);
    // mesh.addFace(&v3, &v1, &v2);
    // mesh.addFace(&v1, &v3, &v4);
 
	// Vertex 	v1  = Vertex(0,  -0.525731,  0.850651),
	// 		v2  = Vertex(0.850651,  0,  0.525731),
	// 		v3  = Vertex(0.850651,  0 , -0.525731),
	// 		v4  = Vertex(-0.850651,  0,  -0.525731),
	// 		v5  = Vertex(-0.850651,  0,  0.525731),
	// 		v6  = Vertex(-0.525731,  0.850651,  0),
	// 		v7  = Vertex(0.525731,  0.850651,  0),
	// 		v8  = Vertex(0.525731 , -0.850651 , 0),
	// 		v9  = Vertex(-0.525731,  -0.850651,  0),
	// 		v10 = Vertex(0,  -0.525731 , -0.850651),
	// 		v11 = Vertex(0,  0.525731,  -0.850651),
	// 		v12 = Vertex(0 , 0.525731 , 0.850651);

	// mesh.addFace(&v2,  &v3, &v7);
	// mesh.addFace(&v2,  &v8, &v3);
	// mesh.addFace(&v4,  &v5, &v6);
	// mesh.addFace(&v5,  &v4, &v9);
	// mesh.addFace(&v7,  &v6, &v12);
	// mesh.addFace(&v6,  &v7, &v11);
	// mesh.addFace(&v10,  &v11, &v3);
	// mesh.addFace(&v11,  &v10, &v4);
	// mesh.addFace(&v8,  &v9, &v10);
	// mesh.addFace(&v9,  &v8, &v1);
	// mesh.addFace(&v12,  &v1, &v2);
	// mesh.addFace(&v1,  &v12, &v5);
	// mesh.addFace(&v7,  &v3, &v11);
	// mesh.addFace(&v2,  &v7, &v12);
	// mesh.addFace(&v4,  &v6, &v11);
	// mesh.addFace(&v6,  &v5, &v12);
	// mesh.addFace(&v3,  &v8, &v10);
	// mesh.addFace(&v8,  &v2, &v1);
	// mesh.addFace(&v4,  &v10, &v9);
	// mesh.addFace(&v5,  &v9, &v1);

	mesh.displayTriangles();

		// std::cout << "plo: "; v1.printvert();

	glutMainLoop();
}