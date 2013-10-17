/*
 *	Created by
 *	Rafael Umino Nakanishi	nusp 6793198
 *	ICMC USP
 *
 */

#include "Meshadj.hpp"
#include <cmath>
#include <cstdio>

#ifndef MARCHINGTETRA_HPP
#define MARCHINGTETRA_HPP

double function(double x, double y, double z){
	// return x*x + y*y + z*z -1;  //sphere
	return (x*x+y*y+z*z+1-0.25*0.25)*(x*x+y*y+z*z+1 - 0.25*0.25) - 4*(x*x+y*y); //torus
}

Vertex gradient(Vertex point){
	double x = point.x;
	double y = point.y;
	double z = point.z;

	// return Vertex( 2*x, 2*y, 2*z ); //sphere
	return Vertex(	4*x*(x*x+y*y+z*z+0.9375)-8*x, 4*y*(x*x+y*y+z*z+0.9375)-8*y, 4*z*(x*x+y*y+z*z+0.9375) ); //torus
}
double function(Vertex p){
	return function(p.x, p.y, p.z);
}
		
class MarchingTetra{

	private:
		Vertex *bblower, *bbupper;	// Bounding Box
		int slices;					// Tamannho dos cubos 
		double intersect[10];
		int itermax;				// Max de iteracoes para a busca binaria para encontrar o zero da funcao
		double dx, dy, dz;			// fator de deslocamento para os slices

		std::vector< Vertex > cubes; //cubos que tem intereseccao com a superficie
		std::vector< Vertex > tetrahedrons;

		Mesh *mesh;


		Vertex calculateZero(Vertex p1, Vertex p2){
			if(p2 < p1){	Vertex aux = p1; p1 = p2; p2 = aux;	}
			Vertex vetor = p2-p1;
			vetor = vetor/2;

			for (int i = 0; i < itermax; i++) {
				if(function(p1)*function(p1+vetor) > 0) p1 = p1+vetor;
				vetor = vetor/2;
			}
			if(function(p1)*function(p1+vetor) > 0) return p1+vetor;
			return p1;
		}

	public:
		MarchingTetra(Vertex *p1, Vertex *p2){
			bblower = (Vertex*)malloc(sizeof(Vertex));
			bbupper = (Vertex*)malloc(sizeof(Vertex));

			mesh = new Mesh();

			bblower = p1;
			bbupper = p2;
			slices = 30;
			dx = fabs( bbupper->x - bblower->x )/slices;
			dy = fabs( bbupper->y - bblower->y )/slices;
			dz = fabs( bbupper->z - bblower->z )/slices;
			itermax = 10;
		}
		MarchingTetra(){
			*this = MarchingTetra( new Vertex(-2.12451,-2.12451,-2.12451), new Vertex(2.12451,2.12451,2.12451));
		}
		MarchingTetra(double x1, double y1, double z1, double x2, double y2, double z2){
			*this = MarchingTetra( new Vertex(x1,y1,z1), new Vertex(x2,y2,z2));
		}

		void print(){
			printf("%lf %lf %lf\n", bblower->x, bblower->y, bblower->z);
			printf("%lf %lf %lf\n", bbupper->x, bbupper->y, bbupper->z);
			printf("dx %lf, dy %lf, dz %lf\n", dx, dy,dz);
			printf("Slices %d\n", slices);
			printf("Itermax %d\n", itermax);
			printf("\n");
		}

		void printCube(){
			for(int i=0;i<8;i++){
				printf("%3.3lf\t", intersect[i]);
			}
			printf("\n");
		}

		void drawCube(int c){
			if(c<0)c=0; if(c>cubes.size())c=cubes.size();
			glColor3f(0.0,1.0,1.0);
			double i=cubes[c].x;
			double j=cubes[c].y;
			double k=cubes[c].z;
			glBegin(GL_LINE_LOOP); //back
				glVertex3f(i,j,k);
				glVertex3f(i+dx,j,k);
				glVertex3f(i+dx,j+dy,k);
				glVertex3f(i,j+dy,k);
			glEnd();
			glBegin(GL_LINE_LOOP); //left
				glVertex3f(i,j,k);
				glVertex3f(i,j,k+dz);
				glVertex3f(i,j+dy,k+dz);
				glVertex3f(i,j+dy,k);
			glEnd();
			glBegin(GL_LINE_LOOP); //right
				glVertex3f(i+dx,j,k);
				glVertex3f(i+dx,j+dy,k);
				glVertex3f(i+dx,j+dy,k+dz);
				glVertex3f(i+dx,j,k+dz);
			glEnd();
			glBegin(GL_LINE_LOOP); //front
				glVertex3f(i,j,k+dz);
				glVertex3f(i+dx,j,k+dz);
				glVertex3f(i+dz,j+dy,k+dz);
				glVertex3f(i,j+dy,k+dz);
			glEnd();
			glBegin(GL_LINE_LOOP); // top
				glVertex3f(i,j+dy,k);
				glVertex3f(i,j+dy,k+dz);
				glVertex3f(i+dx,j+dy,k+dz);
				glVertex3f(i+dx,j+dy,k);
			glEnd();
			glBegin(GL_LINE_LOOP);//bottom
				glVertex3f(i,j,k);
				glVertex3f(i+dx,j,k);
				glVertex3f(i+dx,j,k+dz);
				glVertex3f(i,j,k+dz);
			glEnd();
		}

		void drawTetrahedron(int t){
			if(t<0)t=0; if(t>tetrahedrons.size())t=tetrahedrons.size();
			Vertex i = tetrahedrons[4*t];
			Vertex j = tetrahedrons[4*t+1];
			Vertex k = tetrahedrons[4*t+2];
			Vertex l = tetrahedrons[4*t+3];
			glColor3f(1,0,0);
			glBegin(GL_LINE_LOOP);
				glVertex3f( i.x, i.y, i.z);
				glVertex3f( j.x, j.y, j.z);
				glVertex3f( k.x, k.y, k.z);
			glEnd();
			glBegin(GL_LINE_LOOP); 
				glVertex3f( i.x, i.y, i.z);
				glVertex3f( j.x, j.y, j.z);
				glVertex3f( l.x, l.y, l.z);
			glEnd();
			glBegin(GL_LINE_LOOP); 
				glVertex3f( i.x, i.y, i.z);
				glVertex3f( k.x, k.y, k.z);
				glVertex3f( l.x, l.y, l.z);
			glEnd();
			glBegin(GL_LINE_LOOP); 
				glVertex3f( j.x, j.y, j.z);
				glVertex3f( l.x, l.y, l.z);
				glVertex3f( k.x, k.y, k.z);
			glEnd();
		}

		void drawCubes(){			
			for(int c=0; c<cubes.size(); c++){
				drawCube(c);
			}
		}

		void setSlices(int slices){
			this->slices = slices;
			dx = fabs( bbupper->x - bblower->x )/slices;
			dy = fabs( bbupper->y - bblower->y )/slices;
			dz = fabs( bbupper->z - bblower->z )/slices;
		}

		void setBB(Vertex *p1, Vertex *p2){
			bblower = p1;
			bbupper = p2;
			setSlices(this->slices);
		}

		void setIterMax(int itermax){
			this->itermax = itermax;
		}

		Mesh *getMesh(){
			return mesh;
		}

		Vertex getPointijk(int cubeCorner, Vertex p0){
			double i = p0.x, j = p0.y, k = p0.z;
			// cubeCorner--;
			switch(cubeCorner){
				case 0:	 return Vertex(	i,	 	j, 		k		); 	
				case 1:	 return Vertex(	i,	 	j,		k+dz	); 	
				case 2:	 return Vertex(	i+dx,	j, 		k+dz	); 	
				case 3:	 return Vertex(	i+dx,	j,		k		); 	
				case 4:	 return Vertex(	i,		j+dy, 	k		); 
				case 5:	 return Vertex(	i,	 	j+dy,	k+dz	); 
				case 6:	 return Vertex(	i+dx, 	j+dy, 	k+dz	); 
				case 7:	 return Vertex(	i+dx,	j+dy,	k		);
			}
		}

		void start(){
			for(double k=bblower->z; k < bbupper->z; k+=dz ){
			for(double j=bblower->y; j < bbupper->y; j+=dy ){
			for(double i=bblower->x; i < bbupper->x; i+=dx ){
				int sign = 0;
				for(int c=0;c<8;c++) intersect[c] = function( getPointijk(c, Vertex(i,j,k) ) );
				if( intersect[0] < 0. )	sign += 1; 
				if( intersect[1] < 0. )	sign += 2; 
				if( intersect[2] < 0. )	sign += 4; 
				if( intersect[3] < 0. )	sign += 8;
				if( intersect[4] < 0. )	sign += 16; 
				if( intersect[5] < 0. )	sign += 32; 
				if( intersect[6] < 0. )	sign += 64; 
				if( intersect[7] < 0. )	sign += 128;
				// printCube();
				// Avaliação de interseccao nos tetraedros
				if( sign > 0 && sign < 255){
					// cubes.push_back(Vertex(i,j,k));
					if( ((sign&139)&139) != 139 ){ 
						evaluateTetra(	getPointijk(0, Vertex(i,j,k)),
						 				getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(3, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(0, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(3, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
					if( ((sign&147)&147) != 147 ){ 
						evaluateTetra(	getPointijk(0, Vertex(i,j,k)),
						 				getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(4, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(0, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(4, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
					if( ((sign&142)&142) != 142 ){ 
						evaluateTetra(	getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(2, Vertex(i,j,k)),
						 				getPointijk(3, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(2, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(3, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
					if( ((sign&198)&198) != 198 ){ 
						evaluateTetra(	getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(2, Vertex(i,j,k)),
						 				getPointijk(6, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(2, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(6, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
					if( ((sign&178)&178) != 178 ){ 
						evaluateTetra(	getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(4, Vertex(i,j,k)),
						 				getPointijk(5, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(4, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(5, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
					if( ((sign&226)&226) != 226 ){ 
						evaluateTetra(	getPointijk(1, Vertex(i,j,k)),
						 				getPointijk(5, Vertex(i,j,k)),
						 				getPointijk(6, Vertex(i,j,k)),
						 				getPointijk(7, Vertex(i,j,k)));
						cubes.push_back( Vertex(i,j,k) );
						tetrahedrons.push_back(getPointijk(1, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(5, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(6, Vertex(i,j,k)));
						tetrahedrons.push_back(getPointijk(7, Vertex(i,j,k)));
					}
				}
			}}} // end for i, for j, for k
			// mesh->printadj();
		}

		void evaluateTetra(Vertex v1, Vertex v2, Vertex v3, Vertex v4){
			int pcount=0;
			bool control[] = {false, false, false, false, false, false};
			Vertex points[5];

			if( function(v1)*function(v2) <= 0 ){
				points[pcount++] = calculateZero(v1, v2);
				control[0] = true;
			}
			if( function(v3)*function(v4) <= 0 ){
				points[pcount++] = calculateZero(v3, v4);
				control[1] = true;
			}

			if( function(v1)*function(v3) <= 0 ){
				points[pcount++] = calculateZero(v1, v3);
				control[2] = true;
			}
			if( function(v2)*function(v4) <= 0 ){
				points[pcount++] = calculateZero(v2, v4);
				control[3] = true;
			}

			if( function(v2)*function(v3) <= 0 ){
				points[pcount++] = calculateZero(v2, v3);
				control[4] = true;
			}
			if( function(v1)*function(v4) <= 0 ){
				points[pcount++] = calculateZero(v1, v4);
				control[5] = true;
			}

			if(pcount == 3)
				mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[2]), &gradient);
			if(pcount == 4){
				if(control[0] && control[1]){
					mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[2]), &gradient);
					mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[3]), &gradient);
				}
				if(control[2] && control[3]){
					mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[2]), &gradient);
					mesh->addFace(new Vertex(points[3]), new Vertex (points[1]), new Vertex (points[2]), &gradient);
				}
				if(control[4] && control[5]){
					mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[2]), &gradient);
					mesh->addFace(new Vertex(points[0]), new Vertex (points[1]), new Vertex (points[3]), &gradient);
				}
			}
		}

		int getfcount(){
			return mesh->getfcount();
		}

		void display(){
			mesh->displayTriangles();
		}

		void display(int k){
			if (k>=mesh->getfcount()) k = mesh->getfcount();
			if (k< 0)k=0;
			mesh->displayTriangles(k);
		}

		void printTriangles(){
			mesh->printTriangles();
		}

		void suavize(){
			for(int i=0; i<10; i++)
				mesh->suavize(&gradient);
		}

		void printHist(){
			mesh->calculateQuality();
			mesh->printHist();
		}
};

#endif