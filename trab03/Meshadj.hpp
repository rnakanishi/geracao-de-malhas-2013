/*
 *	Created by
 *	Rafael Umino Nakanishi	nusp 6793198
 *	ICMC USP
 *
 */

#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <fstream>
#include <cmath>
#include <GL/glut.h>

#ifndef MESH_HPP
#define MESH_HPP

#define eps 1e-6

struct Vertex;
struct Halfedge;
struct Face ;

struct Vertex{
	int index;
	double x, y, z;
	
	Vertex(){
		x = 0; y=0; z=0;
	}
	Vertex(double x, double y, double z){
		*this = Vertex();
		this->x = (x);
		this->y = (y);
		this->z = (z);
	}
	Vertex(const Vertex &v){
		*this = Vertex(v.x, v.y, v.z);
		index = v.index;
	}
	double norm(){
		return sqrt(x*x + y*y + z*z);
	}
	double dot(Vertex v){
		return (x*v.x) + (y*v.y) + (z*v.z);
	}
	Vertex cross (Vertex v) {
		return Vertex(y*v.z-v.y*z, z*v.x-x*v.z, x*v.y-y*v.x); 
	}
	void printvert(){
		std::cout << x << "  " << y << "  " << z /*<< std::endl*/;
	}
	bool operator == (const Vertex &v) const{
		return (x-eps>=v.x || x+eps<=v.x) && (y-eps>=v.y || y+eps<=v.y) && (z-eps>=v.z || z+eps<=v.z);
	}
	bool operator < (const Vertex &v) const{
		if(x < v.x) return true;
		if(x == v.x && y < v.y) return true;
		if(x == v.x && y == v.y && z < v.z) return true;
		return false;
	}
	void operator =(const Vertex &v){
		x = v.x;
		y = v.y;
		z = v.z;
		index = v.index;
	}
	Vertex operator -(const Vertex v){
		return Vertex(x-v.x, y-v.y, z-v.z);
	}
	Vertex operator +(Vertex v){
		return Vertex(x+v.x, y+v.y, z+v.z);	
	}
	Vertex operator /(const double d){
		return Vertex(x/d, y/d, z/d);
	}
	Vertex operator *(const double d){
		return Vertex(x*d, y*d, z*d);
	}
};

struct Face {
	int vertices[3];
	Face(){
	}
};

/***************************/
/*  Mesh Class             */
/***************************/
class Mesh{
	private:	
		std::vector< std::list<int> > adjacent;	// Matrix de ajacencia (topologia)
		std::vector<Vertex*> verts;				// Vetor da posicao dos vertice (geometria)
		std::vector<Face*> faces;				// Controle dos vertices das faces
		std::map< Vertex, int> vertmap;			// Indexacao dos vertices
		std::vector<Vertex*> vcopy;	
		int hist[11];

		int vcount, fcount;
		
		int addVertex(Vertex*);
		void addToAdjacentMatrix(int, int);
		void copyVertices(int);
		Vertex calculateCentroid(int);

	public:
		Mesh();
		int getfcount();
		void addFace(Vertex*, Vertex*, Vertex*, Vertex (*)(Vertex));
		void displayTriangles();
		void displayTriangles(int);
		void printTriangles();
		void suavize(Vertex (*)(Vertex));
		void calculateQuality();
		void printHist();
		void addToHist(double);
		void printvmap(){
			std::map<Vertex, int>::iterator it;
			for(it=vertmap.begin(); it!=vertmap.end(); it++){
				Vertex aux = it->first;
				std::cout<<it->second<< "   ";
				aux.printvert();
			}
		}
		void printadj(){
			std::list<int>::iterator it;
			for(int i=0; i<vcount; i++){
				printf("#Vertice %d: ", i);
				for(it=adjacent[i].begin(); it!=adjacent[i].end(); it++){
					printf("%d -> ", *it);
				}
				printf("\n" );
			}
		}
};

Mesh::Mesh(){
	vcount = fcount = 0;
}

int Mesh::addVertex(Vertex *v){
	std::map<Vertex, int>::iterator it;

	it = vertmap.find(*v);
	if(it == vertmap.end()){

		v->index = vcount;
		vertmap[*v] = vcount++;
		verts.push_back(v);
		return (vcount-1);
	}
	return it->second;
}

int Mesh::getfcount(){
	return fcount;
}

void Mesh::addToAdjacentMatrix(int i, int j){
	// printf("#Adding to adj matrix. %d->%d  size: %d Vcount: %d\n", i,j,adjacent.size(), vcount);
	if( adjacent.size() <= i ){ // Novo vertice que foi inserido
		// printf("# New vertex: %d->%d\n", i, j);
		adjacent.push_back( std::list<int>() );
		adjacent[i].push_back(j);
		return;
	}

	std::list<int>::iterator it;
	for( it=adjacent[i].begin(); it!=adjacent[i].end(); it++){
		
		if( j == *it ) {
			// printf("# Already on adj\n");
			return; // ja era adjacente
		}
	}
	// printf("# Added new adj %d->%d\n",i,j);
	adjacent[i].push_back(j);

}

void Mesh::addFace(Vertex *v1, Vertex *v2, Vertex *v3, Vertex (*gradient)(Vertex)){

	Vertex ntri;
	Vertex dir1, dir2, bar;
	dir1 = (*v1) - (*v2);
	dir2 = (*v2) - (*v3);
	ntri = dir1.cross(dir2);
	bar = ((*v1)+(*v2)+(*v3))/3;
	// Garantir orietacao dos triangulos para fora
	if(ntri.dot(gradient(bar)) < 0){
		Vertex* aux = v2;
		v2 = v3;
		v3 = aux;
	}

	int indv1 = addVertex(v1);
	int indv2 = addVertex(v2);
	int indv3 = addVertex(v3);
	

	// printf("# Face: %d %d %d\n", indv1,indv2,indv3);

	addToAdjacentMatrix(indv1, indv2);
	addToAdjacentMatrix(indv1, indv3);
	addToAdjacentMatrix(indv2, indv1);
	addToAdjacentMatrix(indv2, indv3);
	addToAdjacentMatrix(indv3, indv1);
	addToAdjacentMatrix(indv3, indv2);

	Face *face;
	face = (Face*)malloc(sizeof(Face));

	face->vertices[0] = indv1;
	face->vertices[1] = indv2;
	face->vertices[2] = indv3;

	faces.push_back(face);

	fcount++;
}

void Mesh::displayTriangles(){

	// std::cout << std::endl;
	/**/
	glColor4f(0.3,0.3,0.3,1);
	glBegin(GL_TRIANGLES);
	for(int i=0; i<fcount; i++){
		for(int j=0;j<3;j++){
				Vertex *v = verts[ faces[i]->vertices[j] ];
				double v1 = v->x;
				double v2 = v->y;
				double v3 = v->z;
				glVertex3f(v1,v2,v3);
		}
	}
	glEnd();

	glColor4f(0.8,0.8,0.8,1);
	for(int i=0; i<fcount; i++){
		glBegin(GL_LINE_LOOP);
		for(int j=0;j<3;j++){
			Vertex *v = verts[ faces[i]->vertices[j] ];
				double v1 = v->x;
				double v2 = v->y;
				double v3 = v->z;
				glVertex3f(v1,v2,v3);
			}
		}
		glEnd();
	}/**/

void Mesh::displayTriangles(int k){
	glColor4f(0.3,0.3,0.3,1);
	glBegin(GL_TRIANGLES);
	for(int i=0; i<k; i++){
		for(int j=0;j<3;j++){
				Vertex *v = verts[ faces[i]->vertices[j] ];
				double v1 = v->x;
				double v2 = v->y;
				double v3 = v->z;
				glVertex3f(v1,v2,v3);
		}
	}
	glEnd();

	glColor4f(0.8,0.8,0.8,1);
	for(int i=0; i<k; i++){
		glBegin(GL_LINE_LOOP);
		for(int j=0;j<3;j++){
			Vertex *v = verts[ faces[i]->vertices[j] ];
				double v1 = v->x;
				double v2 = v->y;
				double v3 = v->z;
				glVertex3f(v1,v2,v3);
		}
		glEnd();
	}
}

void Mesh::printTriangles(){
	std::cout << "#Total vertices: " << vcount << std::endl;
	std::cout << "#Total face: " << fcount << std::endl;

	for(int i=0; i<vcount; i++){
		std::cout/* << verts[i]->index*/ << "v " ;
		verts[i]->printvert(); printf("\n");
	}

	std::cout << std::endl;
	for(int i=0; i<fcount; i++){
		printf("f ");
		for(int j=0; j<3; j++){
			printf("%d ", faces[i]->vertices[j]+1);
		}
		printf("\n");
	}
}

void Mesh::copyVertices(int target){
	switch(target){
		case 0:		//from verts to vcopy
			for(int i=0; i<vcount; i++){
				vcopy.push_back( new Vertex( *verts[i] ) );
			}
			break; 
		case 1:		// from vcopy to verts
			verts.clear();
			for(int i=0; i<vcount; i++){
				verts.push_back( new Vertex( *vcopy[i] ) );
			}
	}

}

Vertex Mesh::calculateCentroid( int index ){
	std::list<int>::iterator it;
	Vertex b;
	int count = 0;

	it = adjacent[index].begin();
	// printf("#Centroide de: \n");
	for( ; it!=adjacent[index].end(); it++){
		// printf("\t#");b.printvert();printf("\n");
		b = b + *verts[*it];
		count++;
	}
	b = b/count;
	// printf("#\tResultadndo em: ");b.printvert();printf("---\n\n");
	return b;
}

void Mesh::suavize(Vertex (*gradient)(Vertex)){
	copyVertices(0);
	Vertex normal, baricentro, desloc;
	Vertex *toBeChanged;
	double alpha = 0.1;

	for (int i = 0; i < vcount; ++i) {
		normal = gradient( *verts[i] );
		baricentro = calculateCentroid(i);
		desloc = baricentro - *verts[i];

		// printf("#Normal: "); normal.printvert(); printf(" -> %lf ", sqrt(normal.dot(normal)));		printf("\n");
		// printf("#Baricentro: "); baricentro.printvert();  printf(" -> %lf ", sqrt(baricentro.dot(baricentro)));		printf("\n");
		// printf("#Deslocamento: "); desloc.printvert();  printf(" -> %lf ", sqrt(desloc.dot(desloc)));		printf("\n");
		// printf("\n");

		normal = normal / sqrt(normal.dot(normal));
		// baricentro = baricentro / sqrt(baricentro.dot(baricentro));
		// desloc = desloc / sqrt(desloc.dot(desloc));

		toBeChanged = vcopy[i];
		*toBeChanged = (*verts[i]) + ( desloc - normal*(desloc.dot(normal)) )*alpha;
	}
	copyVertices(1);
}

void Mesh::addToHist(double value){
	printf("%lf ", value);
	if( value > 0.9 ){ hist[9]++; return;}
	if( value > 0.8 ){ hist[8]++; return;}
	if( value > 0.7 ){ hist[7]++; return;}
	if( value > 0.6 ){ hist[6]++; return;}
	if( value > 0.5 ){ hist[5]++; return;}
	if( value > 0.4 ){ hist[4]++; return;}
	if( value > 0.3 ){ hist[3]++; return;}
	if( value > 0.2 ){ hist[2]++; return;}
	if( value > 0.1 ){ hist[1]++; return;}
	hist[0]++;
}

void Mesh::printHist(){
	// printf("#\tHistograma:\n#\t");
	// for(int i=0; i<10; i++){
	// 	printf("%d ", hist[i]);
	// }
	// printf("\n");
}

void Mesh::calculateQuality(){
	double k = 4*sqrt(3);
	double a,b,c,s,area,ratio;
	std::ofstream file;
	Vertex *v1,*v2,*v3;

	file.open("histogram.data", std::ios::out);
	for(int i=0; i<10; i++){
		hist[i] = 0;
	}

	for(int i=0; i<fcount; i++){
		v1 = verts[faces[i]->vertices[0]];
		v2 = verts[faces[i]->vertices[1]];
		v3 = verts[faces[i]->vertices[2]];

		// printf("#Vertices: \n");
		// printf("# ");v1->printvert(); printf("\n");
		// printf("# ");v2->printvert(); printf("\n");
		// printf("# ");v3->printvert(); printf("\n");

		a = ((*v1) - (*v2)).norm();
		b = ((*v2) - (*v3)).norm();
		c = ((*v3) - (*v1)).norm();
		s = (a+b+c)/2;
		area = sqrt( s*(s-a)*(s-b)*(s-c) );
		ratio = (k*area)/(a*a + b*b + c*c);

		file << ratio << " ";

		// addToHist( ratio );
	}
	file.close();
}

#endif