#ifndef PURKINJE_H
#define PURKINJE_H

#include <iostream>
#include <random>			// normal_distribution()
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

#define Iter 1								// Número de iterações de crescimento da árvore.
#define fiber_length 30000.0	// Tamanho das fibras (um)
#define cell_length 100.0			// Tamanho das celulas (um)
#define diameter_cell 25.0		// Diametro da celula (um) (intervalo: 20um <= diameter_cell <= 80um)
#define R_m 25.0							// Raio medio de uma celula.
#define fib_in_bundle 3				// Número de fibras contidas no 'bundle' da raiz.
#define prob_gapJ 0.1					// Probabilidade de ter uma gap junction entre fibras vizinhas (paralela).
#define D 100.0								// Diâmetro do bundle da raiz (um)
#define sigma_c 0.4						// Condutividade da gap junction citoplasmática (default = 0.4 uS/um)
#define G_p 0.5								// Condutância da gap junction plicate (default = 0.5 uS)
#define G_t 1.0								// Condutância da gap junction tranversal (default = 0.5 uS)
#define l 2500								// Tamanho de uma gap junction -> em (um) (intervalo: 2nm < l < 3nm)
#define cc 500								// Cooldown máximo
//#define DEBUG								// Variável de debug

class Edge;
class Node;

// Classe de um volume de controle. É utilizada para armazenar dados da solução da equação do monodomínio.
struct Volume
{
	//int type;					// Tipo de célula.   			### 0 = normal, 1 = gap junction, 2 = gap junction fibra###
	int cooldown;				// Numero de turnos de tempo que a célula ficará inativa se for de estímulo
	double *y_old;			// Resolução da equação no instante n.
	double *y_new;			// Resolução da equação no instante n+1
	double V_star;			// V*

}typedef Volume;

// Classe de um ponto 3D
class Point
{
private:
	double x, y, z;
public:
	Point ();
	Point (double x, double y, double z) { this->x = x; this->y = y; this->z = z;}
	void setPoint (double x, double y, double z) { this->x = x; this->y = y; this->z = z;}
	void setX (double x) { this->x = x; }
	void setY (double y) { this->y = y; }
	void setZ (double z) { this->z = z; }
	double getX () { return (x); }
	double getY () { return (y); }
	double getZ () { return (z); }
};

// Estrutura de nó do grafo
class Node
{
private:
	int id;						// Identificador do nó
	Point *point;			// Coordenadas (x,y,z) do nó
	int fiber_id;			// Identicador da fibra a qual o nó pertence
	int num_edges;		// Contador do número de arestas
	bool stimulus;		// Flag que diz se o nó é de estímulo
	bool gapJ;				// Flag que diz se o nó é de gap junction
	bool retroProg;		// Flag para marcar um nó que houve retroprogração
	Volume *vol;			// Ponteiro para o volume de controle, contém dados sobre a solução da equação do monodomínio
	Node *next;				// Ponteiro para o próximo nó na lista de nós
	Edge *edges;			// Ponteiro para a lista de arestas

public:
// Construtores
	Node ();
	Node (int id, double x, double y, double z, int f_id, bool gap, bool stim);
// Inline
	void setId (int id) { this->id = id; }
	void setPoint (Point *point) { this->point = point; }
	void setIdFiber (int fiber_id) { this->fiber_id = fiber_id; }
	void setNumEdges (int num_edges) { this->num_edges = num_edges; }
	void setStimulus (bool stimulus) { this->stimulus = stimulus; }
	void setGapJ (bool gapJ) { this->gapJ = gapJ; }
	void setRetro (bool retroProg) { this->retroProg = retroProg; }
	void setVolume (Volume *vol) { this->vol = vol; }
	void setNext (Node *next) { this->next = next; }
	void setEdges (Edge *edges) { this->edges = edges; }
	int getId () { return (id); }
	Point* getPoint () { return (point); }
	double getX () { return (point->getX()); }
	double getY () { return (point->getY()); }
	double getZ () { return (point->getZ()); }
	int getIdFiber () { return (fiber_id); }
	int getNumEdges () { return (num_edges); }
	bool getStimulus () { return (stimulus); }
	bool getGapJ () { return (gapJ); }
	bool getRetro () { return (retroProg); }
	Volume* getVolume () { return (vol); }
	Node* getNext () { return (next); }
	Edge* getEdges () { return (edges); }

};

// Estrutura de uma aresta do grafo
class Edge
{
private:
	Node *dest;					// Ponteiro para o nó destino
	double size;				// Tamanho da aresta, norma euclidiana
	double sigma;				// Constante de condutividade da equação de difusão entre dois volumes de controle
	Edge *next;					// Ponteiro para a próxima aresta

public:
	Edge ();
// Inline
	void setDest (Node *dest) { this->dest = dest; }
	void setSize (double size) { this->size = size; }
	void setSigma (double sigma) { this->sigma = sigma; }
	void setNext (Edge *next) { this->next = next; }
	Node* getDest () { return (dest); }
	double getSize () { return (size); }
	double getSigma () { return (sigma); }
	Edge* getNext () { return (next); }

};

// Estrutura do grafo
class Graph
{
private:
	Node *listNodes;			// Ponteiro para a lista de nós
	Node *last_node;			// Ponteiro para último nó da lista de nós
	int total_nodes;			// Contador de nós
	int total_edges;			// Contador de arestas
	double Xmax;					// Valor máximo do domínio no eixo x

// Funções privadas
	Point* setRootFibers ();
	double makeRandom ();
	bool alreadyMake (Point *p, double y, double z, int k);
	double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2);

public:
	Graph (double *y0, int n_eq, double Dx, int exec_number);
// Inline
	void setListNodes (Node *n) { listNodes = n; }
	void setTotalNodes (int tn) { total_nodes = tn; }
	void setTotalEdges (int te) { total_edges = te; }
	void setXmax (int xm) { Xmax = xm; }
	Node* getListNodes () { return (listNodes); }
	int getTotalNodes () { return (total_nodes); }
	int getTotalEdges () { return (total_edges); }
	double getXmax () { return (Xmax); }
// Funções
	void makeRoot (double Dx);
	double* initializeRandomGenerator ();
	Node** setPlotPointers ();
	Node* searchNode (int id);
	void insertNodeGraph (Node *p);
	void insertEdgeGraph (int id_1, int id_2);
	void setLongitidinalGapJunction (double Dx);
	void setTransversalGapJunction (double Dx);
	void setTransversalGapJunction2 (double Dx);
	void calculateVelocity (double t, Node *ptr1, Node *ptr2, double Dt, double Dx, double *maxV, int exec_number);
	void printGraph (int exec_number);
};
#endif
