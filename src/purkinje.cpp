#include "purkinje.h"
#include "plot.h"

Node::Node () {}
Edge::Edge () {}

bool flag;				// Flag utilizado no cálculo da velocidade de propagação
int noChange;				// Contador da velocidade de propagação
int c_nodes;				// Contador de nós
int c_segments;				// Contador de segmentos
double *y_0;				// Vetor das condições inicias
int num_eq;				// Número de equações do modelo celular
Node **start_fiber;			// Vetor que armazena os volumes inicias de cada fibra

// Construtor de um Node
Node::Node (int id, double x, double y, double z, int f_id, bool gap, bool stim)
{
	int j;

	this->id = id;
	this->x = x; this->y = y; this->z = z;
	this->fiber_id = f_id;
	this->gapJ = gap;
	this->stimulus = stim;
	this->num_edges = 0;
	this->retroProg = false;
	this->next = NULL;
	this->edges = NULL;
	cell = new Cell();
	// É célula de gap junction longitudinal ?
	if (gap)
	{
		cell->type = 1;
		cell->sigma = G_p;
	}
	else
	{
		cell->type = 0;
		cell->sigma = sigma_c;
	}
	retroProg = false;
	// Aloca o vetor de solução do volume de controle
	cell->y = new double[num_eq]();
	cell->y_new = new double[num_eq]();
	// Se for volume de estímulo colocar as condições iniciais normalmente
	if (stim)
	{
		for (j = 0; j < num_eq; j++)
			cell->y[j] = y_0[j];
		cell->cooldown = rand() % cc;
	}
	// Senão colocar condições inicias somente nas variáveis gate
	else
	{
		cell->y[0] = y_0[num_eq];
		for (j = 1; j < num_eq; j++)
			cell->y[j] = y_0[j];
	}
	c_nodes++;
}

// Função que calcula a a distância euclidiana entre dois pontos do grafo
double Graph::calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)
{
	return (sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2)));
}

// Função que busca um nó no grafo
Node* Graph::searchNode (int id)
{
	Node *ptr = getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getId() == id)
			return (ptr);
		else
			ptr = ptr->getNext();
	}
	cout << "[-]ERROR! In searchNode(), node was not found. Index = " << id << endl;
	return (ptr);
}

// Procedimento que insere uma aresta no grafo
void Graph::insertEdgeGraph (int id_1, int id_2)
{
	Node *ptr1 = searchNode(id_1);
	Node *ptr2 = searchNode(id_2);

	// Setar os dados da aresta
	Edge *edge = new Edge();
	edge->setId(ptr2->getId());
	edge->setSize(calcNorm(ptr1->getX(),ptr1->getY(),ptr1->getZ(),ptr2->getX(),ptr2->getY(),ptr2->getZ()));
	edge->setX(ptr2->getX());
	edge->setY(ptr2->getY());
	edge->setZ(ptr2->getZ());
	edge->setCell(ptr2->getCell());
	edge->setNext(NULL);
	edge->setDest(ptr2);

	// Incrementa o contador de arestas do nó origem
	ptr1->setNumEdges(ptr1->getNumEdges()+1);

	// Insere a aresta no grafo
	if (ptr1->getEdges() == NULL)
		ptr1->setEdges(edge);
	else
	{
		Edge *ptrl = ptr1->getEdges();
		while (ptrl->getNext() != NULL)
			ptrl = ptrl->getNext();
		ptrl->setNext(edge);
	}
	c_segments++;
}

// Procedimento que insere um nó no grafo
void Graph::insertNodeGraph (Node *p)
{
	if (p != NULL)
	{
		// Primeiro nó do grafo
		if (getListNodes() == NULL)
		{
			setListNodes(p);
			last_node = p;
		}
		// Insere  no último elemento da lista
		else
		{
			Node *ptr = last_node;
			ptr->setNext(p);
			last_node = ptr->getNext();
		}
	}
	else
		cout << "[-]ERROR! In insertNodeGraph(), pointer 'p' is invalid." << endl;
}

// Gera um número aleatório no quadrado entre -D e D
double Graph::makeRandom ()
{
	int signal;
	double r;
	signal = rand() % 2;
	if (signal)
		r = rand() / (double)RAND_MAX*D;
	else
		r = -1*(rand() / (double)RAND_MAX*D);
	return (r);
}

// Verifica se o ponto da raiz inicial já foi feito
bool Graph::alreadyMake (Point *p, double y, double z, int k)
{
	int i;
	for (i = 0; i < k; i++)
	{
		if (p[i].y == y && p[i].z == z)
			return (true);
	}
	return (false);
}

// Seta os pontos que irão ser a raiz do Bundle
Point* Graph::setRootFibers ()
{
	int i, k;
	double norm;
	Point *point = new Point[fib_in_bundle-1];
	srand(time(NULL));
	for (i = 0, k = 0; i < fib_in_bundle-1; i++)
	{
		do
		{
			point[i].y = makeRandom();
			point[i].z = makeRandom();
			norm = sqrt(point[i].y*point[i].y + point[i].z*point[i].z);
		}while (alreadyMake(point,point[i].y,point[i].z,k) && norm >= D);
		k++;
	}
	return (point);
}

// Função para setar as gap junctions transversais
void Graph::setGapJunctionFibers ()
{
	int i;
	int id_1, id_2, f_id;
	double r;
	cout << "======== Fiber starts in =========" << endl;
	for (i = 0; i < fib_in_bundle; i++)
		cout << start_fiber[i]->getId() << " ";
	cout << endl;
	// Pulamos o primeiro Bundle
	Node *ptr, *ptr2;
	ptr = getListNodes()->getNext();
	// Percorrer a fibra do meio
	while (ptr->getIdFiber() == 1)
	{
		// Gap junction entre fibras ?
		r = (double)rand() / (double)RAND_MAX;
		if (r <= prob_gapJ)
		{
			id_1 = ptr->getId();
			// Qual fibra terá gap junction ?
			do {
				f_id = rand() % fib_in_bundle;
			}
			while (f_id == 0);
			ptr2 = start_fiber[f_id];
			while (ptr2 != NULL && ptr2->getId() != start_fiber[f_id]->getId() + id_1 - 1)
				ptr2 = ptr2->getNext();
			if (ptr2 != NULL)
			{
				id_2 = ptr2->getId();
				cout << "Gap Junction between " << id_1 << " -- " << id_2 << endl;
				// Alterar a condutividade dos volumes de controle envolvidos
				ptr->getCell()->sigma = ptr->getCell()->sigmaY = G_t;
				ptr2->getCell()->sigma = ptr->getCell()->sigmaY = G_t;
				ptr->getCell()->type = ptr2->getCell()->type = 2;
				insertEdgeGraph(id_1,id_2);
				insertEdgeGraph(id_2,id_1);
			}
		}
		ptr = ptr->getNext();
	}
}

// Função para setar as gap junctions transversais
void Graph::setGapJunctionFibers2 ()
{
	int i;
	int id_1, id_2, f_id;
	double r;
	cout << "======== Fiber starts in =========" << endl;
	for (i = 0; i < fib_in_bundle; i++)
		cout << start_fiber[i]->getId() << " ";
	cout << endl;
	// Pulamos o primeiro Bundle
	Node *ptr, *ptr2;
	ptr = getListNodes()->getNext();
	// Percorrer a fibra do meio
	while (ptr->getIdFiber() == 1)
	{
		// Gap junction entre fibras ?
		r = (double)rand() / (double)RAND_MAX;
		if (r <= prob_gapJ)
		{
			id_1 = ptr->getId();
			// Qual fibra terá gap junction ?
			do {
				f_id = rand() % (fib_in_bundle-1);
			}
			while (f_id == 0);
			ptr2 = start_fiber[f_id];
			while (ptr2 != NULL && ptr2->getId() != start_fiber[f_id]->getId() + id_1 - 1)
				ptr2 = ptr2->getNext();
			if (ptr2 != NULL)
			{
				id_2 = ptr2->getId();
				cout << "Gap Junction between " << id_1 << " -- " << id_2 << endl;
				// Alterar a condutividade dos volumes de controle envolvidos
				ptr->getCell()->sigma = ptr->getCell()->sigmaY = G_t;
				ptr2->getCell()->sigma = ptr->getCell()->sigmaY = G_t;
				ptr->getCell()->type = ptr2->getCell()->type = 2;
				insertEdgeGraph(id_1,id_2);
				insertEdgeGraph(id_2,id_1);
			}
		}
		// Adicionar uma gap jucntion no final da fibra 5 com a fibra 1, para gerar uma corrente de reentrada
		if (ptr->getId() == start_fiber[1]->getId()-10)
		{
			f_id = fib_in_bundle-1;
			ptr2 = start_fiber[f_id];
			while (ptr2 != NULL && ptr2->getId() != start_fiber[f_id]->getId() + id_1 - 1)
				ptr2 = ptr2->getNext();
			if (ptr2 != NULL)
			{
				id_2 = ptr2->getId();
				cout << "Gap Junction between " << id_1 << " -- " << id_2 << endl;
				// Alterar a condutividade dos volumes de controle envolvidos
				ptr->getCell()->sigma = ptr->getCell()->sigmaY = sigma_c;
				ptr2->getCell()->sigma = ptr->getCell()->sigmaY = sigma_c;
				ptr->getCell()->type = ptr2->getCell()->type = 2;
				insertEdgeGraph(ptr->getId(),ptr2->getId());
				insertEdgeGraph(ptr2->getId(),ptr->getId());
			}
		}
		ptr = ptr->getNext();
	}
}

double* Graph::initializeRandomGenerator ()
{
	int i;
	double *aleatorios = new double[1000];
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,27.8/2.0);
	//cout << "==============================================" << endl;
	for (i = 0; i < 1000; i++)
	{
		aleatorios[i] = distribution(generator);
		//cout << aleatorios[i] << endl;
	}
	return (aleatorios);
	//cout << "\n==============================================" << endl << endl;
}

void Graph::setPlotPointers (int *plot)
{
	// Get the id of the initial and final cell on the central fiber of the bundle
	plot[0] = (start_fiber[1]->getId() - start_fiber[0]->getId()) / 2 - 80;
	plot[1] = (start_fiber[1]->getId() - start_fiber[0]->getId()) / 2 + 80;
}

void Graph::makeRoot (double Dx)
{
	int i, j, k;
	int num_vol;
	double x;
	double fiber_length, cell_length;
	double fl;
	// Inicializa o vetor de números aleatórios
	double *aleatorios = initializeRandomGenerator();
	Node *ptr;
	// Verifica o número de fibras no bundle e constrói um vetor com pontos randômicos
	Point *roots = setRootFibers();
	start_fiber = new Node*[fib_in_bundle];
	c_nodes = 1;
	c_segments = 0;
	//fiber_length = (rand() % 5 + 5)*scale_size; 	// Tamanho de fibra variável
	fiber_length = scale_size*10;										// Tamanho de fibra fixo
	// Montar cada fibra
	for (i = 0; i < fib_in_bundle; i++)
	{
		x = 0.0;
		fl = fiber_length;
		// Cria um primeiro volume de controle da fibra (será sempre de estímulo)
		if (i == 0)
			ptr = new Node(c_nodes,0.0,0.0,0.0,1,false,true);
		else
			ptr = new Node(c_nodes,0,roots[i-1].y,roots[i-1].z,i+1,false,true);
		insertNodeGraph(ptr);
		// Salva o ponteiro para o início da fibra
		start_fiber[i] = ptr;
		// Retira o pedaço do volume de controle do tamanho da fibra
		fl -= Dx;
		// Enquanto tiver espaço na fibra vamos preenchendo com células
		while (fl > 0.0)
		{
			// Avança um volume de controle
			x += Dx;
			// Sorteia tamanho da célula com base no vetor de aleatorios
			k = rand() % 1000;

			// Tamanho da célula
			//cell_length = 120.9 + aleatorios[k];
			cell_length = 100.0;

			// Calcula número de volumes de controle para a célula
			num_vol = nearbyint(cell_length/Dx);
			if (i == 0)
			{
				for (j = 0; j < num_vol; j++)
				{
					ptr = new Node(c_nodes,x,0.0,0.0,1,false,false);
					insertNodeGraph(ptr);
					insertEdgeGraph(c_nodes-2,c_nodes-1);
					insertEdgeGraph(c_nodes-1,c_nodes-2);
					x += Dx;
					//fl -= Dx;
				}
				// Ultimo volume eh gap junction plicate (de ligacao)
				ptr = new Node(c_nodes,x,0.0,0.0,1,true,false);
				insertNodeGraph(ptr);
				insertEdgeGraph(c_nodes-2,c_nodes-1);
				insertEdgeGraph(c_nodes-1,c_nodes-2);
				x += Dx;
				fl -= Dx*(num_vol+1);
			}
			else
			{
				for (j = 0; j < num_vol; j++)
				{
					ptr = new Node(c_nodes,x,roots[i-1].y,roots[i-1].z,i+1,false,false);
					insertNodeGraph(ptr);
					insertEdgeGraph(c_nodes-2,c_nodes-1);
					insertEdgeGraph(c_nodes-1,c_nodes-2);
					x += Dx;
					//fl -= Dx;
				}
				ptr = new Node(c_nodes,x,roots[i-1].y,roots[i-1].z,i+1,true,false);
				insertNodeGraph(ptr);
				insertEdgeGraph(c_nodes-2,c_nodes-1);
				insertEdgeGraph(c_nodes-1,c_nodes-2);
				x += Dx;
				fl -= Dx*(num_vol+1);
			}
		}
	}
	Xmax = x;
}

void Graph::printGraph (int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/graph.txt",exec_number);
	FILE *out = fopen(filename,"w+");
	Node *ptr = getListNodes();
	Edge *ptrl;
	fprintf(out,"================== Printing graph ======================\n");
	fprintf(out,"|| Id x y z Number_edges Fiber_id [Type Conductance Retropropagation? Stimulus?] || --> || Id Size x y z ||\n\n");
	while (ptr != NULL)
	{
		fprintf(out,"|| %d %.2f %.2f %.2f E%d F%d [T %d Sigma %.2f Re? %d St? %d] || ",ptr->getId(),ptr->getX(),ptr->getY(),ptr->getZ(),ptr->getNumEdges(),ptr->getIdFiber(),ptr->getCell()->type,ptr->getCell()->sigma, ptr->getRetro(), ptr->getStimulus());
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(out," --> || %d %.2f %.2f %.2f %.2f ||",ptrl->getId(),ptrl->getSize(),ptrl->getX(),ptrl->getY(),ptrl->getZ());
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
		fprintf(out,"\n");
	}
	fclose(out);
}

// MUDAR ESSA FUNÇÃO PARA UTILIZAR A DERIVADA MAXIMA
void Graph::calculateVelocity (double t, Node *ptr1, Node *ptr2, double Dt, double Dx, double *maxV, int exec_number)
{
	static double t1;
	static bool calcVelocity;
	double delta_t, delta_x, velocity;
	// Achar o pico do potencial de ação e marcar o tempo em que ele ocorreu no volume 1
	if (!flag)
	{
		if (ptr1->getCell()->y_new[0] != *maxV )
		{
			if (ptr1->getCell()->y_new[0] > *maxV)
			{
				*maxV = ptr1->getCell()->y_new[0];
				noChange = 0;
			}
			else
				noChange++;
			if (noChange > 10)
			{
				flag = true;
				calcVelocity = false;
				t1 = t-(10*Dt);
			}
		}
	}
	// Acha no volume 2 quando que o pico do potencial de ação aconteceu e realiza o cálculo da velocidade de propagação
	else
	{
		//cout << fabs(*maxV-ptr2->getCell()->y_new[0]) << endl;
		// Se a diferenca entre o valor do potencial do volume 1 com o volume 2 for menor que uma tolerancia de 1.5, calcular a velocidade
		if (!calcVelocity && fabs(*maxV-ptr2->getCell()->y_new[0]) < 10.0)
		{
			delta_t = t-t1;
			delta_x = (ptr2->getId()-ptr1->getId())*Dx;
			velocity = delta_x / delta_t * 1.0e-04;
			writeVelocity(velocity,delta_x,delta_t,t1,t,maxV,exec_number);
			calcVelocity = true;
			cout << "Volume 1 = " << ptr1->getId() << endl;
			cout << "Volume 2 = " << ptr2->getId() << endl;
		}
	}
}

Graph::Graph (double *y0, int n_eq, double Dx, int exec_number)
{
	y_0 = y0;
	num_eq = n_eq;
	flag = false;
	makeRoot(Dx);
	setGapJunctionFibers();
	//setGapJunctionFibers2();
	total_nodes = c_nodes - 1;
	total_edges = c_segments;
	printGraph(exec_number);
	cout << "[+] Grafo criado com sucesso!" << endl;
	//printCellsVTK(this,exec_number);
}
