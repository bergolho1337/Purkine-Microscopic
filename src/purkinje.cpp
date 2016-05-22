#include "purkinje.h"
#include "plot.h"

Node::Node () {}
Edge::Edge () {}
Point::Point () {}

int noChange;						// Contador da velocidade de propagação
int c_nodes;						// Contador de nós
int c_segments;					// Contador de segmentos
double *y_0;						// Vetor das condições inicias
int num_eq;							// Número de equações do modelo celular
Node **start_fiber;			// Vetor que armazena os volumes inicias de cada fibra
bool descendo, subindo;	// Flags para velocidade de propagacao
bool flag;
bool calcVelocity;


// Construtor de um Node
Node::Node (int id, double x, double y, double z, int f_id, bool gap, bool stim)
{
	this->id = id;
	this->point = new Point(x,y,z);
	this->fiber_id = f_id;
	this->gapJ = gap;
	this->stimulus = stim;
	this->num_edges = 0;
	this->retroProg = false;
	this->next = NULL;
	this->edges = NULL;
	// Aloca os vetores da solução da equacao do monodominio
	this->vol = new Volume();
	this->vol->y_old = new double[num_eq];
	this->vol->y_new = new double[num_eq];
	retroProg = false;
	if (stim)
	{
		for (int i = 0; i < num_eq; i++)
			this->vol->y_old[i] = y_0[i];
	}
	else
	{
		this->vol->y_old[0] = y_0[num_eq];
		for (int i = 1; i < num_eq; i++)
			this->vol->y_old[i] = y_0[i];
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
	edge->setSize(calcNorm(ptr1->getX(),ptr1->getY(),ptr1->getZ(),ptr2->getX(),ptr2->getY(),ptr2->getZ()));
	edge->setNext(NULL);
	edge->setDest(ptr2);
	// Inicialmente todas arestas são ligações citoplasmaticas
	edge->setSigma(sigma_c);

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
		if (p[i].getY() == y && p[i].getZ() == z)
			return (true);
	}
	return (false);
}

// Seta os pontos que irão ser a raiz do Bundle
Point* Graph::setRootFibers ()
{
	int i, k;
	double y, z, norm;
	Point *point = new Point[fib_in_bundle-1]();
	srand(time(NULL));
	for (i = 0, k = 0; i < fib_in_bundle-1; i++)
	{
		do
		{
			y = makeRandom();
			z = makeRandom();
			norm = sqrt(y*y + z*z);
		}while (alreadyMake(point,y,z,k) && norm >= D);
		point[i].setX(0);
		point[i].setY(y);
		point[i].setZ(z);
		k++;
	}
	return (point);
}

// Função para setar as gap junctions longitudinais (na mesma fibra)
void Graph::setLongitidinalGapJunction (double Dx)
{
	int num_vol, cont;
	Node *ptr;
	Edge *ptrl;
	num_vol = nearbyint(cell_length/Dx);
	cont = 0;
	ptr = getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getStimulus())
			cont = 0;
		if (cont == num_vol)
		{
			cont = 1;
			ptrl = ptr->getEdges();
			while (ptrl->getNext() != NULL)
				ptrl = ptrl->getNext();
			ptrl->setSigma(G_p);
			if (ptr->getNext() != NULL && ptr->getNext()->getIdFiber() == ptr->getIdFiber())
			{
				ptr = ptr->getNext();
				ptrl = ptr->getEdges();
				ptrl->setSigma(G_p);
			}
		}
		cont++;
		ptr = ptr->getNext();
	}
	//cout << "[+] setLongitidinalGapJunction is ok!" << endl;
}

// Função para setar as gap junctions transversais. Aqui toda fibra central pode possuir gapJ transversais.
void Graph::setTransversalGapJunction (double Dx)
{
	Node *ptr1, *ptr2;
	Edge *ptrl;
	int id_1, id_2, f_id;
	double r;
	// Pulamos o primeiro bundle
	ptr1 = getListNodes()->getNext();
	// Percorrer somente a fibra do meio
	while (ptr1->getIdFiber() == 1)
	{
		r = (double)rand()/(double)RAND_MAX;
		if (r <= prob_gapJ)
		{
			id_1 = ptr1->getId();
			// Qual fibra vou me ligar ?
			do {
				f_id = rand() % fib_in_bundle;
			}while(f_id == 0);
			// Acha o volume que vou me ligar
			ptr2 = start_fiber[f_id];
			while (ptr2 != NULL && ptr2->getId() != start_fiber[f_id]->getId() + id_1 - 1)
				ptr2 = ptr2->getNext();
			if (ptr2 != NULL)
			{
				id_2 = ptr2->getId();
				insertEdgeGraph(id_1,id_2);
				insertEdgeGraph(id_2,id_1);
				// Inicializa o sigma da ligação
				ptrl = ptr1->getEdges();
				while (ptrl->getNext() != NULL)
					ptrl = ptrl->getNext();
				ptrl->setSigma(G_t);
				ptrl = ptr2->getEdges();
				while (ptrl->getNext() != NULL)
					ptrl = ptrl->getNext();
				ptrl->setSigma(G_t);
			}
		}
		ptr1 = ptr1->getNext();
	}
	//cout << "[+] setTransversalGapJunction is ok!" << endl;
}

// Atribui os volumes que possuem gapJ transversais. Aqui so tentamos colocar a partir da metade da fibra para frente
void Graph::setTransversalGapJunction2 (double Dx)
{
	Node *ptr1, *ptr2;
	Edge *ptrl;
	int id_1, id_2, f_id;
	double r;
	ptr1 = getListNodes();
	// Pulamos até chegar no meio da fibra
	while (ptr1->getId() < start_fiber[1]->getId()/2)
		ptr1 = ptr1->getNext();
	// Percorrer somente a fibra do meio
	while (ptr1->getIdFiber() == 1)
	{
		r = (double)rand()/(double)RAND_MAX;
		if (r <= prob_gapJ)
		{
			id_1 = ptr1->getId();
			// Qual fibra vou me ligar ?
			do {
				f_id = rand() % fib_in_bundle;
			}while(f_id == 0);
			// Acha o volume que vou me ligar
			ptr2 = start_fiber[f_id];
			while (ptr2 != NULL && ptr2->getId() != start_fiber[f_id]->getId() + id_1 - 1)
				ptr2 = ptr2->getNext();
			if (ptr2 != NULL)
			{
				id_2 = ptr2->getId();
				insertEdgeGraph(id_1,id_2);
				insertEdgeGraph(id_2,id_1);
				// Inicializa o sigma da ligação
				ptrl = ptr1->getEdges();
				while (ptrl->getNext() != NULL)
					ptrl = ptrl->getNext();
				ptrl->setSigma(G_t);
				ptrl = ptr2->getEdges();
				while (ptrl->getNext() != NULL)
					ptrl = ptrl->getNext();
				ptrl->setSigma(G_t);
			}
		}
		ptr1 = ptr1->getNext();
	}
}

double* Graph::initializeRandomGenerator ()
{
	int i;
	double *aleatorios = new double[1000];
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,27.8/2.0);
	for (i = 0; i < 1000; i++)
		aleatorios[i] = distribution(generator);
	return (aleatorios);
}

// Inicializa os ponteiros para os volumes de referencia do calculo da velocidade
Node** Graph::setPlotPointers ()
{
	// Set the plot to be at the end of the central fiber
	Node **plot = new Node*[2];
	cout << "[+] Setting plot pointers to volumes:" << endl;
	plot[0] = searchNode(start_fiber[1]->getId() - 600);
	plot[1] = searchNode(start_fiber[1]->getId() - 80);
	cout << "\tVolume 1: " << plot[0]->getId() << endl;
	cout << "\tVolume 2: " << plot[1]->getId() << endl;
	return (plot);
}

// Constroi a estrutura da rede
void Graph::makeRoot (double Dx)
{
	int i, j;
	int num_vol;
	double x;
	double fl;
	// Inicializa o vetor de números aleatórios
	Node *ptr;
	// Verifica o número de fibras no bundle e constrói um vetor com pontos randômicos
	Point *roots = setRootFibers();

	start_fiber = new Node*[fib_in_bundle];
	c_nodes = 1;
	c_segments = 0;

	// Montar cada fibra
	for (i = 0; i < fib_in_bundle; i++)
	{
		x = 0.0;
		fl = fiber_length;
		// Cria um primeiro volume de controle da fibra (será sempre de estímulo)
		if (i == 0)
			ptr = new Node(c_nodes,0.0,0.0,0.0,1,false,true);
		else
			ptr = new Node(c_nodes,0,roots[i-1].getY(),roots[i-1].getZ(),i+1,false,true);
			//ptr = new Node(c_nodes,0,roots[i-1].getY(),roots[i-1].getZ(),i+1,false,false);
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

			// Calcula número de volumes de controle para a célula
			num_vol = nearbyint(cell_length/Dx);
			if (i == 0)
			{
				for (j = 0; j < num_vol-1; j++)
				{
					ptr = new Node(c_nodes,x,0.0,0.0,1,false,false);
					insertNodeGraph(ptr);
					insertEdgeGraph(c_nodes-2,c_nodes-1);
					insertEdgeGraph(c_nodes-1,c_nodes-2);
					x += Dx;
				}
				// Ultimo volume eh gap junction plicate (de ligacao)
				ptr = new Node(c_nodes,x,0.0,0.0,1,true,false);
				insertNodeGraph(ptr);
				insertEdgeGraph(c_nodes-2,c_nodes-1);
				insertEdgeGraph(c_nodes-1,c_nodes-2);
				//x += Dx;
				fl -= Dx*(num_vol);
			}
			else
			{
				for (j = 0; j < num_vol-1; j++)
				{
					ptr = new Node(c_nodes,x,roots[i-1].getY(),roots[i-1].getZ(),i+1,false,false);
					insertNodeGraph(ptr);
					insertEdgeGraph(c_nodes-2,c_nodes-1);
					insertEdgeGraph(c_nodes-1,c_nodes-2);
					x += Dx;
					//fl -= Dx;
				}
				// Insere o último volume da célula
				ptr = new Node(c_nodes,x,roots[i-1].getY(),roots[i-1].getZ(),i+1,true,false);
				insertNodeGraph(ptr);
				insertEdgeGraph(c_nodes-2,c_nodes-1);
				insertEdgeGraph(c_nodes-1,c_nodes-2);
				//x += Dx;
				fl -= Dx*(num_vol);
			}
		}
	}
	Xmax = x;
	//cout << "[+] makeRoot is ok!" << endl;
}

// Imprime o grafo em arquivo .txt
void Graph::printGraph (int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/graph.txt",exec_number);
	FILE *out = fopen(filename,"w+");
	Node *ptr = getListNodes();
	Edge *ptrl;
	fprintf(out,"================== Printing graph ======================\n");
	fprintf(out,"|| Id x y z Number_edges Fiber_id [Type Retropropagation? Stimulus?] || --> || Id x y z Sigma Size ||\n\n");
	while (ptr != NULL)
	{
		fprintf(out,"|| %d %.2f %.2f %.2f E%d F%d [Re? %d St? %d] || ",ptr->getId(),ptr->getX(),ptr->getY(),ptr->getZ(),ptr->getNumEdges(),ptr->getIdFiber(), ptr->getRetro(), ptr->getStimulus());
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(out," --> || %d %.2f %.2f %.2f %.2f %.2f ||",ptrl->getDest()->getId(),ptrl->getDest()->getX(),ptrl->getDest()->getY(),ptrl->getDest()->getZ(),ptrl->getSigma(),ptrl->getSize());
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
		fprintf(out,"\n");
	}
	fclose(out);
}

// Calcula a velocidade de propagacao do segundo potencial de acao
void Graph::calculateVelocity (double t, Node *ptr1, Node *ptr2, double Dt, double Dx, double *maxV, int exec_number)
{
	static double t1;
	double delta_t, delta_x, velocity;
	// Polarizando o primeiro potencial de ação ?
	if (ptr1->getVolume()->y_new[0]-ptr1->getVolume()->y_old[0] < -0.1)
			descendo = true;
	// 2o potencial de ação foi gerado se for verdadeiro ...
	if (descendo && ptr1->getVolume()->y_new[0]-ptr1->getVolume()->y_old[0] > 0.1)
		subindo = true;

	if (subindo)
	{
		if (!flag)
		{
			if (ptr1->getVolume()->y_new[0] != *maxV )
			{
				if (ptr1->getVolume()->y_new[0] > *maxV)
				{
					*maxV = ptr1->getVolume()->y_new[0];
					noChange = 0;
				}
				else
					noChange++;
				if (noChange > 15)
				{
					flag = true;
					calcVelocity = false;
					t1 = t-(15*Dt);
				}
			}
		}
		else
		{
			//cout << "maxV = " << *maxV << " happens at t = " << t1 << endl;
			if (!calcVelocity && fabs(*maxV-ptr2->getVolume()->y_new[0]) < 0.1)
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
}

// Construtor do grafo
Graph::Graph (double *y0, int n_eq, double Dx, int exec_number)
{
	y_0 = y0;
	num_eq = n_eq;
	descendo = subindo = false;
	flag = false;
	makeRoot(Dx);
	setLongitidinalGapJunction(Dx);
	setTransversalGapJunction(Dx);
	//setTransversalGapJunction2(Dx);
	total_nodes = c_nodes - 1;
	total_edges = c_segments;
	printGraph(exec_number);
	printCellsVTK(this,exec_number);
	cout << "[+] Graph builded with sucess!" << endl;
}
