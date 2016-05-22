#include "monodomain.h"

// Parâmetros fixados
double dt;										// Tamanho do intervalo no tempo
double dx;										// Tamanho do intervalo no espaço
double C_m;										// Valor da capacitância do modelo
double BETA;									// Relação área/volume. Considerando o volume um cilindro
double radius;								// Raio de uma celula
double **A;										// Matriz do metodo implicito
double *b;										// Vetor do metodo implicito
double *pivot;								// Vetor dos pivos (LU)

// Setando os parâmetros do modelo nas variáveis globais
void setParameters (double Dx, double Dt, double cm)
{
	dx = Dx;
	dt = Dt;
	C_m = cm;
	radius = diameter_cell / 2;
	BETA = (2/cell_length)+(2/radius);
}


// Constrói a matriz do método implícito e já gera a decomposição LU.
void makeMatrix_A (Graph *g)
{
	Node *ptr;
	Edge *ptrl;
	int i, total_nodes;
	double alfa, phi, teta, sum, sigma;

	// Inicializa variaveis
	total_nodes = g->getTotalNodes();
	phi = (sigma_c*dt)/(BETA*C_m*dx*dx);
	teta = (G_p*dt)/(BETA*C_m*PI*radius*radius*dx);
	ptr = g->getListNodes();

	// Aloca memoria
	A = new double*[total_nodes]();
	for (i = 0; i < total_nodes; i++)
		A[i] = new double[total_nodes]();
	pivot = new double[total_nodes]();

	// Percorrer os volumes
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		sum = 0;
		while (ptrl != NULL)
		{
			if (ptrl->getSigma() == sigma_c)
				sigma = phi;
			else
				sigma = teta;
			sum += sigma;
			A[ptr->getId()-1][ptrl->getDest()->getId()-1] = -sigma;
			ptrl = ptrl->getNext();
		}
		sum += 1;
		// Condição de Neumann nas folhas
		if (ptr->getNumEdges() == 1 && !ptr->getStimulus())
		{
			A[ptr->getId()-1][ptr->getId()-1] = 1;
			A[ptr->getId()-1][ptr->getEdges()->getDest()->getId()-1] = -1;
		}
		else
			A[ptr->getId()-1][ptr->getId()-1] = sum;
		ptr = ptr->getNext();
	}

	// Decomposição LU
	int j, k, p;
	double Amax, t, m, r, Mult;
	// 1 PASSO: Transformar a matriz A do problema em duas matrizes triangulares L e U.
	for (i = 0; i < total_nodes; i++)
		pivot[i] = i;
	for (j = 0; j < total_nodes-1; j++)
	{
		// Escolher pivot
		p = j;
		Amax = abs(A[j][j]);
		// Verifica na coluna a ser eliminada qual elemento possui o maior valor absoluto, este elemento será o pivô.
		for (k = j+1; k < total_nodes; k++)
		{
			if (abs(A[k][j]) > Amax)
			{
				Amax = abs(A[k][j]);
				p = k;
			}
		}
		// Se (p != j) então deve-se trocar de linhas
		if (p != j)
		{
			for (k = 0; k < total_nodes; k++)
			{
				t = A[j][k];
				A[j][k] = A[p][k];
				A[p][k] = t;
			}
			m = pivot[j];
			pivot[j] = pivot[p];
			pivot[p] = m;
		}
		if (abs(A[j][j]) != 0)
		{
			// Eliminação de Gauss
			r = 1 / A[j][j];
			for (i = j+1; i < total_nodes; i++)
			{
				Mult = A[i][j]*r;
				A[i][j] = Mult;
				for (k = j+1; k < total_nodes; k++)
					A[i][k] = A[i][k] - Mult*A[j][k];
			}
		}
	}
}

/*
bool printRetropropagation (Graph *g, int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/info_simulation",exec_number);
	FILE *out = fopen(filename,"a");
	int cont;
	int retro_fiber[fib_in_bundle];
	set<int>::iterator it;
	Node *ptr;
	double fiber_lenght;
	// Inicializa o vetor que conta o numero de volumes que tenham retropropagado nas fibras
	for (int i = 0; i < fib_in_bundle; i++)
		retro_fiber[i] = 0;
	fprintf(out,"======== Retro-propagation Control Volumes ======= [size,percentage] ======\n");
	fiber_lenght = g->getTotalNodes()*dx/fib_in_bundle;
	for (it = retropropaga.begin(); it != retropropaga.end(); ++it)
	{
		cont = 0;
		ptr = g->searchNode(*it);
		// Checar em qual fibra houve a retropropagação e incrementa o numero de volumes em que houve retropro.
		retro_fiber[ptr->getIdFiber()-1]++;
		// Printa a area em que houve a retropropagacao
		fprintf(out,"%d",*it);
		while (ptr->getRetro())
		{
			cont++;
			ptr = ptr->getNext();
		}
		fprintf(out," -> [%lf um %lf %%]\n",cont*dx,cont*dx*100.0/fiber_lenght);
	}
	// ** FALTA DESVIO PADRAO
	fprintf(out,"\n======= Average of Retropropagation Control Volumes per Fibers =============== [ %% ] ==========\n");
	for (int i = 0; i < fib_in_bundle; i++)
	{
		fprintf(out,"Fiber %d -> %lf %%\n",i+1,retro_fiber[i]*100.0/(g->getTotalNodes()/fib_in_bundle));
	}
	fprintf(out,"================================================================================================\n");
	fclose(out);
	if (retropropaga.size() > 0)
		return (true);
	else
		return (false);
}

void checkRetropropagation (Node *ptr)
{
	Node *pg;
	Edge *ptrl;
	ptrl = ptr->getEdges();
	while (ptrl != NULL)
	{
		pg = ptrl->getDest();
		// Trocar essa condição
		if (pg->getId() > ptr->getId() && pg->getCell()->y[0] > ptr->getCell()->y[0] + 5.0)
		{
			ptr->setRetro(true);
			retropropaga.insert(ptr->getId());
		}
		ptrl = ptrl->getNext();
	}
}
*/

// (Reação): Resolve as EDOs associadas para cada volume de controle da rede
void solveEDO (Graph *g, double t, int k, Func *func, int num_eq)
{
	int i;
	double f[num_eq];
	Volume *v;
	Node *ptr = g->getListNodes();
	while (ptr != NULL)
	{
		v = ptr->getVolume();
		// Verifica se a celula se polarizou novamente
		if (v->y_old[0] < -80)
			ptr->setRetro(false);
		for (i = 0; i < num_eq; i++)
		{
			// Calcular o potencial transmembranico intermediário -> V_{i/2} = V*
			if (i == 0)
			{
				f[0] = func[0](k,t,v->y_old)*dt;
				v->V_star = v->y_old[0] + f[0];
			}
			// gate -> n -> n+1
			else
			{
				f[i] = func[i](k,t,v->y_old)*dt;
				v->y_new[i] = v->y_old[i] + f[i];
			}
		}
		ptr = ptr->getNext();
	}
}

// Avança no tempo n -> n+1
void nextTimestep (Graph *g, int num_eq)
{
	int i;
	Node *ptr = g->getListNodes();
	Volume *v;
	while (ptr != NULL)
	{
		v = ptr->getVolume();
		for (i = 0; i < num_eq; i++)
			v->y_old[i] = v->y_new[i];
		ptr = ptr->getNext();
	}
}

void solveEDP_Imp (Graph *g)
{
	Node *ptr;
	// Resolve o sistema do metodo implicito
	double *Vnext = LU(g->getTotalNodes());
	// Seta os valores de V_n+1 nos volumes de controle
	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		ptr->getVolume()->y_new[0] = Vnext[ptr->getId()-1];
		ptr = ptr->getNext();
	}
	delete [] Vnext;
	delete [] b;
}

void makeVector_b (Graph *g)
{
	Node *ptr;
	b = new double[g->getTotalNodes()];
	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		// Condição de Neumann
		if (ptr->getNumEdges() == 1 && !ptr->getStimulus())
			b[ptr->getId()-1] = 0.0;
		else
			b[ptr->getId()-1] = ptr->getVolume()->V_star;
		ptr = ptr->getNext();
	}
}

// Resolve um sistema linear em que a matriz A já foi decomposta.
double* LU (int n)
{
	int i, j, k;
	double *y = new double[n];
	double soma;
	k = pivot[0];
	y[0] = b[k];
	// Realizar substituições sucessivas para resolver o sistema triangular inferior: Ly = b
	for (i = 1; i < n; i++)
	{
		soma = 0;
		for (j = 0; j <= i-1; j++)
			soma += A[i][j]*y[j];
		k = pivot[i];
		y[i] = b[k] - soma;
	}
	// Realizar substituições retroativas para resolver o sistema triangular superior: Ux = y
	double *x = new double[n];
	x[n-1] = y[n-1] / A[n-1][n-1];
	for (i = n-2; i >= 0; i--)
	{
		soma = 0;
		for (j = i+1; j < n; j++)
			soma += A[i][j]*x[j];
		x[i] = (y[i] - soma) / A[i][i];
	}
	delete [] y;
	return (x);
}
