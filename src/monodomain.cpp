#include "monodomain.h"

// Parâmetros fixados dos métodos
double dt;
double dx;
double C_m;
double BETA;									// Considerado como um cubo de dimensão (dx,dx,dx)
bool pacing;
bool retrop;
bool cooldown;
set<int> retropropaga;				// Conjunto que armazena os nós de retropropagação

void setParameters (double Dx, double Dt, double cm, double R, int PACING, int RETROP, int COOLDOWN)
{
	dx = Dx;
	dt = Dt;
	pacing = PACING;
	retrop = RETROP;
	cooldown = COOLDOWN;
	C_m = cm;
	BETA = 4/dx;
	cout << "Beta = " << BETA << endl;
	cout << "Cm = " << C_m << endl;
}

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

void solveEDO (Graph *g, double t, int k, Func *func, int num_eq)
{
	int i;
	double f[num_eq];
	Cell *c;
	Node *ptr = g->getListNodes();
	while (ptr != NULL)
	{
		c = ptr->getCell();
		// Verifica se a celula se polarizou novamente
		if (c->y[0] < -80)
			ptr->setRetro(false);
		for (i = 0; i < num_eq; i++)
		{
			// Calcular o potencial transmembranico intermediário -> V_{i/2} = V*
			if (i == 0)
			{
				f[0] = func[0](k,t,c->y)*dt;
				c->V_star = c->y[0] + f[0];
			}
			// gate -> n -> n+1
			else
			{
				f[i] = func[i](k,t,c->y)*dt;
				c->y_new[i] = c->y[i] + f[i];
			}
		}
		ptr = ptr->getNext();
	}
}

void solveEDP (Graph *g, Func *func, double t, int k, int num_eq)
{
	Node *ptr;
	Cell *c1, *c2;
	double r;
	double f[num_eq];
	int j;

	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		c1 = ptr->getCell();
		// Verifica quantas arestas estao ligados ao volume de controle atual
		switch (ptr->getNumEdges())
		{
			case 1: {
					// Se for celula de estimulo
					if (ptr->getStimulus())
					{
						// Realizar pacing ?
						if (pacing && k % CYCLE_L == 0 && k != 0)
							c1->y_new[0] = -70;
						// Senão resolver a EDO normalmente para a V
						else
						{
							// O volume está em cooldown ?
							if (cooldown && k < c1->cooldown)
							{
								for (j = 0; j < num_eq; j++)
									c1->y_new[j] = c1->y[j];
							}
							// Resolver normalmente
							else
							{
								for (j = 0; j < num_eq; j++)
								{
									f[j] = func[j](k,t,c1->y)*dt;
									c1->y_new[j] = c1->y[j] + f[j];
								}
							}
						}
					}
					// Senao a celula eh folha, logo aplica-se condicao de Neumann dv/dt=0
					else
					{
						c2 = ptr->getEdges()->getDest()->getCell();
						// Citoplasma ? Resolver usando sigma_c
						if (c2->type == 0)
							r = (c2->sigma*dt)/(dx*dx*C_m*BETA);
						// Senão é gap junctions: resolver usando G_p
						else
							r = (c2->sigma*dt)/(dx*dx*dx*C_m*BETA);
						//cout << r << endl;
						checkCFL(r);

						c1->y_new[0] = r*(c2->V_star - c1->V_star) + c1->V_star;
					}
					break;
				}
			default: {
					// Celula de bifurcacao caso geral
					double V_star = 0.0;
					double sigma_med;
					Edge *ptrl = ptr->getEdges();
					c1 = ptr->getCell();
					while (ptrl != NULL)
					{
						c2 = ptrl->getDest()->getCell();
						//sigma_med = (c1->sigma + c2->sigma)/2.0;
						//r = (dt*sigma_med)/(dx*dx*C_m*BETA);
						// Citoplasma ? Resolver usando sigma_c
						if (c2->type == 0)
							r = (c2->sigma*dt)/(dx*dx*C_m*BETA);
						// Senão é gap junctions: resolver usando G_p
						else
							r = (c2->sigma*dt)/(dx*dx*dx*C_m*BETA);
						checkCFL(r);
						// Corrente entrando no volume
						if (ptrl->getId() < ptr->getId())
							V_star += r*(c2->V_star - c1->V_star);
						// Corrente saindo do volume
						else
						{
							// Checa se uma corrente de reentrada esta prestes a acontecer, nesse caso liberar a corrente alterando a condutividade
							//if (c1->y[0] < -80 && c2->y[0] > -79)
							//{
							//	r = (dt*sigma_c)/(dx*dx*C_m*BETA);
							//	checkCFL(r);
							//}
							V_star -= r*(c1->V_star - c2->V_star);
						}
						ptrl = ptrl->getNext();
					}
					checkRetropropagation(ptr);
					c1->y_new[0] = V_star + c1->V_star;
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
	Cell *c;
	while (ptr != NULL)
	{
		c = ptr->getCell();
		for (i = 0; i < num_eq; i++)
			c->y[i] = c->y_new[i];
		ptr = ptr->getNext();
	}
}

void checkCFL (double r)
{
	if (r >= 0.5)
	{
		cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
		exit(1);
	}
}
