#include "monodomain.h"

// Parâmetros fixados dos métodos
double dt;
double dx;
double C_m;
double BETA;
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
	BETA = 2/R;
	cout << "Beta = " << BETA << endl;
	cout << "Cm = " << C_m << endl; 
}

bool printRetropropagation (Graph *g, int exec_number)
{
	char filename[50];
	sprintf(filename,"Runs/Run%d/info_simulation",exec_number);
	FILE *out = fopen(filename,"a");
	int cont;
	int retro_fiber[fib_in_bundle];
	set<int>::iterator it;
	Node *ptr;
	double fiber_lenght;
	// Inicializa o vetor que conta o numero de volumes que tenham retropropagado nas fibras
	for (int i = 0; i < fib_in_bundle; i++)
		retro_fiber[i] = 0;
	fprintf(out,"\nPacing = %d ms\n",CYCLE_L/10);
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
						r = (c1->sigma*dt)/(dx*dx*C_m*BETA);
						//cout << r << endl;
						checkCFL(r);
						// Realizar pacing ?
						
						c2 = ptr->getEdges()->getDest()->getCell();
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
						sigma_med = (c1->sigma + c2->sigma)/2.0;
						r = (dt*sigma_med)/(dx*dx*C_m*BETA);
						checkCFL(r);
						// Corrente entrando no volume
						if (ptrl->getId() < ptr->getId())
							V_star += r*(c2->V_star - c1->V_star);
						// Corrente saindo do volume
						else
						{
							// Checa se uma corrente de reentrada esta prestes a acontecer, nesse caso liberar a corrente alterando a condutividade
							if (c1->y[0] < -80 && c2->y[0] > -79)
							{
								r = (dt*sigma_c)/(dx*dx*C_m*BETA);
								checkCFL(r);
							  }
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

// ** TESTE ***
void solveEDP2 (Graph *g, double t, int k, int N, int M, Func *func, int num_equacoes, double v_gate)
{
	Cell *c, *c1, *c2, *c3;
	int j, connect;
	double r;
	double f[num_equacoes];
	// Resolve os pontos da malha para um instante de tempo
	Node *ptr = g->getListNodes();
	while (ptr != NULL)
	{
		connect = ptr->getNumEdges();
		c = ptr->getCell();
		switch (connect)
		{
			// Volume inicial ou volume das folhas
			case 1:
			{
				// Se o volume for de estímulo então é o inicial, com isso basta resolver a EDO para todas variáveis
				if (ptr->getStimulus())
				{
					c1 = ptr->getCell();
					r = (c1->sigma*dt)/(dx*dx*C_m*BETA);
					if (r >= 0.5)
					{
						cout << "[-] ERROR! Problem in the CFL condition!" << endl;
						cout << "On node " << ptr->getId() << " || r = " << r << endl;
						exit(1);
					}
					else
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
								for (j = 0; j < num_equacoes; j++)
									c1->y_new[j] = c1->y[j];
							}
							// Resolver normalmente
							else
							{
								for (j = 0; j < num_equacoes; j++)
								{
									f[j] = func[j](k,t,c1->y)*dt;
									c1->y_new[j] = c1->y[j] + f[j];
								}
							}
						}
					}
				}
				// Senão é folha, e aplica-se uma condição de Neumann dv/dt = 0
				else
				{
					c1 = ptr->getCell();
					c2 = ptr->getEdges()->getCell();
					r = (c1->sigma*dt)/(dx*dx*C_m*BETA);
					if (r >= 0.5)
					{
						cout << "[-] ERROR! Problem in the CFL condition!" << endl;
						cout << "On node " << ptr->getId() << " || r = " << r << endl;
						exit(1);
					}
					else
					{
						// V* -> V^n+1
						c1->y_new[0] = r*(c2->V_star-c1->V_star) + c1->V_star;
						// gate^n -> gate^n+1
						for (j = 1; j < num_equacoes; j++)
						{
							// Ultrapassou o limite de v_gate para o volume ser estimulado ?
							if (c1->y[0] > v_gate)
							{
								f[j] = func[j](k,t,c1->y)*dt;
								c1->y_new[j] = c1->y[j] + f[j];
							}
							// Senão mantém o volume no mesmo estado
							else
								c1->y_new[j] = c1->y[j];
						}
					}
				}
				break;
			}
			// Volumes no interior da malha devem ser resolvidas por FTCP (Forward Time, Centered Space)
			case 2:
			{
				Node *ptr2, *ptr3;
				c1 = ptr->getCell();
				c2 = ptr->getEdges()->getCell();
				c3 = ptr->getEdges()->getNext()->getCell();

				ptr2 = ptr->getEdges()->getDest();
				ptr3 = ptr->getEdges()->getNext()->getDest();

				// ***** TESTE *******
				// Checa se uma corrente de reentrada esta prestes a acontecer
				if (ptr2->getId() > ptr->getId() && c->y[0] < -80 && c2->y[0] > -79 )
					r = (sigma_c*dt)/(dx*dx*C_m*BETA); 
				else if (ptr3->getId() > ptr->getId() && c->y[0] < -80 && c3->y[0] > -79)
					r = (sigma_c*dt)/(dx*dx*C_m*BETA);
				else
					r = (c->sigma*dt)/(dx*dx*C_m*BETA);

				if (r >= 0.5)
				{
					cout << "[-] ERROR! Problem in the CFL condition!" << endl;
					cout << "On node " << ptr->getId() << " || r = " << r << endl;
					exit(1);
				}
				else
				{
					// *** NESSA PARTE DEVE-SE CHECAR A RETROPROPRAGAÇÃO
					checkRetropropagation(ptr);
					// V* -> V^n+1
					c1->y_new[0] = r*(c2->V_star - 2*c1->V_star + c3->V_star) + c1->V_star;
					// gate^n -> gate^n+1
					for (j = 1; j < num_equacoes; j++)
					{
						// Ultrapassou o limite de v_gate para o volume ser estimulado ?
						if (c1->y[0] > v_gate)
						{
							f[j] = func[j](k,t,c1->y)*dt;
							c1->y_new[j] = c1->y[j] + f[j];
						}
						// Senão mantém o volume no mesmo estado
						else
							c1->y_new[j] = c1->y[j];
					}
				}
				break;	
			}
			// Senão é célula de bifurcação
			default:
			{
				Edge *ptrl;
				int j;
				double Vstar_viz;
				double Vstar_x, Vstar_y, rx, ry;
				Vstar_x = Vstar_y = 0.0; Vstar_viz = 0.0;
				// Aloca um vetor de ponteiros para cada vizinho
				Node **pg = new Node*[connect];
				c1 = ptr->getCell();
				ptrl = ptr->getEdges();
				// Verifica se uma corrente de reentrada esta prestes a ocorrer
				while (ptrl != NULL)
				{
					c = ptrl->getDest()->getCell();
					// Celula que estou nao esta estimulada, porem a da frente esta, nesse caso nao deve ocorrer bloqueio
					if (ptrl->getId() > ptr->getId() && c1->y[0] < -80 && c->y[0] > -79 )
					{
						rx = ry = (sigma_c*dt)/(dx*dx*C_m*BETA);
						ptr->setRetro(true);
						//cout << "Reentry corrent at " << ptr->getId() << endl;
						break;
					}
					// Se uma das minhas vizinhas estiver retropropagando deixo o estimulo passar, pois o bloqueio ocorre em uma direcao
					else if (ptrl->getDest()->getRetro())
					{
						rx = ry = (sigma_c*dt)/(dx*dx*C_m*BETA);
						ptr->setRetro(true);
						break;
					}
					else
						rx = ry = (c1->sigma*dt)/(dx*dx*C_m*BETA);

					ptrl = ptrl->getNext();
				}

				ptrl = ptr->getEdges();

				if (rx >= 0.5 || ry >= 0.5)
				{
					cout << "[-] ERROR! Problem in the CFL condition!" << endl;
					cout << "On node " << ptr->getId() << " || rx = " << rx << " and ry = " << ry << endl;
					exit(1);
				}
				else
				{
					// Setar os ponteiros para os nós das ligações
					j = 0;
					//cout << "---> Solving cell " << ptr->getId() << endl;
					while (ptrl != NULL)
					{
						pg[j] = ptrl->getDest();
						ptrl = ptrl->getNext();
						j++;
					}
					// Somar todas as contribuições de V_star
					for (j = 0; j < connect; j++)
					{	
						// No eixo x, estão na mesma fibra
						if (pg[j]->getIdFiber() == ptr->getIdFiber())
						{
							c2 = pg[j]->getCell();
							Vstar_x += c2->V_star;
						}
						// No eixo y, estão em fibras diferentes
						else
						{
							c2 = pg[j]->getCell();
							Vstar_y += c2->V_star;
						}
					}
					Vstar_viz = Vstar_x + Vstar_y;
					// V* -> V^n+1
					c1->y_new[0] = rx*(Vstar_viz - connect*c1->V_star) + c1->V_star;
					// gate^n -> gate^n+1
					for (j = 1; j < num_equacoes; j++)
					{
						// Ultrapassou o limite de v_gate para o volume ser estimulado ?
						if (c1->y[0] > v_gate)
						{
							f[j] = func[j](k,t,c1->y)*dt;
							c1->y_new[j] = c1->y[j] + f[j];
						}
						// Senão mantém o volume no mesmo estado
						else
							c1->y_new[j] = c1->y[j];
					}
				}
				//delete [] *pg;
				break;
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
