#include "plot.h"

void printCellsVTK (Graph *g, int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/celulas.vtk",exec_number);
	FILE *fileVTK = fopen(filename,"w+");
	Node *ptr = g->getListNodes();
	Edge *ptrl;
	fprintf(fileVTK,"# vtk DataFile Version 3.0\n");
	fprintf(fileVTK,"Monodomain\n");
	fprintf(fileVTK,"ASCII\n");
	fprintf(fileVTK,"DATASET POLYDATA\n");
	fprintf(fileVTK,"POINTS %d float\n",g->getTotalNodes());
	// Imprime os pontos
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e %e %e\n",ptr->getX(),ptr->getY(),ptr->getZ());
		ptr = ptr->getNext();
	}
	fprintf(fileVTK,"LINES %d %d\n",g->getTotalEdges(),g->getTotalEdges()*3);
	ptr = g->getListNodes();
	// Imprime as linhas
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(fileVTK,"2 %d %d\n",ptr->getId()-1,ptrl->getId()-1);
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
	}
	// Imprime o valor da solução da equação do monodomínio
	fprintf(fileVTK,"POINT_DATA %d\n",g->getTotalNodes());
	fprintf(fileVTK,"SCALARS vm float 1\n");
	fprintf(fileVTK,"LOOKUP_TABLE default\n");
	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getCell()->type == 1 || ptr->getCell()->type == 2)
			fprintf(fileVTK,"%e\n",1.0);
		else
			fprintf(fileVTK,"%e\n",0.0);
		ptr = ptr->getNext();
	}
	fclose(fileVTK);
}

// Imprime informações do modelo
void printInfoModel (Graph *g, double *y0, double Dx, double Dt, double t_max, int exec_number, int num_fibers, double probGap, 
			double diameter_bundle, double s_c, double s_p, double s_t, int pac, int ret, int cool)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/info_simulation",exec_number);
	FILE *out = fopen(filename,"a");
	fprintf(out,"=================== PARAMETERS OF THE MODEL ===================\n");
	fprintf(out,"Number of segments = %d\n",g->getTotalEdges());
	fprintf(out,"Number of nodes = %d\n",g->getTotalNodes());
	fprintf(out,"V_stim = %.2lf mV\n",y0[0]);
	fprintf(out,"t_max = %.2lf s\n",t_max);
	fprintf(out,"N = %d\n",(int)(g->getXmax()/Dx));
	fprintf(out,"M = %d\n",(int)(t_max/Dt));
	fprintf(out,"Dx = %e um\n",Dx);
	fprintf(out,"Dt = %e s\n",Dt);
	fprintf(out,"Fiber length = %.2f um\n",g->getTotalNodes()*Dx/fib_in_bundle);
	fprintf(out,"Number of fibers = %d\n",num_fibers);
	fprintf(out,"Probability of gap junction = %.2f %%\n",probGap*100.0);
	fprintf(out,"Diameter of the bundle = %.2f um\n",diameter_bundle);
	fprintf(out,"\n===== Conductivities: =====\n");
	fprintf(out,"Citoplasmatic gap junction = %.2f uS/um\n",s_c);
	fprintf(out,"Plicate gap junction = %.2f uS/um\n",s_p);
	fprintf(out,"Transversal gap junction = %.2f uS/um\n",s_t);
	fprintf(out,"\n===== FLAGS: ======\n");
	fprintf(out,"Pacing = %s\n",pac?"True":"False");
	fprintf(out,"Retroprogation = %s\n",ret?"True":"False");
	fprintf(out,"Cooldown = %s\n",cool?"True":"False");
	fclose(out);
}

// Escreve um arquivo de dados da solução da equação para um vetor de volumes de controle
void writeGraphic (Graph *g, double t, const int *id, int size, int num_eq, int exec_number)
{
	int k;
	char filename[50];
	FILE *file;
	for (k = 0; k < size; k++)
	{
		sprintf(filename,"../Runs/Run%d/data%d.dat",exec_number,id[k]);
		file = fopen(filename,"a");
	
		Node *ptr = g->searchNode(id[k]);
		fprintf(file,"%e ",t);
		for (int i = 0; i < num_eq; i++)
			fprintf(file,"%e ", ptr->getCell()->y[i]);
		fprintf(file,"\n");
		fclose(file);
	}
}

// Escreve arquivo VTK com as informações do modelo no instante de tempo j
void writeVTK (Graph *g, int j, int exec_number)
{
	char nameFile[50];
	FILE *fileVTK;
	Node *ptr = g->getListNodes();
	Edge *ptrl;
	sprintf(nameFile,"../Runs/Run%d/VTK/monodomain%d.vtk",exec_number,j);
	fileVTK = fopen(nameFile,"w+");
	fprintf(fileVTK,"# vtk DataFile Version 3.0\n");
	fprintf(fileVTK,"Monodomain\n");
	fprintf(fileVTK,"ASCII\n");
	fprintf(fileVTK,"DATASET POLYDATA\n");
	fprintf(fileVTK,"POINTS %d float\n",g->getTotalNodes());
	// Imprime os pontos
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e %e %e\n",ptr->getX(),ptr->getY(),ptr->getZ());
		ptr = ptr->getNext();
	}
	fprintf(fileVTK,"LINES %d %d\n",g->getTotalEdges(),g->getTotalEdges()*3);
	ptr = g->getListNodes();
	// Imprime as linhas
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(fileVTK,"2 %d %d\n",ptr->getId()-1,ptrl->getId()-1);
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
	}
	// Imprime o valor da solução da equação do monodomínio
	fprintf(fileVTK,"POINT_DATA %d\n",g->getTotalNodes());
	fprintf(fileVTK,"SCALARS vm float 1\n");
	fprintf(fileVTK,"LOOKUP_TABLE default\n");
	ptr = g->getListNodes();
		
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e\n",ptr->getCell()->y[0]);
		ptr = ptr->getNext();
	}			
	fclose(fileVTK);
}

// Procedimento para imprimir em VTK a estrutura das fibras destacando os nós em que ocorreu retropropagação
void printCellsRetroprogVTK (Graph *g, int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/retropropagation.vtk",exec_number);
	FILE *fileVTK = fopen(filename,"w+");
	Node *ptr = g->getListNodes();
	Edge *ptrl;
	fprintf(fileVTK,"# vtk DataFile Version 3.0\n");
	fprintf(fileVTK,"Monodomain\n");
	fprintf(fileVTK,"ASCII\n");
	fprintf(fileVTK,"DATASET POLYDATA\n");
	fprintf(fileVTK,"POINTS %d float\n",g->getTotalNodes());
	// Imprime os pontos
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e %e %e\n",ptr->getX(),ptr->getY(),ptr->getZ());
		ptr = ptr->getNext();
	}
	fprintf(fileVTK,"LINES %d %d\n",g->getTotalEdges(),g->getTotalEdges()*3);
	ptr = g->getListNodes();
	// Imprime as linhas
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(fileVTK,"2 %d %d\n",ptr->getId()-1,ptrl->getId()-1);
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
	}
	// Imprime o valor da solução da equação do monodomínio
	fprintf(fileVTK,"POINT_DATA %d\n",g->getTotalNodes());
	fprintf(fileVTK,"SCALARS vm float 1\n");
	fprintf(fileVTK,"LOOKUP_TABLE default\n");
	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getRetro())
			fprintf(fileVTK,"%e\n",1.0);
		else
			fprintf(fileVTK,"%e\n",0.0);
		ptr = ptr->getNext();
	}
	fclose(fileVTK);
	//FILE *file = fopen("retro.txt","a");
	//fprintf(file,"RETROPROPAGA!\n");
	//fclose(file);
}

// Plota os gráficos dos volumes passados como parâmetros
void makePlot (const int *plot, int size, int exec_number)
{
	FILE *arq;
	char filename[30];
	int k;
	for (k = 0; k < size; k++)
	{
		sprintf(filename,"../Runs/Run%d/data%d.dat",exec_number,plot[k]);
		// Imprime o gráfico do potencial transmembrânico (V)
		arq = fopen("graph.plt","w+");
		fprintf(arq,"set grid\n");
		fprintf(arq,"set terminal png\n");
		fprintf(arq,"set output \"../Runs/Run%d/Graphics/dvdt%d.png\"\n",exec_number,k);
		fprintf(arq,"set title \"Vol = %d\"\n",plot[k]);
		fprintf(arq,"plot \"%s\" using 1:2 title \"Vm\" w l\n",filename);
		fclose(arq);
		if (!system("gnuplot graph.plt"))
			cout << "[+] Grafico plotado com sucesso!" << endl;
		else
			cout << "[-] ERRO! Plotando grafico!" << endl;
		// Imprime o gráfico da variável gate (m)
		arq = fopen("graph.plt","w+");
		fprintf(arq,"set grid\n");
		fprintf(arq,"set terminal png\n");
		fprintf(arq,"set output \"../Runs/Run%d/Graphics/dmdt%d.png\"\n",exec_number,k);
		fprintf(arq,"set title \"Vol = %d\"\n",plot[k]);
		fprintf(arq,"plot \"%s\" using 1:3 title \"m\" w l\n",filename);
		fclose(arq);
		if (!system("gnuplot graph.plt"))
			cout << "[+] Grafico plotado com sucesso!" << endl;
		else
			cout << "[-] ERRO! Plotando grafico!" << endl;
		// Imprime o gráfico da variável gate (h)
		arq = fopen("graph.plt","w+");
		fprintf(arq,"set grid\n");
		fprintf(arq,"set terminal png\n");
		fprintf(arq,"set output \"../Runs/Run%d/Graphics/dhdt%d.png\"\n",exec_number,k);
		fprintf(arq,"set title \"Vol = %d\"\n",plot[k]);
		fprintf(arq,"plot \"%s\" using 1:4 title \"m\" w l\n",filename);
		fclose(arq);
		if (!system("gnuplot graph.plt"))
			cout << "[+] Grafico plotado com sucesso!" << endl;
		else
			cout << "[-] ERRO! Plotando grafico!" << endl;
		// Imprime o gráfico da variável gate (n)
		arq = fopen("graph.plt","w+");
		fprintf(arq,"set grid\n");
		fprintf(arq,"set terminal png\n");
		fprintf(arq,"set output \"../Runs/Run%d/Graphics/dndt%d.png\"\n",exec_number,k);
		fprintf(arq,"set title \"Vol = %d\"\n",plot[k]);
		fprintf(arq,"plot \"%s\" using 1:4 title \"m\" w l\n",filename);
		fclose(arq);
		if (!system("gnuplot graph.plt"))
			cout << "[+] Grafico plotado com sucesso!" << endl;
		else
			cout << "[-] ERRO! Plotando grafico!" << endl;	
	}
}

void writeVelocity (double velocity, double delta_x, double delta_t, double t1, double t2, double *maxV, int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/info_simulation",exec_number);
	FILE *file = fopen(filename,"a");
	fprintf(file,"\n==========================================================\n");
	fprintf(file,"Propagation velocity = %.2f cm/s\n",velocity);
	fprintf(file,"delta_x = %.2f\tdelta_t = %.2f\n",delta_x,delta_t);
	fprintf(file,"t1 = %.2f\tt2 = %.2f\n",t1,t2);
	fprintf(file,"maxV = %.2f\n",*maxV);
	fprintf(file,"==========================================================\n");
	fclose(file);
}

void writeCylinderVTK (Graph *g, int exec_number)
{
	int i, k;
	double PI = 3.14159265359, delta, teta;
	delta = 0.1;
	double R = (D/2)*4.0;
	k = floor(2*PI/delta);
	cout << "Intervalos = " << k << endl;
	char filename[50];
	sprintf(filename,"../Runs/Run%d/cylinder.vtk",exec_number);
	FILE *file = fopen(filename,"w+");
	fprintf(file,"# vtk DataFile Version 3.0\n");
	fprintf(file,"Cylinder\n");
	fprintf(file,"ASCII\n");
	fprintf(file,"DATASET POLYDATA\n");
	fprintf(file,"POINTS %d float\n",2*k);	
	for (i = 0; i < k; i++)
	{
		teta = i*delta;
		fprintf(file,"%e %e %e\n",0.0,R*sin(teta),R*cos(teta));
	}
	for (i = 0; i < k; i++)
	{
		teta = i*delta;
		fprintf(file,"%e %e %e\n",g->getXmax(),R*sin(teta),R*cos(teta));
	}
	fprintf(file,"LINES %d %d\n",2*(k-1)+2,6*(k-1)+6);
	for (i = 0; i < k-1; i++)
	{
		if (i == k-2)
			fprintf(file,"2 %d %d\n",i+1,1);
		else	
			fprintf(file,"2 %d %d\n",i+1,i+2);
	}
	for (i = k-1; i < 2*(k-1); i++)
	{
		if (i == 2*(k-1)-1)
			fprintf(file,"2 %d %d\n",i+1,k+1);
		else	
			fprintf(file,"2 %d %d\n",i+1,i+2);
	}
	fprintf(file,"2 %d %d\n",1,k);
	fprintf(file,"2 %d %d\n",k/2,3*k/2);
	fclose(file);
}

void printCooldown (Graph *g, int exec_number)
{
	char filename[50];
	sprintf(filename,"../Runs/Run%d/info_simulation",exec_number);
	FILE *file = fopen(filename,"a");
	Node *ptr = g->getListNodes();		
	while (ptr != NULL)
	{
		if (ptr->getStimulus())
			fprintf(file,"\nCell %d - Cooldown = %d",ptr->getId(),ptr->getCell()->cooldown);
		ptr = ptr->getNext(); 
	}
	fclose(file);	
}

void printExecution (int i, bool retro)
{
	FILE *file;
	char filename[10];
	sprintf(filename,"../Runs/Run%d/run%d.txt",i);
	file = fopen(filename,"w+");
	if (retro)
		fprintf(file,"[!] RETROPROPAGATION happen!\n");
	else
		fprintf(file,"[+] Normal activity.\n");
	fclose(file);
}
