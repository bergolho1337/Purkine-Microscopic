#include <iostream>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include "purkinje.h"
#include "plot.h"
#include "monodomain.h"

// Deve-se incluir um modelo celular
#include "noble.h"
//#include "mitchell.h"

using namespace std;

#define PACING 0			// Flag para realizar ou não pacing
#define RETRO 0				// Flag para calcular ou não retropropagação
#define COOLDOWN 0			// Flag que controla o momento de ativação do estímulo das células

const int plot[2] = {1,42};		// Índice do volume de controle que iremos plotar ao longo do tempo

// Define o tipo Func, ponteiro para função que recebe t, um vetor y e o intervalo de tempo k como parâmetro e retorna a avaliação da função
typedef
	double (*Func) (int k, double t, double *y);

/* ============================ FUNÇAO PRINCIPAL ================================================================================================ */
int main (int argc, char *argv[])
{
	if (argc-1 < 3)
	{
		cout << "================ MONODOMAIN =====================" << endl;
		cout << "Solving with -> Euler Explicit (ODE and PDE)" << endl;
		cout << "Usage:> ./monodomain <t_max> <Dx> <Dt> (execution_number)" << endl;
		cout << "<t_max> = Maximum time" << endl;
		cout << "<Dx> = Size of the subinterval in space (um)" << endl;
		cout << "<Dt> = Size of the subinterval in time (s)" << endl;
		cout << "(execution number) = used for the shell script output" << endl;
		cout << "Try this for example: ./monodomain 2.0 50 1.0e-04" << endl;
		cout << "--> After running, open with Paraview monodomain.vtk in the VTK folder" << endl;
		exit(1);
	}
	else
	{
		Node *ptr1, *ptr2;
		FILE *out;
		clock_t begin, end;
		char cmd[50];
		double elapsed;
		double t_max, Dx, Dt, t, maxV;
		double *y0;
		int N, M;
		int j, tick;
		int exec_number;
		bool retro;

		begin = clock();
		// Recebe os parâmetros de entrada
		t_max = atof(argv[1]);
		Dx = atof(argv[2]);
		Dt = atof(argv[3]);
		// Rodando com o shell script ?
		if (argc-1 == 4)
			exec_number = atoi(argv[4]);
		// Rodando simplesmente ?
		else
			exec_number = 0;
		// Criando os diretorios para armazenar os dados da simulacao
		if (exec_number == 1 || exec_number == 0)
			tick = system("mkdir ../Runs");
		sprintf(cmd,"mkdir ../Runs/Run%d",exec_number);
		tick = system(cmd);
		sprintf(cmd,"mkdir ../Runs/Run%d/VTK",exec_number);
		tick = system(cmd);
		sprintf(cmd,"mkdir ../Runs/Run%d/Graphics",exec_number);
		tick = system(cmd);

		// Captura as condições inicias do modelo
		y0 = getInitialConditions();
		maxV = y0[0];

		// Funções da EDO
		Func *f = getFunctions();

		// Cria o grafo que representa as fibras de Purkinje microscópica
		Graph *g = new Graph(y0,num_equacoes,Dx,exec_number);
		// Seta os ponteiros para o cálculo da velocidade de propagação
		ptr1 = g->searchNode(plot[0]);
		ptr2 = g->searchNode(plot[1]);

		printCellsVTK(g,exec_number);
		printInfoModel(g,y0,Dx,Dt,t_max,exec_number,fib_in_bundle,prob_gapJ,D,sigma_c,sigma_p,sigma_t,PACING,RETRO,COOLDOWN);

		// Setar os parâmetros para resolver as EDOs e as EDPs
		setParameters(Dx,Dt,Cm,D/2,PACING,RETRO,COOLDOWN);

		// Calcula o número de subintervalos no tempo
		M = nearbyint(t_max/Dt);

		// Calcula o número de subintervalos no espaço
		N = nearbyint(g->getXmax()/Dx);

		cout << "==== Solving monodomain equation ... =====" << endl;
		// Iterar para resolver a equação do monodomínio
		for (j = 0; j < M; j++)
		{
			t = j*Dt;
			writeGraphic(g,t,plot,2,num_equacoes,exec_number);

			// Escreve a iteração atual em VTK
			if (M > 1000)
			{
				// Printa só algumas iterações
				if (j % 5 == 0)
					writeVTK(g,j,exec_number);
			}
			else
				writeVTK(g,j,exec_number);

			// Para cada volume de controle resolver a EDO associada
			solveEDO(g,t,j,f,num_equacoes);
			// Agora cada volume possui um V* e calculamos o V_n+1 de cada volume pela EDP
			solveEDP(g,f,t,j,num_equacoes);
			//solveEDP2(g,t,j,N,M,f,num_equacoes,v_gate);
			// Velocidade de propagação
			g->calculateVelocity(t,ptr1,ptr2,Dt,Dx,&maxV,exec_number);
			// Avança o tempo n -> n+1
			nextTimestep(g,num_equacoes);

		}
		// Libera memória
		delete [] y0;
		delete [] f;

		// Plotar gráficos
		makePlot(plot,2,exec_number);

		// Printa informação sobre cooldown
		//printCooldown(g,exec_number);
		// Printa retropropagação
		//retro = printRetropropagation(g,exec_number);
		//printCellsRetroprogVTK(g,exec_number);
		// Imprime cilindro
		writeCylinderVTK(g,exec_number);

		// Se estamos executando via shell script, salvar a execução dizendo se houve ou não retropropagação
		//if (argc-1 == 4)
		//	printExecution(exec_number,retro);

		end = clock();
		elapsed = (double)(end-begin) / (double)CLOCKS_PER_SEC;
		sprintf(cmd,"../Runs/Run%d/info_simulation",exec_number);
		out = fopen(cmd,"a");
		fprintf(out,"\n[+] Simulation completed with sucess!\n");
		fprintf(out,"Time elapsed: %e s\n",elapsed);
		fclose(out);
		cout << "[+] Simulation completed with sucess!" << endl;
		cout << "Time elapsed: " << elapsed << " s" << endl << endl;
	}
}
