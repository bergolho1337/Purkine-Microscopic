#include "microscopic.h"

// Cria os diretorios necessarios
void MicroscopicModel::createDirectory (int exec_number)
{
  int tick;
  char cmd[50];
  if (exec_number < 2)
    tick = system("mkdir ../Runs");
  sprintf(cmd,"mkdir ../Runs/Run%d",exec_number);
  tick = system(cmd);
  sprintf(cmd,"mkdir ../Runs/Run%d/VTK",exec_number);
  tick = system(cmd);
  sprintf(cmd,"mkdir ../Runs/Run%d/Graphics",exec_number);
  tick = system(cmd);
  cout << "[+] Directories builded with sucess!" << endl;
}

// Construtor do modelo
MicroscopicModel::MicroscopicModel (int argc, char *argv[])
{
  this->t_max = atof(argv[1]);
  this->dx = atof(argv[2]);
  this->dt = atof(argv[3]);
  this->y0 = getInitialConditions();
  this->f = getFunctions();
  // Rodando normal ou com shell script ?
  (argc-1 == 4) ? this->exec_number = atoi(argv[4]) : this->exec_number = 0;
  createDirectory(this->exec_number);
  this->g = new Graph(y0,num_equacoes,this->dx,this->exec_number);
  this->plot = this->g->setPlotPointers();
  printInfoModel(this->g,this->y0,this->dx,this->dt,this->t_max,this->exec_number,fib_in_bundle,prob_gapJ,D,cell_length,diameter_cell,fiber_length,sigma_c,G_p,G_t,Cm);
  cout << "[+] Model builded with sucess!" << endl;
}

// Resolvedor do modelo (Implicit Euler -> EDP e Explicit Euler -> EDO) -> MVF
void MicroscopicModel::solve ()
{
  clock_t begin, end;
  int N, M;
  double t, maxV;
  maxV = this->y0[4];
  // Numero de subintervalos no espaco e no tempo
  N = nearbyint(this->g->getXmax()/this->dx);
  M = nearbyint(this->t_max/this->dt);
  // Setar os parametros para resolver as EDOs e as EDPs
  setParameters(this->dx,this->dt,Cm);
  cout << "[!] Building matrix of the method (Implicit Euler) ..." << endl;
  // Monta a matriz do metodo implicito
  makeMatrix_A(this->g);
  cout << "=========== Solving monodomain equation ... =============" << endl;
  cout << "With N = " << N << " (space)" << endl;
  cout << "With M = " << M << " (time)" << endl;
  begin = clock();
  // Iterar para resolver a equacao do monodominio
  for (int j = 0; j < M; j++)
  {
    t = j*this->dt;
    // Escreve dados do plot do grafico
    writeGraphic(this->g,t,this->plot,num_equacoes,this->exec_number);
    // Escreve a iteração atual em VTK
    writeIteration(this->g,j,M,this->exec_number);
    // Resolve a iteracao atual
    solveEDO(this->g,t,j,this->f,num_equacoes);
    makeVector_b(this->g);
    solveEDP_Imp(this->g);
    // Velocidade de propagação
    this->g->calculateVelocity(t,this->plot[0],this->plot[1],this->dt,this->dx,&maxV,this->exec_number);
    // Avança o tempo n -> n+1
    nextTimestep(this->g,num_equacoes);
  }
  end = clock();
  // Libera memória
  delete [] this->y0;
  delete [] this->f;
  // Plotar gráficos
  makePlot(this->plot,exec_number);
  // Imprime cilindro
  writeCylinderVTK(this->g,exec_number);
  // Imprime tempo
  writeTime(begin,end,this->exec_number);

}
