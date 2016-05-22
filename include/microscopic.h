#ifndef MICROSCOPIC_H
#define MICROSCOPIC_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "purkinje.h"
#include "plot.h"
#include "monodomain.h"
// Deve-se incluir um modelo celular
#include "noble.h"
//#include "mitchell.h"

#define VTK 0						// Flag para printar os arquivos .vtk

using namespace std;

// Define o tipo Func, ponteiro para função que recebe t, um vetor y e o intervalo de tempo k como parâmetro e retorna a avaliação da função
typedef
	double (*Func) (int k, double t, double *y);

class MicroscopicModel
{
private:
  Graph *g;
  double *y0;
	double dx;
	double dt;
	double t_max;
  Func *f;
  Node **plot;
  int exec_number;

  void createDirectory (int exec_number);
public:
  MicroscopicModel (int argc, char *argv[]);
	void solve ();
};

#endif
