#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include "purkinje.h"

using namespace std;

#define CYCLE_L 4000				// Intervalo de tempo do pacing 

// Define o tipo Func, ponteiro para função que recebe t, e um vetor y como parâmetro e retorna a avaliação da função
typedef
	double (*Func) (int k, double t, double *y);

void solveEDO (Graph *g, double t, int k, Func *func, int num_eq);
void solveEDP (Graph *g, Func *func, double t, int k, int num_eq);
void solveEDP2 (Graph *g, double t, int k, int N, int M, Func *func, int num_equacoes, double v_gate);
void nextTimestep (Graph *g, int num_eq);
void setParameters (double Dx, double Dt, double cm, double R, int PACING, int RETROP, int COOLDOWN);
void checkRetropropagation (Node *ptr);
bool printRetropropagation (Graph *g, int exec_number);
void printCooldown (Graph *g);
void checkCFL (double r);

#endif
