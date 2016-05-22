#include <iostream>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include "microscopic.h"

using namespace std;

/* ============================ FUNÇAO PRINCIPAL ================================================================================================ */
int main (int argc, char *argv[])
{
	if (argc-1 < 3)
	{
		cout << "=============== MICROSCOPIC MODEL =====================" << endl;
		cout << "Usage:> ./microscopic <t_max> <Dx> <Dt> (execution_number)" << endl;
		cout << "<t_max> = Maximum time of the simulation (s)" << endl;
		cout << "<Dx> = Timestep in space (um)" << endl;
		cout << "<Dt> = Timestep in time (s)" << endl;
		cout << "(execution number) = used for the shell script output" << endl;
		cout << "Try this for example: ./microscopic 1.2 25 1.0e-04" << endl;
		cout << "--> After running, open with Paraview monodomain.vtk in the VTK folder" << endl;
		exit(1);
	}
	else
	{
		MicroscopicModel *model = new MicroscopicModel(argc,argv);
		model->solve();
	}
	return 0;
}
