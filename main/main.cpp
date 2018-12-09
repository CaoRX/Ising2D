#ifndef MAIN_CPP
#define MAIN_CPP
#include <lattice/lattice.hpp>
#include <config/Parameters.hpp>
#include <iostream>
#include <cstdio>
#include <vector>


std::vector<double> Energy, Magnet;

void update(ISING_LATTICE &IL, bool heatbath)
{
	if (heatbath)
		IL.Heatbath_Sweep();
	else
		IL.Metro_Sweep();
	Energy.push_back(IL.get_Energy());
	Magnet.push_back(IL.get_Magnet());
	//std::cout << IL.get_Energy() * 100 - IL.Energy_NN() << std::endl;
	//std::cout << IL.get_Energy() * 100 << std::endl;
}
void output_energy(FILE *f)
{
	for (int i = 0; i < Energy.size(); ++i)
		fprintf(f, "%f\n", Energy[i]);
	//std::cout << f << std::endl;
}
void output_magnet(FILE *f)
{
	for (int i = 0; i < Magnet.size(); ++i)
		fprintf(f, "%f\n", Magnet[i]);
}
int main(int argc, char **argv)
{
	int loop_num = atoi(argv[1]);
	double T = atof(argv[2]);
	ISING_LATTICE IL(10, 10, T);
	for (int i = 0; i < loop_num; ++i)
		update(IL, false);
	std::cout << "finish updating" << std::endl;
	FILE *energy_f = fopen("data/energy.dat", "w");
	FILE *magnet_f = fopen("data/magnet.dat", "w");
	//std::cout << "open file" << std::endl;
	output_energy(energy_f);
	output_magnet(magnet_f);
	fclose(energy_f);
	fclose(magnet_f);
	return 0;
}

#endif