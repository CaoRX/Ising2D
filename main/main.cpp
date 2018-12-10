#ifndef MAIN_CPP
#define MAIN_CPP
#include <lattice/lattice.hpp>
#include <config/Parameters.hpp>
#include <iostream>
#include <cstdio>
#include <vector>


std::vector<double> Energy, Magnet;
double *corr;
int lattice_size = 100;

std::vector<double> Temp, corrl, energy, magnet;
std::vector<int> Latt;
extern double corr_length(double *corr, int corr_size, double stop_point = exp(-8.0));

void calc_ising(double temp, int latt, int loop_num = 10000, int equil = 100, bool heatbath = false)
{
	std::cout << "T = " << temp << ", lattice size = " << latt << std::endl;
	Temp.push_back(temp);
	Latt.push_back(latt);
	ISING_LATTICE IL(latt, latt, temp);
	double *corr = new double[IL.corr_size];
	double Energy_sum = 0.0, Magnet_sum = 0.0;
	int corr_size = IL.corr_size;
	for (int i = 0; i < corr_size; ++i)
		corr[i] = 0.0;
	for (int i = 0; i < loop_num; ++i)
	{
		if (heatbath)
		{
			IL.Heatbath_Sweep();
		}	
		else
		{
			IL.Metro_Sweep();
		}
		
		//IL.Wolff_Update();
		if (i < equil)
		{
			continue;
		}
		IL.update_corr();
		Energy_sum += IL.get_Energy();
		Magnet_sum += IL.get_Magnet();
		for (int i = 0; i < corr_size; ++i)
			corr[i] += IL.corr[i];
	}
	energy.push_back(Energy_sum / (loop_num - equil));
	magnet.push_back(Magnet_sum / (loop_num - equil));
	for (int i = 0; i < corr_size; ++i)
		corr[i] /= (loop_num - equil);
	double corr_l = corr_length(corr, corr_size);
	corrl.push_back(corr_l);;
	std::cout << "correlation length = " << corr_l << std::endl;
	std::cout << "average energy = " << (Energy_sum / (loop_num - equil)) << std::endl;
	std::cout << "average magnet = " << (Magnet_sum / (loop_num - equil)) << std::endl;
}

void update(ISING_LATTICE &IL, bool heatbath)
{
	if (heatbath)
		IL.Heatbath_Sweep();
	else
		IL.Metro_Sweep();
	IL.update_corr();
	Energy.push_back(IL.get_Energy());
	Magnet.push_back(IL.get_Magnet());
	for (int i = 0; i < IL.corr_size; ++i)
		corr[i] += IL.corr[i];
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
void output_corr(FILE *f, int corr_size)
{
	for (int i = 0; i < corr_size; ++i)
		fprintf(f, "%f\n", corr[i]);
}
const int lattice_num = 5;
const int lattice_grid[lattice_num] = {10, 30, 50, 70, 100};
const int temp_num = 10;
const double temp_grid[temp_num] = {1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.0, 5.0};
const double detail_grid[temp_num] = {2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29};
// critical temperature \approx 2.27
void output_lists(FILE *f)
{
	for (int i = 0; i < Temp.size(); ++i) {
		fprintf(f, "%f %d %f %f %f\n", Temp[i], Latt[i], corrl[i], energy[i], magnet[i]);
	}
}
void Initialize()
{
	Temp.clear();
	corrl.clear();
	energy.clear();
	magnet.clear();
	Latt.clear();
}
int main(int argc, char **argv)
{
	//int loop_num = atoi(argv[1]);
	//double T = atof(argv[2]);
	//ISING_LATTICE IL(lattice_size, lattice_size, T);
	//corr = new double[IL.corr_size];

	//for (int i = 0; i < loop_num; ++i)
	//	update(IL, false);
	//for (int i = 0; i < IL.corr_size; ++i)
	//{
	//	corr[i] /= loop_num;
	//}
	//std::cout << "finish updating" << std::endl;
	//FILE *energy_f = fopen("data/energy.dat", "w");
	//FILE *magnet_f = fopen("data/magnet.dat", "w");
	//FILE *corr_f = fopen("data/corr.dat", "w");
	//std::cout << "open file" << std::endl;
	//output_energy(energy_f);
	//output_magnet(magnet_f);
	//output_corr(corr_f, IL.corr_size);
	//fclose(energy_f);
	//fclose(magnet_f);
	//fclose(corr_f);
	//std::cout << corr_length(corr, IL.corr_size) << std::endl;
	for (int i = 0; i < lattice_num; ++i)
		for (int j = 0; j < temp_num; ++j)
			calc_ising(detail_grid[j], lattice_grid[i]);
	FILE *result_file = fopen("data/result.dat", "w");
	output_lists(result_file);
	return 0;
}

#endif