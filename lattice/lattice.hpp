#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <config/Parameters.hpp>
#include <cmath>
#include <algorithm>
const int nn_dx[4] = {0, 0, 1, -1};
const int nn_dy[4] = {1, -1, 0, 0};
const int nnn_dx[4] = {1, 1, -1, -1};
const int nnn_dy[4] = {1, -1, 1, -1};
class ISING_LATTICE
{
private:
	bool **a;
	int n, m;
	double Energy;
	double Beta;

public:
	ISING_LATTICE(int N, int M): n(N), m(M)
	{
		a = new bool*[N];
		for (int i = 0; i < N; ++i)
		{
			a[i] = new bool[M];
			for (int j = 0; j < M; ++j)
				a[i][j] = 0;
		}
	}
	int getS(int x, int y)
	{
		return (a[(x + n) % n][(y + m) % m] ? 1 : -1);
	}
	double local_Energy(int x, int y)
	{
		x = (x + n) % n;
		y = (y + m) % m;
		int sxy = getS(x, y);
		double res = -Config::Magnet_Field * sxy;
		for (int i = 0; i < 4; ++i)
			res -= Config::J_NN * sxy * getS(x + nn_dx[i], y + nn_dy[i]);
		if (Config::USE_NNN)
		for (int i = 0; i < 4; ++i)
			res -= Config::J_NNN * sxy * getS(x + nnn_dx[i], y + nnn_dy[i]);
		return res;
	}
	double flip_delta_Energy(int x, int y)
	{
		return -2.0 * local_Energy(x, y);
	}
	void flipS(int x, int y)
	{
		Energy += flip_delta_Energy(x, y);
		a[(x + n) % n][(y + m) % m] = (!a[(x + n) % n][(y + m) % m]);
	}
	void setS(int x, int y, bool value)
	{
		a[(x + n) % n][(y + m) % m] = (value == 1);
	}
	double Energy_Field()
	{
		double res = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res -= Config::Magnet_Field * (a[i][j] ? 1 : -1);
		return res;
	}
	double Energy_NN()
	{
		double res = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				for (int dir = 0; dir < 4; ++dir)
					res -= Config::J_NN * (a[i][j] ? 1 : -1) * getS(i + nn_dx[dir], j + nn_dy[dir]);
		return res;
	}
	double Energy_NNN()
	{
		double res = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				for (int dir = 0; dir < 4; ++dir)
					res -= Config::J_NNN * (a[i][j] ? 1 : -1) * getS(i + nnn_dx[dir], j + nnn_dy[dir]);
		return res;
	}

	double Flip_prob_Metro(int x, int y)
	{
		return std::min(1.0, exp(-Beta * flip_delta_Energy(x, y)));
	}
	double Flip_prob_Heatbath(int x, int y)
	{
		return 1.0 / (1.0 + exp(Beta * flip_delta_Energy(x, y)));
	}
	void Metro_Sweep();
	void Heatbath_Sweep();
	//void Metro_Sweep()
	//{
	//	for (int i = 0; i <)
	//}
};

#endif

