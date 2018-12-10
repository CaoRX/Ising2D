#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <config/Parameters.hpp>
#include <cmath>
#include <algorithm>
#include <iostream>
const int nn_dx[4] = {0, 0, 1, -1};
const int nn_dy[4] = {1, -1, 0, 0};
const int nnn_dx[4] = {1, 1, -1, -1};
const int nnn_dy[4] = {1, -1, 1, -1};
class ISING_LATTICE
{
private:
	bool **a;
	bool **vis;
	int *f;
	int *vf;
	int n, m;
	double Beta;
	double Energy;
	double Mag;
	double expJ;

public:
	double *corr;
	int corr_size;
	ISING_LATTICE(int N, int M, double T): n(N), m(M), Beta(1.0 / T)
	{
		vis = new bool*[N];
		a = new bool*[N];
		for (int i = 0; i < N; ++i)
		{
			a[i] = new bool[M];
			vis[i] = new bool[M];
			for (int j = 0; j < M; ++j)
			{
				a[i][j] = (random_double() < 0.5) ? true : false;
				vis[i][j] = false;
			}

		}
		f = new int[N * M];
		vf = new int[N * N];
		Mag = Magnets();
		Energy = Energy_Field();
		Energy += Energy_NN();
		if (Config::USE_NNN)
			Energy += Energy_NNN();
		corr_size = std::min(n, m) / 2;
		corr = new double[corr_size];
		expJ = exp(-Config::J_NN);
	}
	int getS(int x, int y)
	{
		return (a[(x + n) % n][(y + m) % m] ? 1 : -1);
	}
	bool getvis(int x, int y)
	{
		return (vis[(x + n) % n][(y + m) % m]);
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
		Mag -= (2.0 * getS(x, y));
		a[(x + n) % n][(y + m) % m] = (!a[(x + n) % n][(y + m) % m]);
	}
	//void setS(int x, int y, bool value)
	//{
	//	a[(x + n) % n][(y + m) % m] = (value == 1);
	//}
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
		return res * 0.5;
	}
	double Energy_NNN()
	{
		double res = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				for (int dir = 0; dir < 4; ++dir)
					res -= Config::J_NNN * (a[i][j] ? 1 : -1) * getS(i + nnn_dx[dir], j + nnn_dy[dir]);
		return res * 0.5;
	}
	double Magnets()
	{
		double res = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				res += getS(i, j);
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
	void SW_Update();
	void Wolff_Update();
	double get_Energy()
	{
		return Energy / (n * m);
	}
	double get_Magnet()
	{
		return Mag / (n * m);
	}
	void update_corr()
	{
		for (int i = 0; i < corr_size; ++i)
		{
			corr[i] = 0.0;
			for (int x = 0; x < n; ++x)
				for (int y = 0; y < m; ++y)
					corr[i] += getS(x, y) * getS(x + i, y + i);
			corr[i] /= (n * m);
		}
	}
	double get_corr(int x)
	{
		return corr[x];
	}
	int gf(int x)
	{
		if (f[x] == x) {
			return x;
		}
		f[x] = gf(f[x]);
		return f[x];
	}
	void combine(int x, int y)
	{
		f[gf(x)] = gf(y);
	}
	int num_1d(int x, int y)
	{
		x = (x + n) % n;
		y = (y + m) % m;
		return x * m + y;
	}
	~ISING_LATTICE()
	{
		for (int i = 0; i < n; ++i) {
			delete []a[i];
			delete []vis[i];
		}
		delete []a;
		delete []vis;
		delete []f;
		delete []vf;
	}
	// test push
};

#endif

