#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <config/Parameters.hpp>
class ISING_LATTICE
{
private:
	bool **a;
	int n, m;

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
	bool getS(int x, int y)
	{
		return a[(x + n) % n][(y + m) % m];
	}
	void setS(int x, int y, bool value)
	{
		a[(x + n) % n][(y + m) % m] = value;
	}
};

#endif

