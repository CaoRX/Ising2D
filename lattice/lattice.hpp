#ifndef LATTICE_HPP
#define LATTICE_HPP

class ISING_LATTICE
{
private:
	bool **a;
	int n;

public:
	ISING_LATTICE(int N): n(N)
	{
		a = new bool*[N];
		for (int i = 0; i < N; ++i)
		{
			a[i] = new bool[N];
			for (int j = 0; j < N; ++j)
				a[i][j] = 0;
		}
	}
}
#endif

