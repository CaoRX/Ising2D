#ifndef LATTICE_CPP
#define LATTICE_CPP

#include <config/Parameters.hpp>
#include <lattice/lattice.hpp>

void ISING_LATTICE::Metro_Sweep()
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			if (random_double() < Flip_prob_Metro(i, j))
				flipS(i, j);
}
void ISING_LATTICE::Heatbath_Sweep()
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			if (random_double() < Flip_prob_Heatbath(i, j))
				flipS(i, j);
}

#endif