#ifndef LATTICE_CPP
#define LATTICE_CPP

#include <config/Parameters.hpp>
#include <lattice/lattice.hpp>
#include <vector>

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

void ISING_LATTICE::SW_Update() 
{
	for (int i = 0; i < n * m; ++i) {
		f[i] = i;
		vf[i] = 0;
	}
	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < m; ++y) {
			if (getS(x, y) == getS(x, y + 1) && random_double() > expJ) {
				combine(num_1d(x, y), num_1d(x, y + 1));
			}
			if (getS(x, y) == getS(x + 1, y) && random_double() > expJ) {
				combine(num_1d(x, y), num_1d(x + 1, y));
			}
		}
	}
	int tmp = 0, ftmp = 0;
	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < m; ++y) {
			tmp = num_1d(x, y);
			ftmp = gf(tmp);
			if (!vf[ftmp]) {
				vf[ftmp] = (random_double() < 0.5) ? 1 : -1;
			}
			if (vf[ftmp] == -1) {
				flipS(x, y);
			}
		}
	}
}
void ISING_LATTICE::Wolff_Update() 
{
	std::vector<int> qx, qy;
	int bx = random_int(n), by = random_int(m);
	qx.push_back(bx);
	qy.push_back(by);
	vis[bx][by] = true;
	int hd = 0, tl = 1, ts = 0;
	int nx, ny;
	while (hd < tl)
	{
		bx = qx[hd];
		by = qy[hd];
		ts = getS(bx, by);
		for (int i = 0; i < 4; ++i) {
			nx = ((bx + nn_dx[i]) + n)	% n;
			ny = ((by + nn_dy[i]) + m) % m;
			if (ts == getS(nx, ny) && (!getvis(nx, ny)) && random_double() > expJ) {
				qx.push_back(nx);
				qy.push_back(ny);
				vis[nx][ny] = true;
				++tl;
			}
		}
		++hd;
	}
	if (random_double() < 0.5) {
		for (int i = 0; i < qx.size(); ++i) {
			flipS(qx[i], qy[i]);
		}
	}
	for (int i = 0; i < qx.size(); ++i)
		vis[qx[i]][qy[i]] = false;
	qx.clear();
	qy.clear();
}
#endif





