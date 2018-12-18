#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

namespace Config {
const double J_NN = 100.0;
const double J_NNN = 0.0;
const double J_NNN_1 = 1.0;
const double J_NNN_2 = J_NN;
const double Magnet_Field = 6.0 * J_NN;
const bool USE_NNN = true;
const int default_n = 100;
const int default_m = 100;
}

// H = sum J_NN * Si * Sj + sum J_NNN * Si * Sj
extern double random_double();
extern int random_int(int n) ;
#endif
