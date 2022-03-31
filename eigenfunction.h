#ifndef EIGENFUNCTION_H
#define EIGENFUNCTION_H

#include <vector>
#include <fstream>
using namespace std;

class Eigenfunction {
        vector<double> xvalues_n{0};
        vector<double> xvalues_p{0};
        vector<double> psi_b{0};                // wavefunction
        vector<double> d_psi_b{0};	// d (psi_b) / dt
        vector<double> psi_f{0};
        vector<double> d_psi_f{0};
        vector<double> eigenvalues{0};
	vector<vector<double>> eigenfunctions{};
	vector<double> xvalues{};
        double energy_l, energy_h, step_size, h_cut, mass, omega;
        int num_steps;

  public:
	void setValues();
        double slope(double x, double wvfn, double energy);
        void solve();
	void plot(ofstream& file);
};

#endif
