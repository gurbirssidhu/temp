#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "eigenfunction.h"

using namespace std;
const double epsilon = 1e-10;

void Eigenfunction::setValues() {
	cout << "Enter lower limit of x(negative value):	";
	cin >> xvalues_n[0];
	cout << "Upper limit is set to " << -xvalues_n[0] << '\n';
	cout << "Enter step_size:	";
	cin >> step_size;
	num_steps = -(xvalues_n[0] / step_size);
	xvalues_p[0] = -xvalues_n[0];
	psi_b[0] = epsilon;	// wavefunction at -infinity is 0
	d_psi_b[0] = epsilon;	// its derivative is also zero
	psi_f[0] =  epsilon;
	d_psi_f[0] = epsilon;
	h_cut = 1;
	cout << "Enter mass of particle:	";
	cin >> mass;
	cout << "Enter value of omega:		";
	cin >> omega;
	cout << "Enter lower limit of energy:	";
	cin >> energy_l;
	cout << "Enter upper limit of energy:	";
	cin >> energy_h;
}

double Eigenfunction::slope(double x, double wvfn, double energy){
	return (-2*mass*(energy - 0.5*mass*(omega*omega)*x*x)*wvfn/(h_cut*h_cut));
}

void Eigenfunction::solve() {
  	double energy = energy_l;
	bool flag{true};
	for (int j{0}; j < (energy_h - energy_l)/0.05; j++) {
	  	double norm{0};
		for (int i{0}; i < num_steps; i++) {
			double k1, k2, k3, k4, q1, q2, q3, q4;

			q1 = d_psi_b[i];
			k1 = slope(xvalues_n[i], psi_b[i], energy);

			q2 = d_psi_b[i] + k1*step_size/2;
			k2 = slope(xvalues_n[i] + step_size*0.5, psi_b[i] + step_size*q1/2, energy);

			q3 = d_psi_b[i] + k2*step_size/2;
			k3 = slope(xvalues_n[i] + step_size*0.5, psi_b[i] + step_size*q2/2, energy);

			q4 = d_psi_b[i] + k3*step_size;
			k4 = slope(xvalues_n[i] + step_size, psi_b[i] + step_size*q3, energy);
		
			xvalues_n.push_back(xvalues_n[i] + step_size);
			d_psi_b.push_back(d_psi_b[i] + step_size * (k1 + 2*k2 + 2*k3 +k4)/6);
			psi_b.push_back(psi_b[i] + step_size * (q1 +2*q2 + 2*q3 + q4)/6);
			
			norm += psi_b[i]*psi_b[i];
		}
	
		for (int i{0}; i < num_steps; i++) {
			double k1, k2, k3, k4, q1, q2, q3, q4;

			q1 = d_psi_f[0];
			k1 = slope(xvalues_p[0], psi_f[0], energy);

			q2 = d_psi_f[0] - k1*step_size/2;
			k2 = slope(xvalues_p[0] - step_size*0.5, psi_f[0] - step_size*q1/2, energy);

			q3 = d_psi_f[0] - k2*step_size/2;
			k3 = slope(xvalues_p[0] - step_size*0.5, psi_f[0] - step_size*q2/2, energy);

			q4 = d_psi_f[0] - k3*step_size;
			k4 = slope(xvalues_p[0] - step_size, psi_f[0] - step_size*q3, energy);
		
			xvalues_p.insert(xvalues_p.begin(), xvalues_p[0] - step_size);
			d_psi_f.insert(d_psi_f.begin(), d_psi_f[0] - step_size * (k1 + 2*k2 + 2*k3 +k4)/6);
			psi_f.insert(psi_f.begin(), psi_f[0] - step_size * (q1 +2*q2 + 2*q3 + q4)/6);
			norm += psi_f[0]*psi_f[0];
		}
		
		double diff {(d_psi_b[num_steps]/psi_b[num_steps]) - (d_psi_f[0]/psi_f[0])};
		double add {(psi_b[num_steps]/d_psi_b[num_steps]) - (psi_f[0]/d_psi_f[0])};

		if (abs(diff) < epsilon || abs(add) < epsilon){
			eigenvalues.push_back(energy);
			vector<double> psi = psi_b;
			if (flag == false){
				for (int i{0}; i < psi_f.size(); i++){
					psi_f[i] = -psi_f[i];
				}
			}

			flag = !flag;

			psi.insert(psi.end(), psi_f.begin(), psi_f.end());

			for (int i{0}; i < 2*num_steps; i++)
				psi[i] = psi[i] / sqrt(norm);
			eigenfunctions.push_back(psi);
			cout << flag << '\n';
			cout << energy << '\n';
		}
		energy += 0.05;
		if (j == 0){
			xvalues = xvalues_n;
			xvalues.insert(xvalues.end(), xvalues_p.begin(), xvalues_p.end());
		}
	
		psi_b.erase(psi_b.begin() + 1, psi_b.end());
		d_psi_b.erase(d_psi_b.begin() + 1, d_psi_b.end());
		xvalues_n.erase(xvalues_n.begin() + 1, xvalues_n.end());
		psi_f.erase(psi_f.begin(), psi_f.end() - 1);
		d_psi_f.erase(d_psi_f.begin(), d_psi_f.end() - 1);
		xvalues_p.erase(xvalues_p.begin(), xvalues_p.end() - 1);

		psi_f[0] = abs(psi_f[0]);
	}
}


void Eigenfunction::plot(ofstream& file) {
	for (int i{0}; i < 2*num_steps; i++){
		file << xvalues[i] << "\t\t";
		for (auto & elem:eigenfunctions)
			file << elem[i] << "\t\t";
		file << '\n';
	}
}

