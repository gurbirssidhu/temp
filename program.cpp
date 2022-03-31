#include "eigenfunction.h"
#include <fstream>

int main() {
	Eigenfunction o1;
	o1.setValues();
	o1.solve();
	
	ofstream file("energy.dat");
	o1.plot(file);

	return 0;
}
