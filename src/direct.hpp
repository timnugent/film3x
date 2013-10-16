#ifndef DIRECT_HPP
#define DIRECT_HPP

#include <pdb.hpp>	

class Pattern{
public:
	Pattern() : funevals(0), verbose(1), rho(0.9), epsilon(1E-6), itermax(5000) {};
	~Pattern(){};
	void set_target(PDB* t){protein = t;};
	void set_rho(double t){rho = t;};
	void set_epsilon(double t){epsilon = t;};
	void set_itermax(int t){itermax = t;};
	void set_verbose(bool t){verbose = t;};
	int hooke(int, double*, double*);

private:
	double best_nearby(double*, double*, double, int);
	double func(double*);
	int funevals;
	bool verbose;
	double rho, epsilon;
	int itermax;
	PDB* protein;
};

#endif
