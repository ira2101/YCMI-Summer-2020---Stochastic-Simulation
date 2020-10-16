#ifndef PROJECT_H
#define PROJECT_H

#include<ctime>

typedef struct diffusion {
	double diffusion_const;

	double r_rate;
	double l_rate;
} diffusion_t;

typedef struct record {
	bool on;
	double *vector;
} record_t;

typedef struct species {
	double molecules;
	int cost;
	char type;
	char *name;
	diffusion_t diffusion;
	record_t record;
	struct species *next;
} species_t;

typedef struct reaction {
	double forward_reaction_const;
	double reverse_reaction_const;
} reaction_t;

class Compartment {
	public:
		reaction_t reaction;
		species_t *species_head;
		double volume;
		double area;
		double dx;

		Compartment(double volume, double area, double dx);
};

class Simulation {
	public:
		Compartment *compartment;
		int num_compartments;
		double runtime;
		double tau_step;
		time_t time_seed;

		Simulation(double runtime, double tau_step);
};

extern Simulation sim;

extern "C" {
	void create_environment(double runtime, double tau_step);
	void set_compartments(int num_compartments, double volume, double area, double dx);
	void change_compartment(int comp_num, double volume, double area, double dx);
	void set_reaction(double forward_reaction_const, double reverse_reaction_const);
	void add_species(const char *name, double molecules, int cost, char type, double diffusion_const);
	void set_species_molecules(int compartment, char *name, double molecules);

	void start_reaction();
	void start_diffusion(); //sets the rates
	void finish_diffusion(); //changes the values
	void set_record(int compartment, char *spec_name);
	double *get_record(int compartment, char *spec_name);
	void record(int w);
	double *get_time_vector();
	int get_arr_size();
}
#endif
