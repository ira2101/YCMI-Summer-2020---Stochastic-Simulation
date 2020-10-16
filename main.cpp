#include "project.h"

#include<iostream>
#include<cstdlib>
#include<cstring>
#include<random>

#define NA (60200.0)

Simulation sim = Simulation(1.0, 1.0);

void create_environment(double runtime, double tau_step) {
	sim = Simulation(runtime, tau_step);
}

void set_compartments(int num_compartments, double volume, double area, double dx) {
	sim.num_compartments = num_compartments;
	sim.compartment = (Compartment *) malloc(sizeof(Compartment) * num_compartments);

	for (int i = 0; i < num_compartments; i++) {
		sim.compartment[i] = Compartment(volume, area, dx);
	}
}

void change_compartment(int comp_num, double volume, double area, double dx) {
	sim.compartment[comp_num].volume = volume;
	sim.compartment[comp_num].area = area;
	sim.compartment[comp_num].dx = dx;
}

void add_species(const char *name, double molecules, int cost, char type,
								 double diffusion_const) {
	for (int i = 0; i < sim.num_compartments; i++) {
		if (sim.compartment[i].species_head == NULL) {
			sim.compartment[i].species_head = (species_t *) malloc(sizeof(species_t));
			sim.compartment[i].species_head->next = NULL;

			sim.compartment[i].species_head->name = (char *) malloc(strlen(name) + 1);
			strcpy(sim.compartment[i].species_head->name, name);
			sim.compartment[i].species_head->molecules = molecules;
			sim.compartment[i].species_head->cost = cost;
			sim.compartment[i].species_head->type = type;
			sim.compartment[i].species_head->diffusion.diffusion_const = diffusion_const;
			sim.compartment[i].species_head->diffusion.r_rate = 0.0;
			sim.compartment[i].species_head->diffusion.l_rate = 0.0;
			sim.compartment[i].species_head->record.on = false;
			sim.compartment[i].species_head->record.vector = NULL;
		}
		else {
			species_t *tail = sim.compartment[i].species_head;
			while (tail->next != NULL) {
				tail = tail->next;
			}
			tail->next = (species_t *) malloc(sizeof(species_t));
			tail->next->next = NULL;

			tail->next->name = (char *) malloc(strlen(name) + 1);
			strcpy(tail->next->name, name);
			tail->next->molecules = molecules;
			tail->next->cost = cost;
			tail->next->type = type;
			tail->next->diffusion.diffusion_const = diffusion_const;
			tail->next->diffusion.r_rate = 0.0;
			tail->next->diffusion.l_rate = 0.0;
			tail->next->record.on = false;
			tail->next->record.vector = NULL;
		}
	}
}

void set_species_molecules(int compartment, char *name, double molecules) {
	species_t *species_tail = sim.compartment[compartment].species_head;
	while (species_tail != NULL) {
		if (strcmp(species_tail->name, name) == 0) {
			species_tail->molecules = molecules;
			break;
		}

		species_tail = species_tail->next;
	}
}

void set_reaction(double forward_reaction_const, double reverse_reaction_const) {
	for (int i = 0; i < sim.num_compartments; i++) {
		sim.compartment[i].reaction.forward_reaction_const = forward_reaction_const;
		sim.compartment[i].reaction.reverse_reaction_const = reverse_reaction_const;
	}
}

void start_reaction() {
	double runtime = sim.runtime;
	double tau_step = sim.tau_step;
	std::default_random_engine generator(sim.time_seed);
	record(0); //records initial molecules

	for (int w = 1; runtime > 0; w++) {
		if (runtime < tau_step) {
			tau_step = runtime;
		}
		runtime -= tau_step;

		start_diffusion(); //determines the diffusion rates

		for (int i = 0; i < sim.num_compartments; i++) {
			double volume = sim.compartment[i].volume;
			double forward_reaction_rate = sim.compartment[i].reaction.forward_reaction_const *
																		 NA * volume;
			double reverse_reaction_rate = sim.compartment[i].reaction.reverse_reaction_const *
																		 NA * volume;
			species_t *species_tail = sim.compartment[i].species_head;
			while (species_tail != NULL) {
				if (species_tail->type == 'f') {
					for (int q = 0; q < species_tail->cost; q++) {
						forward_reaction_rate *= (species_tail->molecules / (volume * NA));
					}
				}
				else if (species_tail->type == 'r') {
					for (int q = 0; q < species_tail->cost; q++) {
						reverse_reaction_rate *= (species_tail->molecules / (volume * NA)) ;
					}
				}
				species_tail = species_tail->next;
			}

			std::poisson_distribution<int> forward_distribution(forward_reaction_rate * tau_step);
			std::poisson_distribution<int> reverse_distribution(reverse_reaction_rate * tau_step);

			double forward_reactions = forward_distribution(generator);
			double reverse_reactions = reverse_distribution(generator);

			species_tail = sim.compartment[i].species_head;
			while (species_tail != NULL) {
				if (species_tail->type == 'f') {
					species_tail->molecules += species_tail->cost * (reverse_reactions - forward_reactions);
				}
				else if (species_tail->type == 'r') {
					species_tail->molecules += species_tail->cost * (forward_reactions - reverse_reactions);
				}

				species_tail = species_tail->next;
			}
		}

		finish_diffusion(); //calculates the values for diffusion
		record(w); //record after each tau_step
	}
}

void record(int w) {
	for (int i = 0; i < sim.num_compartments; i++) {
		species_t *species_tail = sim.compartment[i].species_head;
		while (species_tail != NULL) {
			if (species_tail->record.on == true) {
				species_tail->record.vector[w] = species_tail->molecules; 
			}

			species_tail = species_tail->next;
		}
	}
}

void start_diffusion() {
	if (sim.num_compartments > 1) {
		for (int i = 0; i < sim.num_compartments; i++) {
			species_t *species_tail = sim.compartment[i].species_head;
			while (species_tail != NULL) {
				if (i == 0) { //leftmost compartment
					species_tail->diffusion.l_rate = 0.0;
					species_tail->diffusion.r_rate = species_tail->diffusion.diffusion_const * species_tail->molecules *
																					 sim.compartment[i].area / (sim.compartment[i].volume * sim.compartment[i].dx);
				}
				else if (i == sim.num_compartments - 1) { //rightmost compartment
					species_tail->diffusion.l_rate = species_tail->diffusion.diffusion_const * species_tail->molecules *
																					 sim.compartment[i].area / (sim.compartment[i].volume * sim.compartment[i].dx);
					species_tail->diffusion.r_rate = 0.0;
				}
				else { //somewhere in the middle
					species_tail->diffusion.l_rate = species_tail->diffusion.diffusion_const * species_tail->molecules *
																					 sim.compartment[i].area / (sim.compartment[i].volume * sim.compartment[i].dx);
					species_tail->diffusion.r_rate = species_tail->diffusion.diffusion_const * species_tail->molecules *
																					 sim.compartment[i].area / (sim.compartment[i].volume * sim.compartment[i].dx);
				}

				species_tail = species_tail->next;
			}
		}
	}
}

void finish_diffusion() {
	if (sim.num_compartments > 1) {
		std::default_random_engine generator(sim.time_seed);
		for (int i = 0; i < sim.num_compartments; i++) {
			species_t *species_tail = sim.compartment[i].species_head;
			while (species_tail != NULL) {
				std::poisson_distribution<int> r_loss_distribution(species_tail->diffusion.r_rate * sim.tau_step);
				std::poisson_distribution<int> l_loss_distribution(species_tail->diffusion.l_rate * sim.tau_step);

				double r_loss = r_loss_distribution(generator);
				double l_loss = l_loss_distribution(generator);

				if (i == 0) { //leftmost compartment
					species_tail->molecules -= r_loss;

					species_t *next_comp_species_tail = sim.compartment[i + 1].species_head;
					while (next_comp_species_tail != NULL) {
						if (strcmp(species_tail->name, next_comp_species_tail->name) == 0) {
							next_comp_species_tail->molecules += r_loss;
							break;
						}

						next_comp_species_tail = next_comp_species_tail->next;
					}
				}
				else if (i == sim.num_compartments - 1) { //rightmost compartment
					species_tail->molecules -= l_loss;

					species_t *prior_comp_species_tail = sim.compartment[i - 1].species_head;
					while (prior_comp_species_tail != NULL) {
						if (strcmp(species_tail->name, prior_comp_species_tail->name) == 0) {
							prior_comp_species_tail->molecules += l_loss;
							break;
						}

						prior_comp_species_tail = prior_comp_species_tail->next;
					}
				}
				else { //middle compartment
					species_tail->molecules -= (r_loss + l_loss);

					species_t *next_comp_species_tail = sim.compartment[i + 1].species_head;
					while (next_comp_species_tail != NULL) {
						if (strcmp(species_tail->name, next_comp_species_tail->name) == 0) {
							next_comp_species_tail->molecules += r_loss;
							break;
						}

						next_comp_species_tail = next_comp_species_tail->next;
					}

					species_t *prior_comp_species_tail = sim.compartment[i - 1].species_head;
					while (prior_comp_species_tail != NULL) {
						if (strcmp(species_tail->name, prior_comp_species_tail->name) == 0) {
							prior_comp_species_tail->molecules += l_loss;
							break;
						}

						prior_comp_species_tail = prior_comp_species_tail->next;
					}
				}
	
				species_tail = species_tail->next;
			}
		}
	}
}

void set_record(int compartment, char *spec_name) {
	species_t *species_tail = sim.compartment[compartment].species_head;
	while (species_tail != NULL) {
		if (strcmp(species_tail->name, spec_name) == 0) {
			species_tail->record.on = true;
			species_tail->record.vector = (double *) malloc(sizeof(double) * get_arr_size());
		}
		
		species_tail = species_tail->next;
	}
}

double *get_record(int compartment, char *spec_name) {
	species_t *species_tail = sim.compartment[compartment].species_head;
	while (species_tail != NULL) {
		if (strcmp(species_tail->name, spec_name) == 0) {
			return species_tail->record.vector;
		}

		species_tail = species_tail->next;
	}

	return NULL;
}

double *get_time_vector() {
	double runtime = sim.runtime;
	double tau_step = sim.tau_step;

	double *time_vector = (double *) malloc(sizeof(double) * get_arr_size());

	for (int i = 1; runtime > 0; i++) {
		if (runtime < tau_step) {
			tau_step = runtime;
		}
		runtime -= tau_step;

		time_vector[i] = time_vector[i - 1] + tau_step;
	}

	return time_vector;
}

int get_arr_size() {
	double runtime = sim.runtime;
	double tau_step = sim.tau_step;
	int size = 1;

	while (runtime > 0.0) {
		runtime -= tau_step;
		size++;
	}

	return size;
}
