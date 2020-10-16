#include "project.h"

#include<iostream>
#include<cstdlib>
#include<cstring>
#include<random>

#define NA (60200.0)

Simulation::Simulation(double runtime, double tau_step) {
	this->runtime = runtime;
	this->tau_step = tau_step;
	this->num_compartments = 0;
	this->time_seed = time(NULL);
	this->compartment = NULL;
}
