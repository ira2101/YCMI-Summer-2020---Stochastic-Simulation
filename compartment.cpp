#include "project.h"

#include<cstdlib>
#include<iostream>

Compartment::Compartment(double volume, double area, double dx) {
	//this->reaction = {};
	this->species_head = NULL;
	this->volume = volume;
	this->area = area;
	this->dx = dx;
}
