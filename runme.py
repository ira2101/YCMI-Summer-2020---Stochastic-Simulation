import ctypes
import sys
import numpy
#import numpy.ctypeslib
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


def numpy_from_pointer(cpointer, size):
    if sys.version_info.major < 3:
        return numpy.frombuffer(numpy.core.multiarray.int_asbuffer(
            ctypes.addressof(cpointer.contents),
            size * numpy.dtype(float).itemsize))
    else:
        buf_from_mem = ctypes.pythonapi.PyMemoryView_FromMemory
        buf_from_mem.restype = ctypes.py_object
        buf_from_mem.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        cbuffer = buf_from_mem(cpointer, size * numpy.dtype(float).itemsize, 0x200)
        return numpy.ndarray((size,), numpy.float, cbuffer, order='C')

double_p = numpy.ctypeslib.ndpointer(dtype=float, ndim=1, flags='CONTIGUOUS')

a = ctypes.c_char_p("a".encode('utf8'))
b = ctypes.c_char_p("b".encode('utf8'))
ab = ctypes.c_char_p("ab".encode('utf8'))

f = ctypes.c_wchar('f')
r = ctypes.c_wchar('r')

try:
    libsim = ctypes.cdll.LoadLibrary('simulate.so')
except OSError:
    print('Missing simulate.so; aborting.')
    print('Did you forget to run "make" first?')
    sys.exit()

create_environment = libsim.create_environment
create_environment.argtypes = [ctypes.c_double, ctypes.c_double]
create_environment.restype = None

set_compartments = libsim.set_compartments
set_compartments.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
set_compartments.restype = None

set_reaction = libsim.set_reaction
set_reaction.argtypes = [ctypes.c_double, ctypes.c_double]
set_reaction.restype = None

add_species = libsim.add_species
add_species.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int, ctypes.c_wchar, ctypes.c_double]
add_species.restype = None

start_reaction = libsim.start_reaction
start_reaction.artypes = []
start_reaction.restype = None

start_diffusion = libsim.start_diffusion
start_diffusion.argtypes = []
start_diffusion.restype = None

set_record = libsim.set_record
set_record.argtypes = [ctypes.c_int, ctypes.c_char_p]
set_record.restype = None

get_record = libsim.get_record
get_record.argtypes = [ctypes.c_int, ctypes.c_char_p]
get_record.restype = double_p

get_arr_size = libsim.get_arr_size
get_arr_size.argtypes = []
get_arr_size.restype = ctypes.c_int

get_time_vector = libsim.get_time_vector
get_time_vector.argtypes = []
get_time_vector.restype = double_p

set_species_molecules = libsim.set_species_molecules
set_species_molecules.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_double]
set_species_molecules.restype = None

change_compartment = libsim.change_compartment
change_compartment.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
change_compartment.restype = None
#------------------------------------------------------------------------------------

create_environment(6.0, .005); #runtime, tau_step
set_compartments(1, 1.0, 1.0, 1.0); #num_comp, volume, area, dx
set_reaction(1.0, 0.1); #forw_react_const, rev_react_const

add_species(a, 60200, 2, f, 1.0); #name, init_molecules, coeff, type, dif_const
add_species(b, 60200, 1, f, 1.0);
add_species(ab, 0, 1, r, 1.0);

compartment_list = range(1)
species_list = [a, b, ab]

for compartment in compartment_list:
	for species in species_list:
		set_record(compartment, species) #compartment, spec_name

start_reaction();


arr_size = get_arr_size()
time_vector = numpy_from_pointer(get_time_vector(), arr_size)

for i in compartment_list:
	ax = plt.axes()

	for q, species, color, label in zip(range(len(species_list)), species_list, ['blue', 'red', 'green'], ['a', 'b', 'ab']):
		ax.plot(time_vector, numpy_from_pointer(get_record(i, species), arr_size), color=color, label=label)

	plt.title('Compartment %d' % (i))
	plt.xlabel('Time (sec)')
	plt.ylabel('Molecules')
	plt.legend()
	plt.ylim([0, 60000])

	if i != (len(compartment_list) - 1):
		plt.figure()

plt.show()
