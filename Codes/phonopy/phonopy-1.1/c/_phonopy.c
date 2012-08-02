#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <numpy/arrayobject.h>

#define KB 8.6173382568083159E-05

/* Build dynamical matrix */
static PyObject * get_dynamical_matrix(PyObject *self, PyObject *args);

static int get_dynamical_matrix_at_q(double *dynamical_matrix_real,
				     double *dynamical_matrix_image,
				     const int d_prim, const int d_super,
				     const double *force_constants, 
				     const double *q,
				     const double *r,
				     const int *multi,
				     const double *mass,
				     const int *s2p_map, 
				     const int *p2s_map);

/* Thermal properties */
static double get_free_energy_omega(double temperature, double omega);
static double get_entropy_omega(double temperature, double omega);
static double get_heat_capacity_omega(double temperature, double omega);
/* static double get_energy_omega(double temperature, double omega); */
static PyObject * get_thermal_properties(PyObject *self, PyObject *args);


static PyMethodDef functions[] = {
  {"dynamical_matrix", get_dynamical_matrix, METH_VARARGS, "Dynamical matrix"},
  {"thermal_properties", get_thermal_properties, METH_VARARGS, "Thermal properties"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_phonopy(void)
{
  Py_InitModule3("_phonopy", functions, "C-extension for phonopy\n\n...\n");
  return;
}

static PyObject * get_dynamical_matrix(PyObject *self, PyObject *args)
{
  PyArrayObject* dynamical_matrix_real;
  PyArrayObject* dynamical_matrix_image;
  PyArrayObject* force_constants;
  PyArrayObject* r_vector;
  PyArrayObject* q_vector;
  PyArrayObject* multiplicity;
  PyArrayObject* mass;
  PyArrayObject* s2p_map;
  PyArrayObject* p2s_map;

  if (!PyArg_ParseTuple(args, "OOOOOOOOO", &dynamical_matrix_real,
			&dynamical_matrix_image,
			&force_constants, &q_vector,
			&r_vector,
			&multiplicity,
			&mass,
			&s2p_map,
			&p2s_map ) )
    return NULL;

  int i;
  double* dm_r = (double*)dynamical_matrix_real->data;
  double* dm_i = (double*)dynamical_matrix_image->data;
  const double* fc = (double*)force_constants->data;
  const double* q = (double*)q_vector->data;
  const double* r = (double*)r_vector->data;
  const double* m = (double*)mass->data;
  const long* multi_long = (long*)multiplicity->data;
  const long* s2p_map_long = (long*)s2p_map->data;
  const long* p2s_map_long = (long*)p2s_map->data;
  const int d_prim = p2s_map->dimensions[0];
  const int d_super = s2p_map->dimensions[0];

  int *multi, *s2p_map_int, *p2s_map_int;

  multi = (int*) malloc( d_prim * d_super * sizeof( int ) );
  for (i = 0; i < d_prim*d_super; i++) {
    multi[i] = (int)multi_long[i];
  }

  s2p_map_int = (int*) malloc( d_super * sizeof( int ) );
  for (i = 0; i < d_super; i++) {
    s2p_map_int[i] = (int)s2p_map_long[i];
  }

  p2s_map_int = (int*) malloc( d_prim * sizeof( int ) );
  for (i = 0; i < d_prim; i++) {
    p2s_map_int[i] = (int)p2s_map_long[i];
  }

  get_dynamical_matrix_at_q( dm_r,
			     dm_i,
			     d_prim,
			     d_super,
			     fc,
			     q,
			     r,
			     multi,
			     m,
			     s2p_map_int,
			     p2s_map_int );

  free( multi );
  free( s2p_map_int );
  free( p2s_map_int );

  Py_RETURN_NONE;
}

static int get_dynamical_matrix_at_q(double *dynamical_matrix_real,
				     double *dynamical_matrix_image,
				     const int d_prim, 
				     const int d_super,
				     const double *force_constants,
				     const double *q,
				     const double *r,
				     const int *multi,
				     const double *mass,
				     const int *s2p_map, 
				     const int *p2s_map )
{
  int i, j, k, l, m;
  double phase, cos_phase, sin_phase, mass_sqrt;
  double dm_real[3][3], dm_imag[3][3];
 

#pragma omp parallel for private( j, k, l, m, phase, cos_phase, sin_phase, mass_sqrt, dm_real, dm_imag ) 
  for ( i = 0; i < d_prim; i++ ) { /* left index of dm */

    for ( j = 0; j < d_prim; j++ ) { /* right index of dm */
      mass_sqrt = sqrt( mass[i] * mass[j] );

      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++ ) {
	  dm_real[k][l] = 0;
	  dm_imag[k][l] = 0;
	}
      }

      for ( k = 0; k < d_super; k++ ) { /* Lattice points of right index of fc */
	if ( ! ( s2p_map[k] == p2s_map[j] ) ) {
	  continue;
	}

	cos_phase = 0;
	sin_phase = 0;
	for ( l = 0; l < multi[ k * d_prim + i ]; l++ ) {
	  phase = 0;
	  for ( m = 0; m < 3; m++ ) {
	    phase += q[m] * r[  k * d_prim*81 + i*81 + l*3 + m ];
	  }
	  cos_phase += cos(phase * 2 * M_PI) / multi[ k * d_prim + i ];
	  sin_phase += sin(phase * 2 * M_PI) / multi[ k * d_prim + i ];
	}

	for ( l = 0; l < 3; l++ ) {
	  for ( m = 0; m < 3; m++ ) {
	    dm_real[l][m] += force_constants[ p2s_map[i] * d_super*9 + k*9 + l*3 + m ] * cos_phase / mass_sqrt;
	    dm_imag[l][m] += force_constants[ p2s_map[i] * d_super*9 + k*9 + l*3 + m ] * sin_phase / mass_sqrt;
	  }
	}
      }
      
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++ ) {
	  dynamical_matrix_real[ ( i*3 + k ) * d_prim * 3 + j*3 + l ] += dm_real[k][l];
	  dynamical_matrix_image[ ( i*3 + k ) * d_prim * 3 + j*3 + l ] += dm_imag[k][l];
	}
      }
    }
  }

  return 0;
}


/* Thermal properties */
static PyObject * get_thermal_properties(PyObject *self, PyObject *args)
{
  double temperature;
  double factor; /* conversion factor to eV */
  double cutoff; /* Cutoff eigenvalue not to be calculated */
  PyArrayObject* eigenvalues;
  PyArrayObject* weights;

  if (!PyArg_ParseTuple( args, "dOOdd",
			 &temperature,
			 &eigenvalues,
			 &weights,
			 &factor,
			 &cutoff ))
    return NULL;

  const double* eigs = (double*)eigenvalues->data;
  const long* w = (long*)weights->data;
  const int num_qpoints = eigenvalues->dimensions[0];
  const int num_bands = eigenvalues->dimensions[1];

  int i, j;
  long sum_weights = 0;
  double free_energy = 0;
  double entropy = 0;
  double heat_capacity = 0;
  double omega = 0;

  for ( i = 0; i < num_qpoints; i++ ){
    sum_weights += w[i];
    for ( j = 0; j < num_bands; j++ ){
      if ( eigs[i*num_bands+j] > cutoff ) {
	omega = sqrt( eigs[i*num_bands+j] ) * factor;
	free_energy += get_free_energy_omega(temperature, omega) * w[i];
	entropy += get_entropy_omega(temperature, omega) * w[i];
	heat_capacity += get_heat_capacity_omega(temperature, omega)* w[i];
      }
    }
  }

  return PyTuple_Pack( 3,
		       PyFloat_FromDouble(free_energy / sum_weights), 
		       PyFloat_FromDouble(entropy / sum_weights),
		       PyFloat_FromDouble(heat_capacity / sum_weights) );
}

static double get_free_energy_omega(double temperature, double omega){
  /* temperature is defined by T (K) */
  /* omega must be normalized to eV. */
  return KB * temperature * log(1 - exp(- omega / (KB * temperature)));
}

static double get_entropy_omega(double temperature, double omega){
  /* temperature is defined by T (K) */
  /* omega must be normalized to eV. */
  double val;

  val = omega / (2 * KB * temperature);
  return 1 / (2 * temperature) * omega * cosh(val) / sinh(val) - KB * log(2 * sinh(val));
}

static double get_heat_capacity_omega(double temperature, double omega){
  /* temperature is defined by T (K) */
  /* omega must be normalized to eV. */
  /* If val is close to 1. Then expansion is used. */
  double val, val1, val2;

  val = omega / ( KB * temperature );
  val1 = exp( val );
  val2 = ( val ) / ( val1 - 1 );
  return KB * val1 * val2 * val2;
}

/* static double get_energy_omega(double temperature, double omega){ */
/*   /\* temperature is defined by T (K) *\/ */
/*   /\* omega must be normalized to eV. *\/ */
/*   return omega / ( exp( omega / ( KB * temperature ) ) - 1 ); */
/* } */
