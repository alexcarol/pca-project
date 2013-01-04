/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

#ifdef USE_AVX
#include <immintrin.h>
#else
#include <xmmintrin.h>
#endif

#include <malloc.h>

void assign_charges(struct Structure This_Structure)
{

/************/

	/* Counters */

	int residue, atom;

/************/

	for (residue = 1; residue <= This_Structure.length; residue++) {
		for (atom = 1; atom <= This_Structure.Residue[residue].size; atom++) {

			This_Structure.Residue[residue].Atom[atom].charge = 0.0;

			/* peptide backbone */

			if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name, " N  ") == 0) {
				if (strcmp(This_Structure.Residue[residue].res_name, "PRO") == 0) {
					This_Structure.Residue[residue].Atom[atom].charge = -0.10;
				} else {
					This_Structure.Residue[residue].Atom[atom].charge = 0.55;
					if (residue == 1)
						This_Structure.Residue[residue].Atom[atom].charge = 1.00;
				}
			}

			if (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name, " O  ") == 0) {
				This_Structure.Residue[residue].Atom[atom].charge = -0.55;
				if (residue == This_Structure.length)
					This_Structure.Residue[residue].Atom[atom].charge = -1.00;
			}

			/* charged residues */

			if ((strcmp(This_Structure.Residue[residue].res_name, "ARG") == 0) && (strncmp(This_Structure.Residue[residue].Atom[atom].atom_name, " NH", 3) == 0))
				This_Structure.Residue[residue].Atom[atom].charge = 0.50;
			if ((strcmp(This_Structure.Residue[residue].res_name, "ASP") == 0) && (strncmp(This_Structure.Residue[residue].Atom[atom].atom_name, " OD", 3) == 0))
				This_Structure.Residue[residue].Atom[atom].charge = -0.50;
			if ((strcmp(This_Structure.Residue[residue].res_name, "GLU") == 0) && (strncmp(This_Structure.Residue[residue].Atom[atom].atom_name, " OE", 3) == 0))
				This_Structure.Residue[residue].Atom[atom].charge = -0.50;
			if ((strcmp(This_Structure.Residue[residue].res_name, "LYS") == 0) && (strcmp(This_Structure.Residue[residue].Atom[atom].atom_name, " NZ ") == 0))
				This_Structure.Residue[residue].Atom[atom].charge = 1.00;

		}
	}

/************/

}

/************************/

#ifdef USE_AVX
#define IN_NFLOATS 8
#define __mtype __m256
#define _set1_ps(a) _mm256_set1_ps(a)
#define _load_ps(a) _mm256_load_ps(a)
#define _sub_ps(a,b) _mm256_sub_ps(a,b)
#define _add_ps(a,b) _mm256_add_ps(a,b)
#define _mul_ps(a,b) _mm256_mul_ps(a,b)
#define _div_ps(a,b) _mm256_div_ps(a,b)
#define _sqrt_ps(a) _mm256_sqrt_ps(a)
#define _max_ps(a,b) _mm256_max_ps(a,b)
#define _and_ps(a,b) _mm256_and_ps(a,b)
#define _or_ps(a,b) _mm256_or_ps(a,b)
#define _cmpge_ps(a,b) _mm256_cmp_ps(a,b,_CMP_GE_OS)
#define _cmple_ps(a,b) _mm256_cmp_ps(a,b,_CMP_LE_OS)
#define _cmpgt_ps(a,b) _mm256_cmp_ps(a,b,_CMP_GT_OS)
#define _cmpeq_ps(a,b) _mm256_cmp_ps(a,b,_CMP_EQ_OQ)
#else
#define IN_NFLOATS 4
#define __mtype __m128
#define _set1_ps(a) _mm_set1_ps(a)
#define _load_ps(a) _mm_load_ps(a)
#define _sub_ps(a,b) _mm_sub_ps(a,b)
#define _add_ps(a,b) _mm_add_ps(a,b)
#define _mul_ps(a,b) _mm_mul_ps(a,b)
#define _div_ps(a,b) _mm_div_ps(a,b)
#define _sqrt_ps(a) _mm_sqrt_ps(a)
#define _max_ps(a,b) _mm_max_ps(a,b)
#define _and_ps(a,b) _mm_and_ps(a,b)
#define _or_ps(a,b) _mm_or_ps(a,b)
#define _cmpge_ps(a,b) _mm_cmpge_ps(a,b)
#define _cmple_ps(a,b) _mm_cmple_ps(a,b)
#define _cmpgt_ps(a,b) _mm_cmpgt_ps(a,b)
#define _cmpeq_ps(a,b) _mm_cmpeq_ps(a,b)
#endif

// Estructura auxiliar
struct atom_values {
	float xs[IN_NFLOATS];
	float ys[IN_NFLOATS];
	float zs[IN_NFLOATS];
	float charges[IN_NFLOATS];
};

void electric_field(struct Structure This_Structure, float grid_span, int grid_size, fftw_real * grid)
{

/************/

	/* Counters */

	int residue, atom;

	/* Co-ordinates */

	int x, y, z;
	float x_centre, y_centre, z_centre;

	/* Variables */

	float distance;
	float phi, epsilon;

/************/

	for (x = 0; x < grid_size; x++) {
		for (y = 0; y < grid_size; y++) {
			for (z = 0; z < grid_size; z++) {

				grid[gaddress(x, y, z, grid_size)] = (fftw_real) 0;

			}
		}
	}

/************/

	setvbuf(stdout, (char *)NULL, _IONBF, 0);

	printf("  electric field calculations ( one dot / grid sheet ) ");
	
	// Bien, vamos a arrejuntar todo lo interesante...
	int natoms = 0;
	for (residue = 1; residue <= This_Structure.length; residue++)
		natoms += This_Structure.Residue[residue].size;
	
	int natoms_in = natoms % IN_NFLOATS == 0 ? natoms / IN_NFLOATS : natoms / IN_NFLOATS + 1;
	
	// El array debe estar alineado a IN_NFLOATS * 4 bytes
	struct atom_values *atoms = memalign(IN_NFLOATS * 4, natoms_in * sizeof(struct atom_values));
	
	int i = 0;
	for (residue = 1; residue <= This_Structure.length; residue++)
		for (atom = 1; atom <= This_Structure.Residue[residue].size; atom++) {
			if (This_Structure.Residue[residue].Atom[atom].charge == 0) {
				natoms--;
				continue;
			}
			
			atoms[i/IN_NFLOATS].xs[i%IN_NFLOATS] = This_Structure.Residue[residue].Atom[atom].coord[1];
			atoms[i/IN_NFLOATS].ys[i%IN_NFLOATS] = This_Structure.Residue[residue].Atom[atom].coord[2];
			atoms[i/IN_NFLOATS].zs[i%IN_NFLOATS] = This_Structure.Residue[residue].Atom[atom].coord[3];
			atoms[i/IN_NFLOATS].charges[i%IN_NFLOATS] = This_Structure.Residue[residue].Atom[atom].charge;
			i++;
		}
	
	// Me aseguro de que todos los átomos quedan inicializados en caso de que
	// su número no sea múltiplo de IN_NFLOATS. Los átomos con carga 0 no afectan al
	// cálculo ya que la carga se usa para multiplicar el incremento de phi,
	// así que es seguro computar estos "átomos" de más.
	for (; i % IN_NFLOATS != 0; i++) {
		atoms[i/IN_NFLOATS].xs[i%IN_NFLOATS] = 0;
		atoms[i/IN_NFLOATS].ys[i%IN_NFLOATS] = 0;
		atoms[i/IN_NFLOATS].zs[i%IN_NFLOATS] = 0;
		atoms[i/IN_NFLOATS].charges[i%IN_NFLOATS] = 0;
	}
	
	natoms_in = natoms % IN_NFLOATS == 0 ? natoms / IN_NFLOATS : natoms / IN_NFLOATS + 1;

	for (x = 0; x < grid_size; x++) {

		printf(".");

		x_centre = gcentre(x, grid_span, grid_size);
		__mtype mx_centre = _set1_ps(x_centre);

		for (y = 0; y < grid_size; y++) {

			y_centre = gcentre(y, grid_span, grid_size);
			__mtype my_centre = _set1_ps(y_centre);

			for (z = 0; z < grid_size; z++) {

				z_centre = gcentre(z, grid_span, grid_size);
				__mtype mz_centre = _set1_ps(z_centre);

				phi = 0;
				__mtype phis = _set1_ps(0.0);
				
				for (i = 0; i < natoms_in; i++) {
				
					__mtype xs = _load_ps(atoms[i].xs);
					__mtype ys = _load_ps(atoms[i].ys);
					__mtype zs = _load_ps(atoms[i].zs);
					__mtype charges = _load_ps(atoms[i].charges);
					__mtype distances;
					
					// Calculo distancias (el original pythagoras)
					__mtype diffxs = _sub_ps(xs, mx_centre);
					__mtype diffys = _sub_ps(ys, my_centre);
					__mtype diffzs = _sub_ps(zs, mz_centre);
					
					diffxs = _mul_ps(diffxs, diffxs);
					diffys = _mul_ps(diffys, diffys);
					diffzs = _mul_ps(diffzs, diffzs);
					
					distances = _add_ps(diffxs, diffys);
					distances = _add_ps(distances, diffzs);
					
					distances = _sqrt_ps(distances);
					
					// A partir de aquí implemento los if's originales usando solo máscaras de bits
					
					// Trunco a 2 como mínimo
					distances = _max_ps(distances, _set1_ps(2.0));
					
					__mtype epsilons = _set1_ps(0.0);
					__mtype tmp;
					__mtype tmp2;
					
					// if >= 8
					tmp = _cmpge_ps(distances, _set1_ps(8.0));
					epsilons = _and_ps(tmp, _set1_ps(80.0));
					
					// else if <= 6
					tmp = _cmple_ps(distances, _set1_ps(6.0));
					tmp = _and_ps(tmp, _set1_ps(4.0));
					epsilons = _or_ps(epsilons, tmp);
					
					// else
					tmp = _cmpgt_ps(distances, _set1_ps(6.0));
					tmp2 = _cmpeq_ps(epsilons, _set1_ps(0.0));
					tmp = _and_ps(tmp, tmp2);
					tmp2 = _mul_ps(distances, _set1_ps(38.0));
					tmp2 = _sub_ps(tmp2, _set1_ps(224.0));
					tmp = _and_ps(tmp, tmp2);
			
					// Valor final
					epsilons = _or_ps(epsilons, tmp);
					
					// Calculo las phis
					tmp = _mul_ps(epsilons, distances);
					tmp = _div_ps(charges, tmp);
					
					// Acumulo las phis
					phis = _add_ps(phis, tmp);
				}
				#ifdef USE_AVX
				
				phi += phis[0];
				phi += phis[1];
				phi += phis[2];
				phi += phis[3];
				phi += phis[4];
				phi += phis[5];
				phi += phis[6];
				phi += phis[7];
				
				#else
				
				phi += phis[0];
				phi += phis[1];
				phi += phis[2];
				phi += phis[3];
				
				#endif
				grid[gaddress(x, y, z, grid_size)] = (fftw_real) phi;

			}
		}
	}

	printf("\n");

	free(atoms);
/************/

	return;

}

/************************/

void electric_point_charge(struct Structure This_Structure, float grid_span, int grid_size, fftw_real * grid)
{

/************/

	/* Counters */

	int residue, atom;

	/* Co-ordinates */

	int x, y, z;
	int x_low, x_high, y_low, y_high, z_low, z_high;

	float a, b, c;
	float x_corner, y_corner, z_corner;
	float w;

	/* Variables */

	float one_span;

/************/

	for (x = 0; x < grid_size; x++) {
		for (y = 0; y < grid_size; y++) {
			for (z = 0; z < grid_size; z++) {

				grid[gaddress(x, y, z, grid_size)] = (fftw_real) 0;

			}
		}
	}

/************/

	one_span = grid_span / (float)grid_size;

	for (residue = 1; residue <= This_Structure.length; residue++) {
		for (atom = 1; atom <= This_Structure.Residue[residue].size; atom++) {

			if (This_Structure.Residue[residue].Atom[atom].charge != 0) {

				x_low = gord(This_Structure.Residue[residue].Atom[atom].coord[1] - (one_span / 2), grid_span, grid_size);
				y_low = gord(This_Structure.Residue[residue].Atom[atom].coord[2] - (one_span / 2), grid_span, grid_size);
				z_low = gord(This_Structure.Residue[residue].Atom[atom].coord[3] - (one_span / 2), grid_span, grid_size);

				x_high = x_low + 1;
				y_high = y_low + 1;
				z_high = z_low + 1;

				a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre(x_low, grid_span, grid_size) - (one_span / 2);
				b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre(y_low, grid_span, grid_size) - (one_span / 2);
				c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre(z_low, grid_span, grid_size) - (one_span / 2);

				for (x = x_low; x <= x_high; x++) {

					x_corner = one_span * ((float)(x - x_high) + .5);

					for (y = y_low; y <= y_high; y++) {

						y_corner = one_span * ((float)(y - y_high) + .5);

						for (z = z_low; z <= z_high; z++) {

							z_corner = one_span * ((float)(z - z_high) + .5);

							w = ((x_corner + a) * (y_corner + b) * (z_corner + c)) / (8.0 * x_corner * y_corner * z_corner);

							grid[gaddress(x, y, z, grid_size)] += (fftw_real) (w * This_Structure.Residue[residue].Atom[atom].charge);

						}
					}
				}

			}

		}
	}

/************/

	return;

}

/************************/

void electric_field_zero_core(int grid_size, fftw_real * elec_grid, fftw_real * surface_grid, float internal_value)
{

/************/

	/* Co-ordinates */

	int x, y, z;

/************/

	for (x = 0; x < grid_size; x++) {
		for (y = 0; y < grid_size; y++) {
			for (z = 0; z < grid_size; z++) {

				if (surface_grid[gaddress(x, y, z, grid_size)] == (fftw_real) internal_value)
					elec_grid[gaddress(x, y, z, grid_size)] = (fftw_real) 0;

			}
		}
	}

/************/

	return;

}
