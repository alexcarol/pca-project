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
#include <xmmintrin.h>

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

	for (x = 0; x < grid_size; x++) {

		printf(".");

		x_centre = gcentre(x, grid_span, grid_size);
		__m128 mx_centre = _mm_set1_ps(x_centre);

		for (y = 0; y < grid_size; y++) {

			y_centre = gcentre(y, grid_span, grid_size);
			__m128 my_centre = _mm_set1_ps(y_centre);

			for (z = 0; z < grid_size; z++) {

				z_centre = gcentre(z, grid_span, grid_size);
				__m128 mz_centre = _mm_set1_ps(z_centre);

				phi = 0;
				__m128 phis = _mm_set1_ps(0.0);
				
				for (residue = 1; residue <= This_Structure.length; residue++) {
					for (atom = 1; atom <= This_Structure.Residue[residue].size - 3; atom += 4) {
						__m128 charges = {
							This_Structure.Residue[residue].Atom[atom].charge,
							This_Structure.Residue[residue].Atom[atom+1].charge,
							This_Structure.Residue[residue].Atom[atom+2].charge,
							This_Structure.Residue[residue].Atom[atom+3].charge
						};
						
						if (charges[0] == 0 && charges[1] == 0 && charges[2] == 0 && charges[3] == 0)
							continue;
						
						__m128 xs = {
							This_Structure.Residue[residue].Atom[atom].coord[1],
							This_Structure.Residue[residue].Atom[atom+1].coord[1],
							This_Structure.Residue[residue].Atom[atom+2].coord[1],
							This_Structure.Residue[residue].Atom[atom+3].coord[1]
						};
						
						__m128 ys = {
							This_Structure.Residue[residue].Atom[atom].coord[2],
							This_Structure.Residue[residue].Atom[atom+1].coord[2],
							This_Structure.Residue[residue].Atom[atom+2].coord[2],
							This_Structure.Residue[residue].Atom[atom+3].coord[2]
						};
						
						__m128 zs = {
							This_Structure.Residue[residue].Atom[atom].coord[3],
							This_Structure.Residue[residue].Atom[atom+1].coord[3],
							This_Structure.Residue[residue].Atom[atom+2].coord[3],
							This_Structure.Residue[residue].Atom[atom+3].coord[3]
						};
						
						__m128 distances;
						
						// Calculo distancias
						__m128 diffxs = _mm_sub_ps(xs, mx_centre);
						__m128 diffys = _mm_sub_ps(ys, my_centre);
						__m128 diffzs = _mm_sub_ps(zs, mz_centre);
						
						diffxs = _mm_mul_ps(diffxs, diffxs);
						diffys = _mm_mul_ps(diffys, diffys);
						diffzs = _mm_mul_ps(diffzs, diffzs);
						
						distances = _mm_add_ps(diffxs, diffys);
						distances = _mm_add_ps(distances, diffzs);
						distances = _mm_sqrt_ps(distances);
						distances = _mm_max_ps(distances, _mm_set1_ps(2.0));
						
						__m128 epsilons = _mm_set1_ps(0.0);
						__m128 tmp;
						__m128 tmp2;
						
						tmp = _mm_cmpge_ps(distances, _mm_set1_ps(8.0));
						epsilons = _mm_and_ps(tmp, _mm_set1_ps(80.0));
						
						if ((int)tmp[0] == 0 || (int)tmp[1] == 0 || (int)tmp[2] == 0 || (int)tmp[3] == 0) {
						
							tmp = _mm_cmple_ps(distances, _mm_set1_ps(6.0));
							tmp = _mm_and_ps(tmp, _mm_set1_ps(4.0));
							epsilons = _mm_or_ps(epsilons, tmp);
						
							tmp = _mm_cmpgt_ps(distances, _mm_set1_ps(6.0));
							tmp2 = _mm_cmpeq_ps(epsilons, _mm_set1_ps(0.0));
							tmp = _mm_and_ps(tmp, tmp2);
							tmp2 = _mm_mul_ps(distances, _mm_set1_ps(38.0));
							tmp2 = _mm_sub_ps(tmp2, _mm_set1_ps(224.0));
							tmp = _mm_and_ps(tmp, tmp2);
						
							epsilons = _mm_or_ps(epsilons, tmp);
						}
						
						tmp = _mm_mul_ps(epsilons, distances);
						tmp = _mm_div_ps(charges, tmp);
						phi += tmp[0];
						phi += tmp[1];
						phi += tmp[2];
						phi += tmp[3];
					}
					
					// Los restantes
					for (; atom <= This_Structure.Residue[residue].size; atom++) {
						
						if (This_Structure.Residue[residue].Atom[atom].charge == 0)
							continue;
						
						distance =
						    pythagoras(This_Structure.Residue[residue].Atom[atom].coord[1], This_Structure.Residue[residue].Atom[atom].coord[2],
							       This_Structure.Residue[residue].Atom[atom].coord[3], x_centre, y_centre, z_centre);
						
						if (distance < 2.0)
							distance = 2.0;
						
						if (distance >= 8.0) {
							epsilon = 80;
						}
						else if (distance <= 6.0) {
							epsilon = 4;
						}
						else {
							epsilon = (38 * distance) - 224;
						}
						
						phi += (This_Structure.Residue[residue].Atom[atom].charge / (epsilon * distance));

					}
				}

				grid[gaddress(x, y, z, grid_size)] = (fftw_real) phi;

			}
		}
	}

	printf("\n");

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
