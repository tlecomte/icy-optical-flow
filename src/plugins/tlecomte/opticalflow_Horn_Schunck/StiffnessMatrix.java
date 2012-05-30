/*
 * Copyright 2010, 2011 Institut Pasteur.
 * 
 * This file is part of ICY.
 * 
 * ICY is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ICY is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with ICY. If not, see <http://www.gnu.org/licenses/>.
 */

package plugins.tlecomte.opticalflow_Horn_Schunck;

import plugins.tlecomte.fem.Element;
import plugins.tlecomte.fem.Node;

public class StiffnessMatrix {
	public static int build_A_HS_COO(int j,
					   				 int[] row,
					   				 int[] col,
					   				 double[] data,	
					   				 double dxI2_e,
					   				 double dyI2_e,
					   				 double alpha,
					   				 int N,
					   				 Element element,
					   				 double[] A_phiphi,
					   				 double[] A_dxphidxphi,
					   				 double[] A_dyphidyphi)
	{
		int i, i1, i2;

		i = 0;
		for (Node node1 : element.nodes) {
			for (Node node2 : element.nodes) {
				i1 = node1.index;
				i2 = node2.index;
				row[j] =     i1; col[j] =     i2; data[j]= Math.pow(dxI2_e, 2)*A_phiphi[i]; j += 1;
				row[j] = N + i1; col[j] = N + i2; data[j]= Math.pow(dyI2_e, 2)*A_phiphi[i]; j += 1;
				row[j] =     i1; col[j] = N + i2; data[j]= dxI2_e*dyI2_e*A_phiphi[i];       j += 1;
				row[j] = N + i1; col[j] =     i2; data[j]= dxI2_e*dyI2_e*A_phiphi[i];       j += 1;           
				row[j] =     i1; col[j] =     i2; data[j]= alpha*A_dxphidxphi[i];           j += 1;
				row[j] = N + i1; col[j] = N + i2; data[j]= alpha*A_dxphidxphi[i];           j += 1;
				row[j] =     i1; col[j] =     i2; data[j]= alpha*A_dyphidyphi[i];           j += 1;
				row[j] = N + i1; col[j] = N + i2; data[j]= alpha*A_dyphidyphi[i];           j += 1; 
				i += 1;
			}
		}

		return j;
	}
}
