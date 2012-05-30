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

package plugins.tlecomte.rectangularFEM;

public class Quadratures {
	
	// compute the integrals numerically with a 1D two point quadrature, which is
	// exact up to order 3
	public static double quadrature_integral_1D_order3(Integrand1D integrand)
	{
		double u_order3 = 1./Math.sqrt(3.);
		double[] U_1D_order3 = {-u_order3, u_order3};
		double sum = 0;
		for (int i=0; i<U_1D_order3.length; i++) {
			sum += integrand.function(U_1D_order3[i]);
		}
		return sum;
	}
	
	// compute the integrals numerically with a 1D two point quadrature, which is
	// exact up to order 3
	static double quadrature_integral_1D_order5(Integrand1D integrand)
	{
		double u_order5 = Math.sqrt(3./5.);
		double wu_order5 = 5./9.;
		double w0_order5 = 8./9.;	
		double[] U_1D_order5 = {-u_order5, 0., u_order5};
		double[] W_1D_order5 = {wu_order5, w0_order5, wu_order5};
		
		double sum = 0;
		for (int i=0; i<U_1D_order5.length; i++) {
			sum += W_1D_order5[i]*integrand.function(U_1D_order5[i]);
		}
		return sum;
	}
	    
	// compute the integrals numerically with a 2D 2x5 point quadrature, which is
	// exact up to order 7 on each direction
	static double quadrature_integral_2D_order7(Integrand2D integrand)
	{
		double u1_order7 = Math.sqrt((3. - 2.*Math.sqrt(6./5.))/7.);
		double u2_order7 = Math.sqrt((3. + 2.*Math.sqrt(6./5.))/7.);
		double w1_order7 = (18. + Math.sqrt(30.))/36.;
		double w2_order7 = (18. - Math.sqrt(30.))/36.;
		double[][] U1_2D_order7 = {{-u1_order7, -u1_order7, -u1_order7, -u1_order7},
		  		                   { u1_order7,  u1_order7,  u1_order7,  u1_order7},
		  		                   {-u2_order7, -u2_order7, -u2_order7, -u2_order7},
		  		                   { u2_order7,  u2_order7,  u2_order7,  u2_order7}};
		double[][] U2_2D_order7 = {{-u1_order7,  u1_order7, -u2_order7, u2_order7},
		  		                   {-u1_order7,  u1_order7, -u2_order7, u2_order7},
		  		                   {-u1_order7,  u1_order7, -u2_order7, u2_order7},
		  		                   {-u1_order7,  u1_order7, -u2_order7, u2_order7}};
		double[][] W_2D_order7 = {{w1_order7*w1_order7, w1_order7*w1_order7, w1_order7*w2_order7, w1_order7*w2_order7},
		  		                  {w1_order7*w1_order7, w1_order7*w1_order7, w1_order7*w2_order7, w1_order7*w2_order7},
		  		                  {w2_order7*w1_order7, w2_order7*w1_order7, w2_order7*w2_order7, w2_order7*w2_order7},
		  		                  {w2_order7*w1_order7, w2_order7*w1_order7, w2_order7*w2_order7, w2_order7*w2_order7}};

	    int Lx = U1_2D_order7.length;
	    int Ly = U1_2D_order7[0].length;
	    
		double sum = 0;
		for (int i=0; i<Lx; i++) {
			for (int j=0; j<Ly; j++) {
				sum += W_2D_order7[i][j]*integrand.function(U1_2D_order7[i][j], U2_2D_order7[i][j]);
			}
		}
		return sum;
	}
	    
	// compute the integrals numerically with a 2D 2x3 point quadrature, which is
	// exact up to order 5 on each direction
	static double quadrature_integral_2D_order5(Integrand2D integrand)
	{
		double u_order5 = Math.sqrt(3./5.);
		double wu_order5 = 5./9.;
		double w0_order5 = 8./9.;
		double[][] U1_2D_order5 = {{-u_order5, 0., u_order5},
                				   {-u_order5, 0., u_order5},
                				   {-u_order5, 0., u_order5}};
		double[][] U2_2D_order5 = {{-u_order5, -u_order5, -u_order5},
				                   {       0.,         0.,       0.},
				                   { u_order5,  u_order5,  u_order5}};
		double[][] W_2D_order5 = {{ wu_order5*wu_order5, w0_order5*wu_order5, wu_order5*wu_order5},
		  		                  { wu_order5*w0_order5, w0_order5*w0_order5, wu_order5*w0_order5},
				                  { wu_order5*wu_order5, w0_order5*wu_order5, wu_order5*wu_order5}};
	    int Lx = U1_2D_order5.length;
	    int Ly = U1_2D_order5[0].length;
	    
		double sum = 0;
		for (int i=0; i<Lx; i++) {
			for (int j=0; j<Ly; j++) {
				sum += W_2D_order5[i][j]*integrand.function(U1_2D_order5[i][j], U2_2D_order5[i][j]);
			}
		}
		return sum;
	}
	
	// compute the integrals numerically with a 2D 2x2 point quadrature, which is
	// exact up to order 3 on each direction
	static double quadrature_integral_2D_order3(Integrand2D integrand)
	{
		double u_order3 = 1./Math.sqrt(3.);
		double[][] U1_2D_order3 = {{-u_order3,  u_order3},
				                   {-u_order3,  u_order3}};
		double[][] U2_2D_order3 = {{-u_order3, -u_order3},
								   { u_order3,  u_order3}};
	    int Lx = U1_2D_order3.length;
	    int Ly = U1_2D_order3[0].length;
	    
		double sum = 0;
		for (int i=0; i<Lx; i++) {
			for (int j=0; j<Ly; j++) {
				sum += integrand.function(U1_2D_order3[i][j], U2_2D_order3[i][j]);
			}
		}
		return sum;
	}
}

