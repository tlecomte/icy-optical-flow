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

package plugins.tlecomte.opticalflow;

import plugins.tlecomte.opticalflow.BilinearMesh;

public class PrecomputeQuadratures { 
	public static double[] precompute_phi_quadratures(final BilinearMesh mesh, final double detJe) {
		System.out.println("Precompute stiffness matrix phi quadratures");
	    
	    int l = mesh.nodes_per_element();
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<l; j++) {
	        int[] ajbj = mesh.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        
	        Integrand2DFunction integrand = new Integrand2DFunction() {
				public double function(double x, double y) {
					return detJe*mesh.Nbar(aj, x)*mesh.Nbar(bj, y);
				}
			};
	        
	        A_pre[j] = Quadratures.quadrature_integral_2D_order3(integrand);
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_phiphi_quadratures(final BilinearMesh mesh, final double detJe, double xi_x, double eta_y) {
		System.out.println("Precompute stiffness matrix phi_j*phi_k quadratures");
	    
		int m = mesh.nodes_per_element();
	    int l = (int) Math.pow(m, 2);
	    double[] A_pre = new double[l];

	    for (int j=0; j<m; j++) {
	        int[] ajbj = mesh.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<m; k++) {
		        int[] akbk = mesh.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];
		        
		        Integrand2DFunction integrand = new Integrand2DFunction() {
    				public double function(double x, double y) {
    					return detJe*mesh.Nbar(aj, x)*mesh.Nbar(bj, y)*mesh.Nbar(ak, x)*mesh.Nbar(bk, y);
    				}
    			};
	            
	    	    if (mesh.order() == 1) // bilinear
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	    	    else // biquad
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order5(integrand);
	        }
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_dxphidxphi_quadratures(final BilinearMesh mesh, final double detJe, final double xi_x, double eta_y) {
		System.out.println("Precompute stiffness matrix dxphi*dxphi quadratures");

		int m = mesh.nodes_per_element();
	    int l = (int) Math.pow(m, 2);
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<m; j++) {
	        int[] ajbj = mesh.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<m; k++) {
		        int[] akbk = mesh.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];
		        
		        Integrand2DFunction integrand = new Integrand2DFunction() {
    				public double function(double x, double y) {
    					return detJe*xi_x*mesh.Nbar_prime(aj, x)*mesh.Nbar(bj, y)*xi_x*mesh.Nbar_prime(ak, x)*mesh.Nbar(bk, y);
    				}
    			};

	    	    if (mesh.order() == 1) // bilinear
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	    	    else // biquad
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order5(integrand);
	        }
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_dyphidyphi_quadratures(final BilinearMesh mesh, final double detJe, double xi_x, final double eta_y) {
		System.out.println("Precompute stiffness matrix dyphi*dyphi quadratures");

		int m = mesh.nodes_per_element();	    
	    int l = (int) Math.pow(m, 2);
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<m; j++) {
	        int[] ajbj = mesh.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<m; k++) {
		        int[] akbk = mesh.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];

		        Integrand2DFunction integrand = new Integrand2DFunction() {
    				public double function(double x, double y) {
    					return detJe*mesh.Nbar(aj, x)*eta_y*mesh.Nbar_prime(bj, y)*mesh.Nbar(ak, x)*eta_y*mesh.Nbar_prime(bk, y);
    				}
    			};
		        
	    	    if (mesh.order() == 1) // bilinear
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	    	    else // biquad
	    	        A_pre[j*m + k] = Quadratures.quadrature_integral_2D_order5(integrand);
	        }
	    }
	
	    return A_pre;
	}
}