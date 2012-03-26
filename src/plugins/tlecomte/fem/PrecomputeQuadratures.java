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

package plugins.tlecomte.fem;

import icy.math.ArrayMath;

public class PrecomputeQuadratures { 
	public static double[] precompute_phi_quadratures(final Mesh mesh, final double detJe) {
		System.out.println("Precompute stiffness matrix phi quadratures");
	    
	    int l = mesh.nodes_per_element();
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<l; j++) {
	        int[] ajbj = mesh.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        
	        Integrand2D integrand = new Integrand2D() {
				public double function(double x, double y) {
					return detJe*mesh.Nbar(aj, x)*mesh.Nbar(bj, y);
				}
			};
	        
	        A_pre[j] = Quadratures.quadrature_integral_2D_order3(integrand);
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_phipsi_quadratures(final Mesh meshPhi, final Mesh meshPsi, final double detJe, double xi_x, double eta_y) {
		System.out.println("Precompute stiffness matrix phi_j*phi_k quadratures");
	    
	    int l = meshPhi.nodes_per_element()*meshPsi.nodes_per_element();
	    double[] A_pre = new double[l];

	    for (int j=0; j<meshPhi.nodes_per_element(); j++) {
	        int[] ajbj = meshPhi.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<meshPsi.nodes_per_element(); k++) {
		        int[] akbk = meshPsi.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];
		        
		        Integrand2D integrand = new Integrand2D() {
    				public double function(double x, double y) {
    					return detJe*meshPhi.Nbar(aj, x)*meshPhi.Nbar(bj, y)*meshPsi.Nbar(ak, x)*meshPsi.Nbar(bk, y);
    				}
    			};
	            
	    	    if (meshPsi.order() + meshPhi.order() <= 3) // lin-lin or quad-lin
	    	        A_pre[j*meshPsi.nodes_per_element() + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	    	    else // biquad
	    	        A_pre[j*meshPsi.nodes_per_element() + k] = Quadratures.quadrature_integral_2D_order5(integrand);
	        }
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_phiphi_quadratures(final Mesh mesh, final double detJe, double xi_x, double eta_y) {
	    return precompute_phipsi_quadratures(mesh, mesh, detJe, xi_x, eta_y);
	}
	
	public static double[] precompute_dxphidxphi_quadratures(final Mesh mesh, final double detJe, final double xi_x, double eta_y) {
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
		        
		        Integrand2D integrand = new Integrand2D() {
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
	
	public static double[] precompute_dyphidyphi_quadratures(final Mesh mesh, final double detJe, double xi_x, final double eta_y) {
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

		        Integrand2D integrand = new Integrand2D() {
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
	
	public static double[] precompute_2D_gradphigradphi_quadratures(final Mesh mesh, final double detJe, double xi_x, final double eta_y) {
	    double[] dx = precompute_dxphidxphi_quadratures(mesh, detJe, xi_x, eta_y);
	    double[] dy = precompute_dyphidyphi_quadratures(mesh, detJe, xi_x, eta_y);
	    return ArrayMath.add(dx, dy);
	}
	
	public static double[] precompute_dxphipsi_quadratures(final BilinearMesh meshLin, final BiquadraticMesh meshQuad, final double detJe, final double xi_x, final double eta_y) {
		System.out.println("Precompute stiffness matrix dxphi*psi quadratures");

	    int l = meshLin.nodes_per_element()*meshQuad.nodes_per_element();
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<meshLin.nodes_per_element(); j++) {
	        int[] ajbj = meshLin.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<meshQuad.nodes_per_element(); k++) {
		        int[] akbk = meshQuad.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];

		        Integrand2D integrand = new Integrand2D() {
    				public double function(double x, double y) {
    					return detJe*meshLin.Nbar(aj, x)*meshLin.Nbar(bj, y)*xi_x*meshQuad.Nbar_prime(ak, x)*meshQuad.Nbar(bk, y);
    				}
    			};
		        
	    	    A_pre[j*meshQuad.nodes_per_element() + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	        }
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_dyphipsi_quadratures(final BilinearMesh meshLin, final BiquadraticMesh meshQuad, final double detJe, final double xi_x, final double eta_y) {
		System.out.println("Precompute stiffness matrix dxphi*psi quadratures");

	    int l = meshLin.nodes_per_element()*meshQuad.nodes_per_element();
	    double[] A_pre = new double[l];
	
	    for (int j=0; j<meshLin.nodes_per_element(); j++) {
	        int[] ajbj = meshLin.canonical_node_coord_2D(j);
	        final int aj = ajbj[0];
	        final int bj = ajbj[1];
	        for (int k=0; k<meshQuad.nodes_per_element(); k++) {
		        int[] akbk = meshQuad.canonical_node_coord_2D(k);
		        final int ak = akbk[0];
		        final int bk = akbk[1];

		        Integrand2D integrand = new Integrand2D() {
    				public double function(double x, double y) {
    					return detJe*meshLin.Nbar(aj, x)*meshLin.Nbar(bj, y)*meshQuad.Nbar(ak, x)*eta_y*meshQuad.Nbar_prime(bk, y);
    				}
    			};
		        
	    	    A_pre[j*meshQuad.nodes_per_element() + k] = Quadratures.quadrature_integral_2D_order3(integrand);
	        }
	    }
	
	    return A_pre;
	}
	
	public static double[] precompute_1D_dxphidxphi_quadratures(final Mesh mesh, final double detJe, final double xi_x, final double eta_y) {
		System.out.println("Precompute 1D stiffness matrix dxphi*dxphi quadratures");

	    int l = (int) Math.pow(mesh.nodesPerEdge(), 2);
	    double[] A_pre = new double[l];
	
	    final double detJe1D = 1./xi_x;
	    final double xi1D = xi_x;
	    
	    for (int j=0; j<mesh.nodesPerEdge(); j++) {
	        for (int k=0; k<mesh.nodesPerEdge(); k++) {
		        final int aj = mesh.canonical_node_coord_1D(j);
		        final int ak = mesh.canonical_node_coord_1D(k);

		        Integrand1D integrand = new Integrand1D() {
    				public double function(double x) {
    					return detJe1D*xi1D*mesh.Nbar_prime(aj, x)*xi1D*mesh.Nbar_prime(ak, x);
    				}
    			};
		        
	    	    A_pre[j*mesh.nodesPerEdge() + k] = Quadratures.quadrature_integral_1D_order3(integrand);
	        }
	    }
	
	    return A_pre;
	}

	public static double[] precompute_1D_dyphidyphi_quadratures(final Mesh mesh, final double detJe, final double xi_x, final double eta_y) {
		System.out.println("Precompute 1D stiffness matrix dyphi*dyphi quadratures");

	    int l = (int) Math.pow(mesh.nodesPerEdge(), 2);
	    double[] A_pre = new double[l];
	
	    final double detJe1D = 1./eta_y;
	    final double xi1D = eta_y;
	    
	    for (int j=0; j<mesh.nodesPerEdge(); j++) {
	        for (int k=0; k<mesh.nodesPerEdge(); k++) {
		        final int aj = mesh.canonical_node_coord_1D(j);
		        final int ak = mesh.canonical_node_coord_1D(k);

		        Integrand1D integrand = new Integrand1D() {
    				public double function(double x) {
    					return detJe1D*xi1D*mesh.Nbar_prime(aj, x)*xi1D*mesh.Nbar_prime(ak, x);
    				}
    			};
		        
	    	    A_pre[j*mesh.nodesPerEdge() + k] = Quadratures.quadrature_integral_1D_order3(integrand);
	        }
	    }
	
	    return A_pre;
	}
}