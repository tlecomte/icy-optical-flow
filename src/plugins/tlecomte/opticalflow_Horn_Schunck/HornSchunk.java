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

import icy.math.ArrayMath;
import plugins.adufour.filtering.Kernels1D;
import plugins.adufour.filtering.Convolution1D;
import plugins.adufour.filtering.FilterToolbox.Axis;
import plugins.tlecomte.fem.PrecomputeQuadratures;
import plugins.tlecomte.fem.BilinearMesh;
import plugins.tlecomte.fem.Element;
import plugins.tlecomte.fem.Node;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.AbstractDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.*;

public class HornSchunk {
	static SparseRCDoubleMatrix2D assemble_stiffness_matrix(BilinearMesh mesh, double[] dxI2, double[] dyI2, double alpha, double detJe, double xi_x, double eta_y) {
	    double[] A_phiphi = PrecomputeQuadratures.precompute_phiphi_quadratures(mesh, detJe, xi_x, eta_y);
	    double[] A_dxphidxphi = PrecomputeQuadratures.precompute_dxphidxphi_quadratures(mesh, detJe, xi_x, eta_y);
	    double[] A_dyphidyphi = PrecomputeQuadratures.precompute_dyphidyphi_quadratures(mesh, detJe, xi_x, eta_y);

	    int N = mesh.nodes.length; // number of nodes
	    int n = 8*mesh.elements.length*((int) Math.pow(mesh.nodes_per_element(), 2));
	    int[] row = new int[n];
	    int[] col = new int[n];
	    double[] data_coo = new double[n];
	    // build the matrix by summing over all the elements contributions
	    int j = 0;
	    for (int i=0; i<mesh.elements.length; i++) {
	        j = StiffnessMatrix.build_A_HS_COO(j, row, col, data_coo, dxI2[i], dyI2[i], alpha, N, mesh.elements[i], A_phiphi, A_dxphidxphi, A_dyphidyphi);
	    }
	    
	    data_coo = ArrayMath.multiply(data_coo, 0.5);
	    
	    SparseRCDoubleMatrix2D A = new SparseRCDoubleMatrix2D(2*N, 2*N, row, col, data_coo, true /* remove duplicates */, true /* remove zeroes */, false /* sort column indices */);
	    
	    return A;
	}
    
	static void build_b(double[] b, double dxI2_e, double dyI2_e, double It_e, int N, Element element, double[] A_phi)
	{
		int i = 0;
	    for (Node node : element.nodes) {
	            b[node.index]     += 2.*It_e*dxI2_e*A_phi[i];
	            b[N + node.index] += 2.*It_e*dyI2_e*A_phi[i];
	            i += 1;
	    }
	}
	
	static double[] assemble_load_vector(BilinearMesh mesh, double[] dxI2, double[] dyI2, double[] It, double detJe, double xi_x, double eta_y) {   
	    double[] b_phi = PrecomputeQuadratures.precompute_phi_quadratures(mesh, detJe);
		
	    int N = mesh.nodes.length; // number of nodes  
	    double[] b = new double[2*N];
	    
	    for (int i=0; i<mesh.elements.length; i++) {
	        build_b(b, dxI2[i], dyI2[i], It[i], N, mesh.elements[i], b_phi);
	    }
	    
	    b = ArrayMath.multiply(b, 0.5);
	
	    return b;
	}
	
	static double J(SparseRCDoubleMatrix2D A, DenseDoubleMatrix1D b, DenseDoubleMatrix1D u, DenseDoubleMatrix1D It) {
		DenseDoubleMatrix1D Au = new DenseDoubleMatrix1D((int) u.size());
	    return u.zDotProduct(A.zMult(u, Au)) + b.zDotProduct(u) + 0.5*It.zDotProduct(It);
	}
	
	public static double[] HornSchunkAlgorithm(double[] I1, double[] I2, int Nx, int Ny, double alpha, int maxiter, double epsilon) {
		System.out.println("Main gradient algorithm");

	    long t0 = System.currentTimeMillis();
		
	    double Lx = Nx;
	    double Ly = Ny;
	    double Dx = Lx/Nx; // Dx == 1.
	    double Dy = Ly/Ny; // Dy == 1.
	    
	    double detJe = Dx*Dy/4.; 
	    double xi_x = 2./Dx;
	    double eta_y = 2./Dy;
	
	    double[] derivKernel = Kernels1D.GRADIENT.getData();
	    double[][] I2z = {I2};
	    double[][] dxI2z = new double[1][Nx*Ny];
	    double[][] dyI2z = new double[1][Nx*Ny];
	    Convolution1D.convolve1D(I2z, dxI2z, Nx, Ny, derivKernel, Axis.X);
	    Convolution1D.convolve1D(I2z, dyI2z, Nx, Ny, derivKernel, Axis.Y);
	    double[] dxI2 = dxI2z[0];
	    double[] dyI2 = dyI2z[0];
	
	    double[] It = ArrayMath.subtract(I2, I1);
	    
	    // build a mesh for the bilinear finite element basis
	    System.out.println("Build the mesh");
	    BilinearMesh mesh = new BilinearMesh(Nx, Ny);
	    
	    System.out.println("Assemble the stiffness matrix");
	    SparseRCDoubleMatrix2D A = assemble_stiffness_matrix(mesh, dxI2, dyI2, alpha, detJe, xi_x, eta_y);
	    System.out.println("Assemble the load vector");
	    double[] b = assemble_load_vector(mesh, dxI2, dyI2, It, detJe, xi_x, eta_y);
	    
	    // initialize the solutions
	    double[] u = new double[2*mesh.nodes.length]; // velocities

	    // convert to DenseVector so that we can pass them to the MTJ library
	    DenseDoubleMatrix1D b_d = new DenseDoubleMatrix1D(b);
	    DenseDoubleMatrix1D u_d = new DenseDoubleMatrix1D(u);
	    DenseDoubleMatrix1D It_d = new DenseDoubleMatrix1D(It);
	    
	    double Jinit = J(A, b_d, u_d, It_d);
	    
	    System.out.println("Starting the solver");
 
	    DoubleCG solver = new DoubleCG(b_d);

	    // System.out.println("Preconditioner");
//	    // Optionally create a preconditioner
//        DoubleAMG M = new DoubleAMG();
//	    M.setMatrix(A);
//	    // Attach the preconditioner
//	    solver.setPreconditioner(M);    
	    
	    // Set the number of iterations, and report on the progress
	    ResidueIterationMonitor monitor = new ResidueIterationMonitor(maxiter, epsilon);
	    OutputInterationReporter reporter = new OutputInterationReporter();
	    solver.setIterationMonitor(monitor);
	    solver.getIterationMonitor().setIterationReporter(reporter);
	    
	    b = ArrayMath.multiply(b, -0.5);
	    DenseDoubleMatrix1D b2_d = new DenseDoubleMatrix1D(b);
	    
	    try {
	    	solver.solve(A, b2_d, u_d);	    	
	    }
	    catch (IterativeSolverDoubleNotConvergedException e) {
	    	System.err.println("Caught solver exception: " + e.getMessage());
	    }
	    
	    System.out.format("System solved after %d iterations\n", monitor.iterations());
	    
	    double Jend = J(A, b_d, u_d, It_d); 
	    System.out.printf("Zero-velocity energy = %g, Final energy = %g\n", Jinit, Jend);
	    
	    long t1 = System.currentTimeMillis();

	    double fsec = (t1-t0)/1000.;
	    System.out.printf("Flow computed in %.3g seconds.\n", fsec);
	    
	    return u_d.elements();
	}
	
	  /**
	   * Iteration monitor. Stops the iteration after the residue value relative to its initial
	   * value is smaller than epsilon, or the maximum number of iterations was reached.
	   */  
	  private static class ResidueIterationMonitor extends AbstractDoubleIterationMonitor {
	
	    private double epsilon;
	    private double r0;
	    private int max;
	
	    ResidueIterationMonitor(int max, double epsilon) {
	      this.epsilon = epsilon;
	      this.max = max;
	    }

	    protected boolean convergedI(double r, DoubleMatrix1D x) {
		      return convergedI(r);
		    }
	    
	    protected boolean convergedI(double r) {
	    	if (isFirst()) {
	    		r0 = r;
	    	}
	        return (r/r0 < epsilon) || iterations() >= max;
	    }

		public int getMaxIterations() {
			return max;
		}

		public void setMaxIterations(int max) {
			this.max = max;
		}
	  }

	  /**
	   * Iteration Reporter. Prints the iteration number and the current residue value
	   * to the standard output.
	   */ 
	  private static class OutputInterationReporter implements DoubleIterationReporter {
		private final int iterStep = 40;
		  
		public void monitor(double r, int i) {
			if (i%iterStep == 0) {
				System.out.println(String.format("Iteration #%d, residue %.3g", i, r));
			}
		}

		public void monitor(double r, DoubleMatrix1D x, int i) {
			monitor(r, i);
		}
		  
	  }
}