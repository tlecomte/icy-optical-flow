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

import icy.image.IcyBufferedImage;
import icy.math.ArrayMath;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.filtering.Kernels1D;
import plugins.adufour.filtering.Convolution1D;
import plugins.adufour.filtering.FilterToolbox.Axis;
import plugins.tlecomte.opticalflow.PrecomputeQuadratures;

import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.sparse.*;

public class HornSchunk {
	static FlexCompRowMatrix assemble_stiffness_matrix(BilinearMesh mesh, double[] dxI2, double[] dyI2, double alpha, double detJe, double xi_x, double eta_y) {
	    double[] A_phiphi = PrecomputeQuadratures.precompute_phiphi_quadratures(mesh, detJe, xi_x, eta_y);
	    double[] A_dxphidxphi = PrecomputeQuadratures.precompute_dxphidxphi_quadratures(mesh, detJe, xi_x, eta_y);
	    double[] A_dyphidyphi = PrecomputeQuadratures.precompute_dyphidyphi_quadratures(mesh, detJe, xi_x, eta_y);

	    int Nnode = mesh.nodes.length; // number of nodes
	    int n = 8*mesh.elements.length*((int) Math.pow(mesh.nodes_per_element(), 2));
	    int[] row = new int[n];
	    int[] col = new int[n];
	    double[] data_coo = new double[n];
	    // build the matrix by summing over all the elements contributions
	    int j = 0;
	    for (int i=0; i<mesh.elements.length; i++) {
	        j = StiffnessMatrix.build_A_HS_COO(j, row, col, data_coo, dxI2[i], dyI2[i], alpha, Nnode, mesh.elements[i], A_phiphi, A_dxphidxphi, A_dyphidyphi);
	    }
	        
	    FlexCompRowMatrix A = new FlexCompRowMatrix(2*Nnode /* numRows */, 2*Nnode /* numColumns */);

	    for (int i=0; i<n; i++) {
	    	A.add(row[i], col[i], data_coo[i]); // isn't it going to be slow ?
	    }
	    
	    int M = 2*Nnode;
	    int N = 2*Nnode;
	    int[] indptr = new int[M + 1];
	    int nnz = data_coo.length;
	    int[] indices = new int[nnz];
	    // data = np.empty(self.nnz, dtype=upcast(self.dtype))
	    double [] data_csr = new double[nnz];

        coo_tocsr(M, N, nnz, row, col, data_coo, indptr, indices, data_csr);

        // A.sum_duplicates()
        // A.prune
   
//            A = csr_matrix((data, indices, indptr), shape=self.shape)

	        
        // CSR sparse format, better suited to solver algorithms
        CompRowMatrix ACSR = CompRowMatrix(M, N, int[][] nz);
        
	    A.scale(0.5);
	    
	    return A;
	}

    /*
    * Compute B = A for COO matrix A, CSR matrix B
    *
    *
    * Input Arguments:
    * I n_row - number of rows in A
    * I n_col - number of columns in A
    * I nnz - number of nonzeros in A
    * I Ai[nnz(A)] - row indices
    * I Aj[nnz(A)] - column indices
    * T Ax[nnz(A)] - nonzeros
    * Output Arguments:
    * I Bp - row pointer
    * I Bj - column indices
    * T Bx - nonzeros
    *
    * Note:
    * Output arrays Bp, Bj, and Bx must be preallocated
    *
    * Note:
    * Input: row and column indices *are not* assumed to be ordered
    *
    * Note: duplicate entries are carried over to the CSR represention
    *
    * Complexity: Linear. Specifically O(nnz(A) + max(n_row,n_col))
    *
    */

    static void coo_tocsr(int n_row, int n_col, int nnz, int[] Ai, int[] Aj, double[] Ax, int[] Bp, int[] Bj, double[] Bx)
    {
        //compute number of non-zero entries per row of A
        for (int i=0; i<n_row; i++) {
        	Bp[i] = 0;
        }

        for (int n = 0; n < nnz; n++){
            Bp[Ai[n]]++;
        }

        //cumsum the nnz per row to get Bp[]
        for(int i = 0, cumsum = 0; i < n_row; i++){
            int temp = Bp[i];
            Bp[i] = cumsum;
            cumsum += temp;
        }
        Bp[n_row] = nnz;

        //write Aj,Ax into Bj,Bx
        for(int n = 0; n < nnz; n++){
            int row = Ai[n];
            int dest = Bp[row];

            Bj[dest] = Aj[n];
            Bx[dest] = Ax[n];

            Bp[row]++;
        }

        for(int i = 0, last = 0; i <= n_row; i++){
            int temp = Bp[i];
            Bp[i] = last;
            last = temp;
        }

        //now Bp,Bj,Bx form a CSR representation (with possible duplicates)
    }

    /*
    * Sum together duplicate column entries in each row of CSR matrix A
    *
    *
    * Input Arguments:
    * I n_row - number of rows in A (and B)
    * I n_col - number of columns in A (and B)
    * I Ap[n_row+1] - row pointer
    * I Aj[nnz(A)] - column indices
    * T Ax[nnz(A)] - nonzeros
    *
    * Note:
    * The column indicies within each row must be in sorted order.
    * Explicit zeros are retained.
    * Ap, Aj, and Ax will be modified *inplace*
    *
    */
    void csr_sum_duplicates(int n_row, int n_col, int[] Ap, int[] Aj, double[] Ax)
    {
        int nnz = 0;
        int row_end = 0;
        for(int i = 0; i < n_row; i++){
            int jj = row_end;
            row_end = Ap[i+1];
            while( jj < row_end ){
                int j = Aj[jj];
                double x = Ax[jj];
                jj++;
                while( jj < row_end && Aj[jj] == j ){
                    x += Ax[jj];
                    jj++;
                }
                Aj[nnz] = j;
                Ax[nnz] = x;
                nnz++;
            }
            Ap[i+1] = nnz;
        }
    }
	  
    void prune() {
        //Remove empty space after all non-zero elements.	
    	nnz = indptr[-1];
        data = data[:nnz];
        indices = indices[:nnz];
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
		System.out.println("Assembling load vector");
	    
	    double[] b_phi = PrecomputeQuadratures.precompute_phi_quadratures(mesh, detJe);
		
	    int N = mesh.nodes.length; // number of nodes  
	    double[] b = new double[2*N];
	    
	    for (int i=0; i<mesh.elements.length; i++) {
	        build_b(b, dxI2[i], dyI2[i], It[i], N, mesh.elements[i], b_phi);
	    }
	    
	    ArrayMath.multiply(b, 0.5);
	
	    return b;
	}
	
	static double J(Matrix A, Vector b, Vector u, Vector It) {
		DenseVector Au = new DenseVector(u.size());
	    return u.dot(A.mult(u, Au)) + b.dot(u) + 0.5*It.dot(It);
	}
	
	public static double[] HornSchunkAlgorithm(double[] I1, double[] I2, int Nx, int Ny, double alpha, double epsilon) {
		System.out.println("Main gradient algorithm");
	   
	    double Lx = Nx;
	    double Ly = Ny;
	    double Dx = Lx/Nx; // Dx == 1.
	    double Dy = Ly/Ny; // Dy == 1.
	    
	    double detJe = Dx*Dy/4.; 
	    double xi_x = 2./Dx;
	    double eta_y = 2./Dy;
	
	    //double[] dxI2, dyI2 = image_first_derivatives(I2);
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
	    FlexCompRowMatrix A = assemble_stiffness_matrix(mesh, dxI2, dyI2, alpha, detJe, xi_x, eta_y);
	    System.out.println("Assemble the load vector");
	    double[] b = assemble_load_vector(mesh, dxI2, dyI2, It, detJe, xi_x, eta_y);
	    
	    // initialize the solutions
	    double[] u = new double[2*mesh.nodes.length]; // velocities

	    // convert to DenseVector so that we can pass them to the MTJ library
	    DenseVector b_d = new DenseVector(b);
	    DenseVector u_d = new DenseVector(u);
	    DenseVector It_d = new DenseVector(It);
	    
	    double Jnew = J(A, b_d, u_d, It_d);
	    System.out.println("Zero-velocity energy: " + Jnew);
	    
	    System.out.println("Starting the solver");
 
	    System.out.println("Diagonal preconditioner");
	    //    #print np.all(np.abs(A.diagonal()) > 1e-15) 
	    //    M = spsparse.spdiags(1./A.diagonal(), 0, 2*N, 2*N)
	    //    M = spsparse.linalg.aslinearoperator(M)
	    //    u, info = spsparse.linalg.cg(A, -b/2., tol=epsilon, M=M); print info
	    //M = spsparse.spdiags(1./A.diagonal(), 0, 2*N, 2*N);
	    //M = spsparse.linalg.aslinearoperator(M);
	    //u = spsparse.linalg.cg(A, -b/2., tol=epsilon, M=M);
	    
	    CG solver = new CG(b_d);
	    
	    // Set the number of iterations, and report on the progress
	    solver.setIterationMonitor(new ResidueIterationMonitor(epsilon));
	    solver.getIterationMonitor().setIterationReporter(new OutputIterationReporter());
	    
	    try {
	    	solver.solve(A, b_d, u_d);	    	
	    }
	    catch (IterativeSolverNotConvergedException e) {
	    	System.err.println("Caught solver exception: " + e.getMessage());
	    }
	    
	    System.out.println("System solved");
	    
	    Jnew = J(A, b_d, u_d, It_d); 
	    System.out.println("Final energy: " + Jnew);
	    
	    return u_d.getData();
	
	    //u1 = u[:N];
	    //u2 = u[N:];
	    //return (np.reshape(u1, (mesh.N_node_y, mesh.N_node_x)),
	    //        np.reshape(u2, (mesh.N_node_y, mesh.N_node_x)));
	    
//		// Put the data in the output image.
//		Array1DUtil.doubleArrayToArray( outputDataBuffer[1], outputImageData);
//		
//		Sequence uvNormSequence = new Sequence();
//      	uvNormSequence.setName("u v velocities Norm");
//      	
//        // create the image object
//        IcyBufferedImage uvNormImage = new IcyBufferedImage(w, h, 1, TypeUtil.TYPE_DOUBLE);
//        Object uvNormImageData = uvNormImage.getDataXY(0);
//		// Put the velocities data in output images.
//		Array1DUtil.doubleArrayToArray(uvNormDataBuffer, uvNormImageData);
//    	
//        // notify to icy that data has changed to refresh internal state and display
//        uvNormImage.dataChanged();
//
//        // add the new images to the sequences at a new time point
//        uvNormSequence.setImage(uvNormSequence.getSizeT(), 0 /*z*/, uvNormImage);
	}
	
	  /**
	   * Iteration monitor. Stops the iteration after the residue value relative to its initial
	   * value is smaller than epsilon.
	   */  
	  private static class ResidueIterationMonitor extends AbstractIterationMonitor {
	
	    private double epsilon;
	    private double r0;
	
	    ResidueIterationMonitor(double epsilon) {
	      this.epsilon = epsilon;
	    }

	    protected boolean convergedI(double r, Vector x) throws IterativeSolverNotConvergedException {
		      return convergedI(r);
		    }
	    
	    protected boolean convergedI(double r) throws IterativeSolverNotConvergedException {
	    	if (isFirst()) {
	    		r0 = r;
	    	}
	        return r/r0 < epsilon;
	    }
	  }
}