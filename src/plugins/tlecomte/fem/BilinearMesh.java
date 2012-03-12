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

import plugins.tlecomte.fem.FemFunctions;

public class BilinearMesh extends Mesh {
	
	// fields
    //public nodes_pointers;
    //public free_nodes;
    //public boundary_edges;
    //public boundary_normals;
    //self.boundary_free_nodes = boundary_free_nodes;
    //self.boundary_nodes_pointers_to_free = boundary_nodes_pointers_to_free;
    //self.boundary_nodes_pointers_to_grid = boundary_nodes_pointers_to_grid;
    //self.boundary_elements = boundary_elements;
	
	public BilinearMesh(int Nx, int Ny) {   
        N_node_x = Nx + 1;
        N_node_y = Ny + 1;
        	
        // data structure for the rectangles
        
        // node list - three columns : Node number, x coordinate, y coordinate
        nodes = new Node[N_node_x*N_node_y];
        int i =0;
        for (int x=0; x<N_node_x; x++) {
        	for (int y=0; y<N_node_y; y++) {
        		nodes[i] = new Node(i, x, y);
        		i++;
        	}
        }

//        int node_number = N_node_x*N_node_y;
//        
//        int j = 0;
//        for (i=0; i<node_number; i++) {
//        	double[] coords = nodes[i].coords;
//            if (coords[0]==0 || coords[0]==Nx || coords[1]==0 || coords[1]==Ny) {
//                // fixed node
//                nodes_pointers[i] = -1;
//            } else {
//                // free node
//                nodes_pointers[i] = j;
//                free_nodes[j] = i;
//                j += 1;
//            }
//        }
        
        // element list - five columns : Element number, node1, node 2, node3, node4
        int element_number = Nx*Ny;
        elements = new Element[element_number];
        for (i=0; i<element_number; i++) {
            int row = i/Nx;
            int col = i%Nx;
            elements[i] = new BilinearElement(i, nodes[row*N_node_x + col],
                    					         nodes[row*N_node_x + col + 1],
                    					         nodes[(row+1)*N_node_x + col],
                    					         nodes[(row+1)*N_node_x + col + 1]);
        }
        
//        // boundary edges list - three columns : Edge number, node1, node2
//        left_boundary_edge_number = Ny; // as many edges as elements
//        right_boundary_edge_number = Ny;
//        top_boundary_edge_number = Nx;
//        bottom_boundary_edge_number = Nx;
//        
//        boundary_edge_number = (left_boundary_edge_number + right_boundary_edge_number
//                             + top_boundary_edge_number + bottom_boundary_edge_number);
//                             
//        boundary_edges = np.zeros((boundary_edge_number, 3));
//        boundary_normals = np.zeros((boundary_edge_number, 3));
//        boundary_free_nodes = np.zeros(boundary_edge_number);
//        boundary_nodes_pointers_to_free = np.zeros(boundary_edge_number);
//        boundary_nodes_pointers_to_grid = np.zeros(boundary_edge_number);
//        boundary_elements = np.zeros(boundary_edge_number);
//        
//        iterator = xrange(boundary_edge_number).__iter__();
//        k = 0    
//        
//        for i in xrange(bottom_boundary_edge_number):
//            j = iterator.next();
//            boundary_edges[j] = [j, j, j+1];
//            boundary_nodes_pointers_to_grid[j] = i;
//            // free boundary node
//            boundary_nodes_pointers_to_free[j] = k;
//            boundary_free_nodes[k] = j;
//            k += 1;
//            boundary_normals[j] = [j, 0., -1.];
//            boundary_elements[j] = i;
//    
//        for i in xrange(right_boundary_edge_number):
//            j = iterator.next();
//            boundary_edges[j] = [j, j, j+1];
//            boundary_nodes_pointers_to_grid[j] = (i+1)*(Nx+1) - 1;
//            // free boundary node
//            boundary_nodes_pointers_to_free[j] = k;
//            boundary_free_nodes[k] = j;
//            k += 1;
//            boundary_normals[j] = [j, 1., 0.];
//            boundary_elements[j] = Nx-1 + i*Nx;
//            
//        for i in range(top_boundary_edge_number)[::-1]:
//            j = iterator.next();
//            boundary_edges[j] = [j, j, j+1];
//            boundary_nodes_pointers_to_grid[j] = (Nx+1)*Ny + i + 1;
//            // free boundary node
//            boundary_nodes_pointers_to_free[j] = k;
//            boundary_free_nodes[k] = j;
//            k += 1;
//            boundary_normals[j] = [j, 0., 1.];
//            boundary_elements[j] = (Ny-1)*Nx + i;
//        
//        for i in range(left_boundary_edge_number)[::-1]:
//            j = iterator.next();
//            boundary_edges[j] = [j, j, (j+1)%boundary_edge_number];
//            boundary_nodes_pointers_to_grid[j] = (i+1)*(Nx+1);
//            // free boundary node
//            boundary_nodes_pointers_to_free[j] = k;
//            boundary_free_nodes[k] = j;
//            k += 1;
//            boundary_normals[j] = [j, -1., 0.];
//            boundary_elements[j] = i*Nx;
        
        // make a dict out of the various lists
        //nodes = nodes;
        //nodes_pointers = nodes_pointers;
        //free_nodes = free_nodes;
        //elements = elements;
        //self.boundary_edges = boundary_edges;
        //self.boundary_normals = boundary_normals;
        //self.boundary_free_nodes = boundary_free_nodes;
        //self.boundary_nodes_pointers_to_free = boundary_nodes_pointers_to_free;
        //self.boundary_nodes_pointers_to_grid = boundary_nodes_pointers_to_grid;
        //self.boundary_elements = boundary_elements;
	}
        
    // convert a boundary node index between 0 and 1 to 1D coordinates in the linear
    // canonical element [-1,1]
    int canonical_node_coord_1D(int i) {
        int ai = 2*(i%2) - 1;
        return ai;
    }
    
    // convert a node index between 0 and 3 to 2D coordinates in the bilinear
    // canonical element [-1,1]x[-1,1]
    int[] canonical_node_coord_2D(int i) {
        int ai = 2*(i%2) - 1;
        int bi = 2*(i/2) - 1;
        int[] ret = {ai, bi};
        return ret;
    }

    double Nbar(int xi, double x) {
        return FemFunctions.Mbar(xi, x);
    }
        
    double[] Nbar_1D(int xi, double[] x) {
        return FemFunctions.Mbar_1D(xi, x);
    }
        
    //def Nbar_2D(xi, x):
    //    return fem_functions.Mbar_2D(xi, x);

    double Nbar_prime(int xi, double x) {
        return FemFunctions.Mbar_prime(xi, x);
    }

    double[] Nbar_prime_1D(int xi, double[] x) {
        return FemFunctions.Mbar_prime_1D(xi, x);
    }

    //def Nbar_prime_2D(xi, x):
    //    return fem_functions.Mbar_prime_2D(xi, x);
    
    int order() {
        return 1; // linear mesh
    }
}