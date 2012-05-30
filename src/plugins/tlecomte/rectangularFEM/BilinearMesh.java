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

import plugins.tlecomte.rectangularFEM.FemFunctions;

public class BilinearMesh extends Mesh {

	public BilinearMesh(int Nx, int Ny) {   
        N_node_x = Nx + 1;
        N_node_y = Ny + 1;
        int node_number = N_node_x*N_node_y;
        
        // data structure for the rectangles
        
        // node list - three columns : Node number, x coordinate, y coordinate
        nodes = new Node[node_number];
        int i =0;
        int j = 0;
    	for (int y=0; y<N_node_y; y++) {
            for (int x=0; x<N_node_x; x++) {
        		nodes[i] = new Node(i, x, y);
        		
                if (x==0 || x==N_node_x-1 || y==0 || y==N_node_y-1) {
                    // fixed node
                	nodes[i].setPointer(-1);
                } else {
                    // free node
                	nodes[i].setPointer(j);
                    j += 1;
                }
        		
        		i++;
        	}
        }
        
        freeNodes = new Node[j];
        j = 0;
        for (i=0; i<node_number; i++) {
        	int pointer = nodes[i].pointer;
            if (pointer != -1) {
                // free node
                freeNodes[j] = nodes[i];
                j += 1;
            }
        }
        
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
        
        int leftEdgeNumber = Ny; // as many edges as elements
        int rightEdgeNumber = Ny;
        int topEdgeNumber = Nx;
        int bottomEdgeNumber = Nx;
        
        int edgeNumber = (leftEdgeNumber + rightEdgeNumber
                             + topEdgeNumber + bottomEdgeNumber);
        
        boundaryNodes = new BoundaryNode[edgeNumber];

        j = 0;
        for (i=0; i<bottomEdgeNumber; i++) {
        	boundaryNodes[j] = new BoundaryNode(j, nodes[i]);
        	j += 1; 
        }

        for (i=0; i<rightEdgeNumber; i++) {
        	boundaryNodes[j] = new BoundaryNode(j, nodes[(i+1)*(Nx+1) - 1]);
        	j += 1; 
        }

        // reversed order for counter-clockwise order !
        for (i=topEdgeNumber-1; i>=0; i--) {
        	boundaryNodes[j] = new BoundaryNode(j, nodes[(Nx+1)*Ny + i + 1]);
        	j += 1; 
        }           

        // reversed order for counter-clockwise order !
        for (i=leftEdgeNumber-1; i>=0; i--) {
        	boundaryNodes[j] = new BoundaryNode(j, nodes[(i+1)*(Nx+1)]);
        	j += 1; 
        }
        
        edges = new BoundaryElement[edgeNumber];

        j = 0;
        for (i=0; i<bottomEdgeNumber; i++) {
        	double normalX = 0.;
        	double normalY = -1.;
        	edges[j] = new BilinearBoundaryElement(j, boundaryNodes[j], boundaryNodes[j+1], normalX, normalY);
        	j += 1; 
        }

        for (i=0; i<rightEdgeNumber; i++) {
        	double normalX = 1.;
        	double normalY = 0.;
        	edges[j] = new BilinearBoundaryElement(j, boundaryNodes[j], boundaryNodes[j+1], normalX, normalY);
        	j += 1; 
        }

        for (i=0; i<topEdgeNumber; i++) {
        	double normalX = 0.;
        	double normalY = 1.;
        	edges[j] = new BilinearBoundaryElement(j, boundaryNodes[j], boundaryNodes[j+1], normalX, normalY);
        	j += 1; 
        }
        
        for (i=0; i<leftEdgeNumber; i++) {
        	double normalX = -1.;
        	double normalY = 0.;
        	edges[j] = new BilinearBoundaryElement(j, boundaryNodes[j], boundaryNodes[(j+1)%edgeNumber], normalX, normalY);
        	j += 1; 
        }
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