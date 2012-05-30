package plugins.tlecomte.rectangularFEM;

public class BiquadraticMesh extends Mesh {
	
	// fields
    //public nodes_pointers;
    //public free_nodes;
    //public boundary_edges;
    //public boundary_normals;
    //self.boundary_free_nodes = boundary_free_nodes;
    //self.boundary_nodes_pointers_to_free = boundary_nodes_pointers_to_free;
    //self.boundary_nodes_pointers_to_grid = boundary_nodes_pointers_to_grid;
    //self.boundary_elements = boundary_elements;
	
	public BiquadraticMesh(int Nx, int Ny) {   
        N_node_x = 2*Nx + 1;
        N_node_y = 2*Ny + 1;
        int node_number = N_node_x*N_node_y;
        
        // data structure for the rectangles
        
        // node list - three columns : Node number, x coordinate, y coordinate
        nodes = new Node[node_number];
        int i = 0;
        int j = 0;
        for (int ny=0; ny<N_node_y; ny++) {
        	for (int nx=0; nx<N_node_x; nx++) {
        		double x = (double) (nx*Nx) / (double) (N_node_x - 1);
        		double y = (double) (ny*Ny) / (double) (N_node_y - 1);
        		nodes[i] = new Node(i, x, y);
        		
                if (nx==0 || nx==N_node_x-1 || ny==0 || ny==N_node_y-1) {
                    // fixed node
                    nodes[i].setPointer(-1);
                } else {
                    // free node
                	// FIXME the meaning of "pointer" could be more explicit, like pointerToFree
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
        elements = new BiquadraticElement[element_number];
        for (i=0; i<element_number; i++) {
            int row = i/Nx;
            int col = i%Nx;
            elements[i] = new BiquadraticElement(i, nodes[2*row*N_node_x + 2*col],
                    					 			nodes[2*row*N_node_x + 2*col + 1],
                    					 			nodes[2*row*N_node_x + 2*col + 2],
                    					 			nodes[(2*row+1)*N_node_x + 2*col],
                    					 			nodes[(2*row+1)*N_node_x + 2*col + 1],
                    					 			nodes[(2*row+1)*N_node_x + 2*col + 2],
                    					 			nodes[(2*row+2)*N_node_x + 2*col],
                    					 			nodes[(2*row+2)*N_node_x + 2*col + 1],
                    					 			nodes[(2*row+2)*N_node_x + 2*col + 2]);
        }
        
        int leftEdgeNumber = Ny; // as many edges as elements
        int rightEdgeNumber = Ny;
        int topEdgeNumber = Nx;
        int bottomEdgeNumber = Nx;
        
        int edgeNumber = (leftEdgeNumber + rightEdgeNumber
                             + topEdgeNumber + bottomEdgeNumber);
        
        boundaryNodes = new BoundaryNode[2*edgeNumber];

        j = 0;
        for (i=0; i<bottomEdgeNumber; i++) {
        	boundaryNodes[2*j]   = new BoundaryNode(2*j,   nodes[2*i]);
        	boundaryNodes[2*j+1] = new BoundaryNode(2*j+1, nodes[2*i+1]);
        	j += 1; 
        }

        for (i=0; i<rightEdgeNumber; i++) {
        	boundaryNodes[2*j]   = new BoundaryNode(2*j,   nodes[(2*i+1)*(2*Nx+1) - 1]);
        	boundaryNodes[2*j+1] = new BoundaryNode(2*j+1, nodes[(2*i+2)*(2*Nx+1) - 1]);
        	j += 1; 
        }

        // reversed order for counter-clockwise order !
        for (i=topEdgeNumber-1; i>=0; i--) {
        	boundaryNodes[2*j]   = new BoundaryNode(2*j,   nodes[(2*Nx+1)*2*Ny + 2*i + 2]);
        	boundaryNodes[2*j+1] = new BoundaryNode(2*j+1, nodes[(2*Nx+1)*2*Ny + 2*i + 1]);
        	j += 1; 
        }           

        // reversed order for counter-clockwise order !
        for (i=leftEdgeNumber-1; i>=0; i--) {
        	boundaryNodes[2*j]   = new BoundaryNode(2*j,   nodes[(2*i+2)*(2*Nx+1)]);
        	boundaryNodes[2*j+1] = new BoundaryNode(2*j+1, nodes[(2*i+1)*(2*Nx+1)]);
        	j += 1; 
        }

        edges = new BoundaryElement[edgeNumber];

        j = 0;
        for (i=0; i<bottomEdgeNumber; i++) {
        	double normalX = 0.;
        	double normalY = -1.;
        	edges[j] = new BiquadraticBoundaryElement(j, boundaryNodes[2*j], boundaryNodes[2*j+1], boundaryNodes[2*j+2], normalX, normalY);
        	j += 1; 
        }

        for (i=0; i<rightEdgeNumber; i++) {
        	double normalX = 1.;
        	double normalY = 0.;
        	edges[j] = new BiquadraticBoundaryElement(j, boundaryNodes[2*j], boundaryNodes[2*j+1], boundaryNodes[2*j+2], normalX, normalY);
        	j += 1; 
        }

        for (i=0; i<topEdgeNumber; i++) {
        	double normalX = 0.;
        	double normalY = 1.;
        	edges[j] = new BiquadraticBoundaryElement(j, boundaryNodes[2*j], boundaryNodes[2*j+1], boundaryNodes[2*j+2], normalX, normalY);
        	j += 1; 
        }
        
        for (i=0; i<leftEdgeNumber; i++) {
        	double normalX = -1.;
        	double normalY = 0.;
        	edges[j] = new BiquadraticBoundaryElement(j, boundaryNodes[2*j], boundaryNodes[2*j+1], boundaryNodes[(2*j+2)%(2*edgeNumber)], normalX, normalY);
        	j += 1; 
        }
	}
        
    // convert a boundary node index between 0 and 1 to 1D coordinates in the quadratic
    // canonical element [-1,1]
    public int canonical_node_coord_1D(int i) {
        int ai = i%3 - 1;
        return ai;
    }
    
    // convert a node index between 0 and 3 to 2D coordinates in the biquadratic
    // canonical element [-1,1]x[-1,1]
    int[] canonical_node_coord_2D(int i) {
        int ai = i%3 - 1;
        int bi = i/3 - 1;
        int[] ret = {ai, bi};
        return ret;
    }

    public double Nbar(int xi, double x) {
        return FemFunctions.Nbar(xi, x);
    }
        
    double[] Nbar_1D(int xi, double[] x) {
        return FemFunctions.Nbar_1D(xi, x);
    }
        
    //def Nbar_2D(xi, x):
    //    return fem_functions.Nbar_2D(xi, x);

    double Nbar_prime(int xi, double x) {
        return FemFunctions.Nbar_prime(xi, x);
    }

    double[] Nbar_prime_1D(int xi, double[] x) {
        return FemFunctions.Nbar_prime_1D(xi, x);
    }

    //def Nbar_prime_2D(xi, x):
    //    return fem_functions.Nbar_prime_2D(xi, x);
    
    int order() {
        return 2; // quadratic mesh
    }
}
