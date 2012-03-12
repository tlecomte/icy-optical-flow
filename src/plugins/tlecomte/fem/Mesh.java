package plugins.tlecomte.fem;

public abstract class Mesh {
	// fields
    public Node[] nodes;
    public Element[] elements;
    public int N_node_x;
    public int N_node_y;
        
    // convert a boundary node index between 0 and 1 to 1D coordinates in the
    // canonical element [-1,1]
    abstract int canonical_node_coord_1D(int i);
    
    // convert a node index between 0 and 3 to 2D coordinates in the
    // canonical element [-1,1]x[-1,1]
    abstract int[] canonical_node_coord_2D(int i);

    abstract double Nbar(int xi, double x);
        
    abstract double[] Nbar_1D(int xi, double[] x);
        
    //def Nbar_2D(xi, x):
    //    return fem_functions.Nbar_2D(xi, x);

    abstract double Nbar_prime(int xi, double x);

    abstract double[] Nbar_prime_1D(int xi, double[] x);

    //def Nbar_prime_2D(xi, x):
    //    return fem_functions.Nbar_prime_2D(xi, x);
    
    abstract int order();

    public int nodes_per_element() {
        return elements[0].nodes.length;
    }
}
