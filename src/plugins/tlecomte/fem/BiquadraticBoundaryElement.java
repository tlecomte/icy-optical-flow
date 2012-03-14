package plugins.tlecomte.fem;

public class BiquadraticBoundaryElement extends BoundaryElement {
	BiquadraticBoundaryElement(int i, BoundaryNode node0,
									  BoundaryNode node1,
									  BoundaryNode node2,
            						  double normalX,
            						  double normalY) {
		nodes = new BoundaryNode[3];
		index = i;
		nodes[0] = node0;
		nodes[1] = node1;
		nodes[2] = node2;
		this.normalX = normalX;
		this.normalY = normalY;
	}
}
