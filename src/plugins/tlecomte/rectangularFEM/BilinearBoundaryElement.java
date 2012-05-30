package plugins.tlecomte.rectangularFEM;

public class BilinearBoundaryElement extends BoundaryElement{
	BilinearBoundaryElement(int i, BoundaryNode node0,
								   BoundaryNode node1,
								   double normalX,
								   double normalY) {
		nodes = new BoundaryNode[2];
		index = i;
		nodes[0] = node0;
		nodes[1] = node1;
		this.normalX = normalX;
		this.normalY = normalY;
	}
}
