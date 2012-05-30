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

public class BiquadraticElement extends Element {
	BiquadraticElement(int i, Node node0,
							  Node node1,
							  Node node2,
							  Node node3,
							  Node node4,
							  Node node5,
							  Node node6,
							  Node node7,
							  Node node8) {
		index = i;
		nodes = new Node[9];
		nodes[0] = node0;
		nodes[1] = node1;
		nodes[2] = node2;
		nodes[3] = node3;
		nodes[4] = node4;
		nodes[5] = node5;
		nodes[6] = node6;
		nodes[7] = node7;
		nodes[8] = node8;
	}
}