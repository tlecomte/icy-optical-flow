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

/**
 * 
 * @author Timothee Lecomte
 *
 * This plugin 
 *
 * TODO:
 * Validate Horn-Schunck algorithm on middlebury sample videos
 * Implement a more recent algorithm such as "Optical Flow in Harmony" from Zimmer et al. (2011) 
 * 
 */

package plugins.tlecomte.flowdisplay;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import icy.canvas.IcyCanvas;
import icy.gui.viewer.Viewer;
import icy.painter.Painter;
import icy.sequence.Sequence;

class FlowArrow implements Comparable<FlowArrow>
{
	public FlowArrow( double x , double y )
	{
		this.x = x;
		this.y = y;
	}
	double norme;
	double vx = 0;
	double vy = 0;
	double x,y,angle,thickness;
	Color color;
	
	public int compareTo(FlowArrow o) {
		
		if ( o.angle < this.angle ) return -1;
		return 1;
	}
}

public class VectorFlowPainter implements Painter {
	// List of List of flow arrows, one list per time
	ArrayList<ArrayList<FlowArrow>> flowArrowList = new ArrayList<ArrayList<FlowArrow>>();
	
	boolean hideZeroVelocities = false;
	int resolution = 1;
	
    // clear the arrows list
	public void clear() {
		flowArrowList.clear();
	}
	
	public void update_flow_arrows(double[] u, double[] v, int w, int h) {
   		/* resolution is the number of pixels for the square containing the arrow */
    	
   		ArrayList<FlowArrow> currentImageArrowList = new ArrayList<FlowArrow>();
   		
		// generates arrows by averaging instant speeds over a region
		for (int x = 0; x < w - resolution; x+=resolution)
		{
			for (int y = 0; y < h - resolution; y+=resolution)
			{
				FlowArrow flowarrow = new FlowArrow(x + resolution/2., y + resolution/2.);

				flowarrow.vx = 0;
				flowarrow.vy = 0;
				for (int i=0; i<resolution; i++)
				{
					for (int j=0; j<resolution; j++)
					{
						int id = (y+j)*w + x + i;
						flowarrow.vx += u[id];
						flowarrow.vy += v[id];
					}
				}
				flowarrow.vx /= Math.pow(resolution, 2);
				flowarrow.vy /= Math.pow(resolution, 2);

				currentImageArrowList.add( flowarrow );
			}
		}

		// compute norm and angle
		for ( FlowArrow flowarrow : currentImageArrowList )
		{
			flowarrow.norme = Math.sqrt( flowarrow.vx*flowarrow.vx + flowarrow.vy*flowarrow.vy);
			flowarrow.angle = Math.atan2( flowarrow.vy , flowarrow.vx );
		}

		// add the current arrow map to the global list
		flowArrowList.add( currentImageArrowList );
	}
	
	public void normalize() {	
		// normalize
		double max = 0;
		for (ArrayList<FlowArrow> currentImageArrowList : flowArrowList) {
			for ( FlowArrow flowarrow : currentImageArrowList ) {					
				if ( flowarrow.norme > max ) {
					max = flowarrow.norme ;
				}
			}
		}
		
		if ( max == 0. ) {
			max = 1;
		}
		
		double rapport =  resolution / max ;
		
		for (ArrayList<FlowArrow> currentImageArrowList : flowArrowList) {
			for ( FlowArrow flowarrow : currentImageArrowList )
			{					
				flowarrow.norme *= rapport;
				flowarrow.vx *= rapport;
				flowarrow.vy *= rapport;
	
				double norme1 = flowarrow.norme/resolution;
				// flowarrow.color = new Color( norme1 , 0f , 1f - norme1 );
				flowarrow.color = new Color((float) norme1, 1f , 0f);
				
				// make the arrows more visible by doubling their length 
				flowarrow.norme *= 2;
				flowarrow.vx *= 2;
				flowarrow.vy *= 2;					
			}
		}
	}
	
	public void hideZeroVelocities(boolean hide) {
		hideZeroVelocities = hide;
	}
	
	public void setResolution(int resolution) {
		this.resolution = resolution;
	}
	
   	@Override
   	public void keyPressed(KeyEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void keyReleased(KeyEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void mouseClick(MouseEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void mouseDrag(MouseEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void mouseMove(MouseEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void mousePressed(MouseEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void mouseReleased(MouseEvent e, Point2D imagePoint, IcyCanvas canvas) {
   	}

   	@Override
   	public void paint(Graphics2D g2, Sequence sequence, IcyCanvas canvas) {
   		//int resolution = resolutionSelector.getValue();		
   		//float maxWidth= 4;
   		//float minWidth= 1;

   		//if ( !isEnabled() ) return;
   		//if ( !affectSequenceButton.isSelected() ) return;

   		Viewer viewer = sequence.getFirstViewer();
   		int t = viewer.getT();
   		
   		g2.setColor( Color.yellow );

   		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
   		//BasicStroke stroke = new BasicStroke((float) (minWidth + (maxWidth - minWidth)*flowarrow.norme/resolution));
		BasicStroke stroke = new BasicStroke((float) 0.5);

		t = Math.min(t, sequence.getSizeT() - 2); /* replicate flow on last image */
		
		if (!flowArrowList.isEmpty() && t < flowArrowList.size()) {
			ArrayList<FlowArrow> currentImageArrowList = flowArrowList.get(t);
			
	   		for ( FlowArrow flowarrow : currentImageArrowList )
	   		{
	   			if (hideZeroVelocities) 
	   				if (flowarrow.norme < 0.3)
	   					continue;

	   			g2.setStroke(stroke);

	   			g2.setColor( flowarrow.color );
	   			AffineTransform at = g2.getTransform();
	   			
	   			// handle the rotation
	   			g2.translate( (int)flowarrow.x , (int)flowarrow.y );
	   			g2.rotate( flowarrow.angle );
	   			g2.translate( (int)-flowarrow.x , (int)-flowarrow.y );

	   			// line
	   			Line2D l = new Line2D.Double(
	   					flowarrow.x - flowarrow.norme/2,
	   					flowarrow.y,
	   					flowarrow.x + flowarrow.norme/2 - stroke.getLineWidth() * 0.25f,
	   					flowarrow.y
	   			);

	   			// arrow head
	   			Path2D.Float path = new Path2D.Float();
	   			path.moveTo(flowarrow.x - flowarrow.norme / 2 + 3*flowarrow.norme / 4,
	   					flowarrow.y + flowarrow.norme / 4);
	   			path.lineTo(flowarrow.x - flowarrow.norme / 2 + flowarrow.norme,
	   					flowarrow.y);
	   			path.lineTo(flowarrow.x - flowarrow.norme / 2 + 3*flowarrow.norme / 4 ,
	   					flowarrow.y - flowarrow.norme / 4);

	   			g2.draw(l);
	   			g2.draw(path);
	   			g2.setTransform(at);
	   		}
		}
   	}
}