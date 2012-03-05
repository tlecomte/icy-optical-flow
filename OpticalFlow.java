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

package plugins.tlecomte.opticalflow;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;

import plugins.adufour.ezplug.*;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import icy.canvas.IcyCanvas;
import icy.gui.dialog.MessageDialog;
import icy.gui.viewer.Viewer;
import icy.image.IcyBufferedImage;
import icy.painter.Painter;

public class OpticalFlow extends EzPlug implements Painter
{
	public EzVarSequence inputSelector = new EzVarSequence("Input");
	public EzVarInteger	channelSelector	= new EzVarInteger("Channel");
	public EzVarDouble alphaSelector = new EzVarDouble("Regularization param.");
	public EzVarInteger	iterSelector = new EzVarInteger("Max. number of iterations");
	Double[] values = {1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10};
	public EzVarDouble epsilonSelector = new EzVarDouble("Iterations termination param.", values, 0 /* default index */, false /*allowUserInput*/);
	public EzVarBoolean hideZeroVelocitiesSelector = new EzVarBoolean("Hide zero velocities", false);
	public EzVarInteger	resolutionSelector = new EzVarInteger("Pixels between neighbour flow arrows");
	public EzButton axisButton;
	public Sequence inputSequence = null;
	ArrayList<ArrayList<FlowArrow>> flowArrowList = new ArrayList<ArrayList<FlowArrow>>();
	
	@Override
	protected void initialize()
	{
		resolutionSelector.setValue(10);
		alphaSelector.setValue(100000.);
		iterSelector.setValue(10000);
		hideZeroVelocitiesSelector.setValue(true);
		
		EzVarListener<Sequence> inputListener = new EzVarListener<Sequence>()
		{
			@Override
			public void variableChanged(EzVar<Sequence> source, Sequence newSequence)
			{
				channelSelector.setValue(0);
				if (newSequence == null)
				{
					channelSelector.setEnabled(false);
				}
				else
				{
					int sizeC = newSequence.getSizeC();
					channelSelector.setMaxValue(sizeC - 1);
					channelSelector.setEnabled(sizeC == 1 ? false : true);
				}
			}
		};

		inputSelector.addVarChangeListener(inputListener);
		
		ActionListener axisListener = new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent event) {
		      	// compute a separate sequence to illustrate the color code
		      	Sequence axisSequence = compute_coloredAxes();
		      	axisSequence.setName("Hue Saturation color code");
		        addSequence(axisSequence);
			}
		};
		
		axisButton = new EzButton("Display flow color code", axisListener);
		
		// sequence selection
		addEzComponent(inputSelector);
		addEzComponent(channelSelector);
		// model parameter
		addEzComponent(alphaSelector);
		// computation
		addEzComponent(epsilonSelector);
		addEzComponent(iterSelector);
		// display
		addEzComponent(resolutionSelector);
		addEzComponent(hideZeroVelocitiesSelector);
		// display, additional
		addEzComponent(axisButton);		
	}
	
	@Override
	protected void execute()
	{
		// main plugin code goes here, and runs in a separate thread
	
		// first clean previous executions (will detach the painters)
		clean();
		
		super.getUI().setProgressBarMessage("Waiting...");
		
        inputSequence = inputSelector.getValue();
        int channel = channelSelector.getValue();
        
        // Check if sequence exists.
        if ( inputSequence == null )
        {
                   MessageDialog.showDialog("Please open a sequence to use this plugin.", MessageDialog.WARNING_MESSAGE );
                   return;
        }

        int z = 0;
    	
    	int numT = inputSequence.getSizeT();
    	
    	// define empty sequences for the velocities maps.
    	Sequence uSequence = new Sequence();
    	Sequence vSequence = new Sequence();
        uSequence.setName("u velocities");
        vSequence.setName("v velocities");
        
        // clear the arrows list
		flowArrowList.clear();
        
    	for (int t = 0; t<numT-1; t++) {
    		// get frames
        	// TODO multi-frame computation
        	
        	int w = inputSequence.getSizeX();
        	int h = inputSequence.getSizeY();
        	
        	// define flow variables
        	double[] u1 = new double[h*w];
        	double[] u2 = new double[h*w];
        	// initialize them
        	Arrays.fill(u1, 0);
        	Arrays.fill(u2, 0);
        	  	    
    		// Get a a direct reference to the data as doubles.
    		double[] I1 = Array1DUtil.arrayToDoubleArray(inputSequence.getDataXY(t  , z, channel), inputSequence.isSignedDataType());
    		double[] I2 = Array1DUtil.arrayToDoubleArray(inputSequence.getDataXY(t+1, z, channel), inputSequence.isSignedDataType());
    	    
        	// Compute optical flow between the two frames using Horn-Schunck method.
        	int maxiter = iterSelector.getValue();
        	double alpha = alphaSelector.getValue();
        	double epsilon = epsilonSelector.getValue();
        	double[] u = HornSchunk.HornSchunkAlgorithm(I1, I2, w, h, alpha, maxiter, epsilon);
        	for (int j = 0; j<h; j++) {
        		for (int k = 0; k<w; k++) {
            		u1[j*w + k] = u[j*(w+1) + k];
            		u2[j*w + k] = u[j*(w+1) + k + (w+1)*(h+1)];	
        		}
        	}
        	
        	// store the results
        	add_velocities_maps_to_sequences(u1, u2, w, h, uSequence, vSequence);
        	
        	update_flow_arrows(u1, u2, w , h);
        	
        	getUI().setProgressBarValue((double) (t) / (double) (numT));
    	}
    	
    	// compute a map of the velocity norm
      	FlowNorm uvNormSequence = new FlowNorm(uSequence, vSequence);
      	
      	// compute a map of the velocity angle
      	FlowAngle uvAngleSequence = new FlowAngle(uSequence, vSequence);
      	
      	// compute a colored map of the velocity norm+angle, coded with hue and saturation
      	FlowMiddlebury uvColoredSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence);
      	
        // Create viewers to watch the velocities sequences.
        //addSequence(uSequence);
        //addSequence(vSequence);
        addSequence(uvNormSequence);
        //addSequence(uvAngleSequence);
        addSequence(uvColoredSequence);
		
        // add a painter to the sequence to draw the arrows
		inputSequence.addPainter( this );
    }
	
	void add_velocities_maps_to_sequences(double[] u, double[] v, int w, int h, Sequence uSequence, Sequence vSequence) {
        // create the image object
        IcyBufferedImage uImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
        IcyBufferedImage vImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
        Object uImageData = uImage.getDataXY(0);
        Object vImageData = vImage.getDataXY(0);
		// Put the velocities data in output images.
		Array1DUtil.doubleArrayToArray( u, uImageData);
		Array1DUtil.doubleArrayToArray( v, vImageData);
    	
        // notify to icy that data has changed to refresh internal state and display
        uImage.dataChanged();
        vImage.dataChanged();

        // add the new images to the sequences at a new time point
        uSequence.setImage(uSequence.getSizeT(), 0 /*z*/, uImage);
        vSequence.setImage(vSequence.getSizeT(), 0 /*z*/, vImage);
	}
	
	Sequence compute_coloredAxes() {
		Sequence uSequence = new Sequence();
		Sequence vSequence = new Sequence();
		
        // create the image object
		int w = 100;
		int h = w;
        IcyBufferedImage uImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
        IcyBufferedImage vImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
        Object uImageData = uImage.getDataXY(0);
        Object vImageData = vImage.getDataXY(0);

        double[] u = new double[w*h];
        double[] v = new double[w*h];
        
        for (int i=0; i<w; i++) {
        	for (int j=0; j<h; j++) {
        		u[j*w + i] = i - w/2;
        		v[j*w + i] = j - h/2;
        	}
        }
        
		Array1DUtil.doubleArrayToArray(u, uImageData);
		Array1DUtil.doubleArrayToArray(v, vImageData);
		
        // notify to icy that data has changed to refresh internal state and display
        uImage.dataChanged();
        vImage.dataChanged();

        // add the new images to the sequences at a new time point
        uSequence.setImage(uSequence.getSizeT(), 0 /*z*/, uImage);
        vSequence.setImage(vSequence.getSizeT(), 0 /*z*/, vImage);
               
		FlowNorm uvNormSequence = new FlowNorm(uSequence, vSequence);
      	FlowAngle uvAngleSequence = new FlowAngle(uSequence, vSequence);
      	FlowMiddlebury axisSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence);
      	
      	// ask Icy core not to distort the colors
      	for (int i=0; i<3; i++) {
      		axisSequence.setComponentUserBoundsAutoUpdate(false);
      		axisSequence.getColorModel().setComponentUserMinValue(i, 0.);
      	}
      	
      	return axisSequence;
	}
	
	void update_flow_arrows(double[] u, double[] v, int w, int h) {
   		/** number of pixel for the square containing the arrow */
   		int resolution = resolutionSelector.getValue();		
    	
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

		// normalize
		double max = 0;
		for ( FlowArrow flowarrow : currentImageArrowList )
		{					
			if ( flowarrow.norme > max ) max = flowarrow.norme ;
		}				

		if ( max > 0 )
		{
			for ( FlowArrow flowarrow : currentImageArrowList )
			{					
				double rapport =  resolution / max ;
				flowarrow.norme *= rapport ;
				flowarrow.vx *= rapport;
				flowarrow.vy *= rapport;

				float norme1 = (float)flowarrow.norme / (float)resolution ;
				// flowarrow.color = new Color( norme1 , 0f , 1f - norme1 ) ;
				flowarrow.color = new Color( norme1 , 1f , 0f ) ;
			}
		}
		
		// add the current arrow map to the global list
		flowArrowList.add( currentImageArrowList );
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
		BasicStroke stroke = new BasicStroke((float) 2.0);

		t = Math.min(t, sequence.getSizeT() - 2); /* replicate flow on last image */
		
		if (!flowArrowList.isEmpty() && t < flowArrowList.size()) {
			ArrayList<FlowArrow> currentImageArrowList = flowArrowList.get(t);
			
	   		for ( FlowArrow flowarrow : currentImageArrowList )
	   		{
	   			if (hideZeroVelocitiesSelector.getValue()) 
	   				if (flowarrow.norme < 0.3)
	   					continue;

	   			g2.setStroke(stroke);

	   			g2.setColor( flowarrow.color );
	   			AffineTransform at = g2.getTransform();
	   			g2.translate( (int)flowarrow.x , (int)flowarrow.y );
	   			g2.rotate( flowarrow.angle );
	   			g2.translate( (int)-flowarrow.x , (int)-flowarrow.y );

	   			Line2D l1 = new Line2D.Double(
	   					flowarrow.x - flowarrow.norme / 2,
	   					flowarrow.y , 
	   					flowarrow.x - flowarrow.norme / 2 + flowarrow.norme ,
	   					flowarrow.y
	   			);

	   			Line2D l2 = new Line2D.Double(
	   					flowarrow.x - flowarrow.norme / 2 + 3*flowarrow.norme / 4,
	   					flowarrow.y + flowarrow.norme / 4,
	   					flowarrow.x - flowarrow.norme / 2 + flowarrow.norme,
	   					flowarrow.y
	   			);

	   			Line2D l3 = new Line2D.Double(
	   					flowarrow.x - flowarrow.norme / 2 + 3*flowarrow.norme / 4 ,
	   					flowarrow.y - flowarrow.norme / 4 ,
	   					flowarrow.x - flowarrow.norme / 2 + flowarrow.norme ,
	   					flowarrow.y 
	   			);

	   			g2.draw( l1 );
	   			g2.draw( l2 );
	   			g2.draw( l3 );
	   			g2.setTransform( at );
	   		}
			
		}
   	}

   	public double getDistance( double x1 , double y1 , double z1 , double x2 , double y2 , double z2 )
   	{
   		double distance = Math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );		
   		return distance;
   	}

   	@Override
	public void clean()
	{
		// use this method to clean local variables or input streams (if any) to avoid memory leaks
   		if (inputSequence != null) {
   			inputSequence.removePainter( this );
   			inputSequence = null;	
   		}
	}
}

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