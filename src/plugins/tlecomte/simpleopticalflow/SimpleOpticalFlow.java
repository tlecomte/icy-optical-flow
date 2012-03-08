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

package plugins.tlecomte.simpleopticalflow;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import plugins.adufour.ezplug.*;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import plugins.tlecomte.flowdisplay.FlowAngle;
import plugins.tlecomte.flowdisplay.FlowMiddlebury;
import plugins.tlecomte.flowdisplay.FlowNorm;

public class SimpleOpticalFlow extends EzPlug
{
	public EzGroup inputGroup = new EzGroup("Input");
	public EzVarSequence sequenceSelector = new EzVarSequence("Sequence");
	public EzVarInteger	channelSelector	= new EzVarInteger("Channel");
	
	public EzGroup modelGroup = new EzGroup("Flow parameters");
	public EzVarDouble alphaSelector = new EzVarDouble("Regularization parameter");
	
	public EzGroup computationGroup = new EzGroup("Computation parameters");
	Double[] values = {1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10};
	public EzVarDouble epsilonSelector = new EzVarDouble("Tolerance before termination", values, 2 /* default index */, false /*allowUserInput*/);
	public EzVarInteger	iterSelector = new EzVarInteger("Maximum number of iterations");
	
	public EzGroup displayGroup = new EzGroup("Display options for the vector flow overlay");
	public EzVarBoolean hideZeroVelocitiesSelector = new EzVarBoolean("Hide zero velocities", false);
	public EzVarInteger	resolutionSelector = new EzVarInteger("Pixels between neighbour flow arrows");
	
	public EzButton axisButton;
	
	public Sequence inputSequence = null;
	
	VectorFlowPainter flowPainter = new VectorFlowPainter();
	
	@Override
	protected void initialize()
	{
		resolutionSelector.setValue(10);
		alphaSelector.setValue(100000.);
		iterSelector.setValue(10000);
		hideZeroVelocitiesSelector.setValue(true);
		
		EzVarListener<Sequence> sequenceListener = new EzVarListener<Sequence>()
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

		sequenceSelector.addVarChangeListener(sequenceListener);
		
		ActionListener axisListener = new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent event) {
		      	// compute a separate sequence to illustrate the color code
		      	Sequence axisSequence = compute_coloredAxes();
		        addSequence(axisSequence);
			}
		};
		
		axisButton = new EzButton("Display flow color code", axisListener);
		
		// sequence selection
		addEzComponent(sequenceSelector);
		sequenceSelector.setToolTipText(  "<html>Choose a sequence. The optical flow will be computed from consecutive<br>"
										+ "frames of that sequence, in the z=0 plane.</html>");
		addEzComponent(channelSelector);
		channelSelector.setToolTipText("The optical flow will be computed from the data in this channel only.");
		inputGroup.addEzComponent(sequenceSelector, channelSelector);
		addEzComponent(inputGroup);
		
		// model parameter
		addEzComponent(alphaSelector);
		alphaSelector.setToolTipText(  "<html>Choose the value of the regularisation parameter. Choose a larger parameter<br>"
									 + "to make the flow smoother, so that it will be more robust against noise.<br>"
									 + "Instead, choose a smaller parameter if your flow is too uniform.</html>");
		modelGroup.addEzComponent(alphaSelector);
		addEzComponent(modelGroup);
		
		// computation
		addEzComponent(epsilonSelector);
		epsilonSelector.setToolTipText(  "<html>Choose the tolerance to achieve in the flow computation. The algorithm<br>"
									   + "will terminate when the relative residue is smaller than this value.<br>"
									   + "A tolerance of 1E-4 is a good start. Increase the tolerance to make<br>"
									   + "the computation stop earlier, decrease it if the results is not realistic.</html>");
		addEzComponent(iterSelector);
		iterSelector.setToolTipText(  "<html>Choose the maximum number of iterations of the flow computation. The<br>" 
									+ "computation will stop after this number of steps, even if the tolerance (as<br>"
									+ "specified in the previous parameter) is not reached.</html>");
		computationGroup.addEzComponent(epsilonSelector, iterSelector);
		addEzComponent(computationGroup);
		
		// display
		addEzComponent(resolutionSelector);
		addEzComponent(hideZeroVelocitiesSelector);
		hideZeroVelocitiesSelector.setToolTipText(  "<html>If checked, the very smaller flow vectors will not be<br>"
												  + "displayed on top of the sequence, so that the visualization is clearer.</html>");
		displayGroup.addEzComponent(resolutionSelector, hideZeroVelocitiesSelector);
		addEzComponent(displayGroup);
		
		// display, additional
		addEzComponent(axisButton);	
		//axisButton.setToolTipText("Click here to display a reference image of the color code used" +
		//	" to display the 2D flow.");
	}
	
	@Override
	protected void execute()
	{
		// main plugin code goes here, and runs in a separate thread
	
		// first clean previous executions (will detach the painters)
		clean();
		
		super.getUI().setProgressBarMessage("Waiting...");
		
        inputSequence = sequenceSelector.getValue();
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
		flowPainter.clear();
        
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
        	
        	flowPainter.hideZeroVelocities(hideZeroVelocitiesSelector.getValue());
        	flowPainter.update_flow_arrows(u1, u2, w , h, resolutionSelector.getValue());
        	
        	getUI().setProgressBarValue((double) (t) / (double) (numT));
    	}
    	
    	// compute a map of the velocity norm
      	FlowNorm uvNormSequence = new FlowNorm(uSequence, vSequence, inputSequence.getName());
      	
      	// compute a map of the velocity angle
      	FlowAngle uvAngleSequence = new FlowAngle(uSequence, vSequence, inputSequence.getName());
      	
      	// compute a colored map of the velocity norm+angle, coded with hue and saturation
      	FlowMiddlebury uvColoredSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence, inputSequence.getName());
      	
        // Create viewers to watch the velocities sequences.
        //addSequence(uSequence);
        //addSequence(vSequence);
        addSequence(uvNormSequence);
        //addSequence(uvAngleSequence);
        addSequence(uvColoredSequence);
		
        // add a painter to the sequence to draw the arrows
		inputSequence.addPainter(flowPainter);
    }
	
	void add_velocities_maps_to_sequences(double[] u, double[] v, int w, int h, Sequence uSequence, Sequence vSequence) {
        // create the image object
        IcyBufferedImage uImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
        IcyBufferedImage vImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));

        // Put the velocities data in output images.
		Array1DUtil.doubleArrayToArray( u, uImage.getDataXY(0));
		Array1DUtil.doubleArrayToArray( v, vImage.getDataXY(0));
    	
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
               
		FlowNorm uvNormSequence = new FlowNorm(uSequence, vSequence, "Reference");
      	FlowAngle uvAngleSequence = new FlowAngle(uSequence, vSequence, "Reference");
      	FlowMiddlebury axisSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence, "Reference");
      	
      	// ask Icy core not to distort the colors
      	for (int i=0; i<3; i++) {
      		axisSequence.setComponentUserBoundsAutoUpdate(false);
      		axisSequence.getColorModel().setComponentUserMinValue(i, 0.);
      	}
      	
      	return axisSequence;
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
   			inputSequence.removePainter(flowPainter);
   			inputSequence = null;	
   		}
	}
}
