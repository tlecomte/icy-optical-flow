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

package plugins.tlecomte.opticalFlowHornSchunck;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.List;

import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.*;
import plugins.adufour.vars.lang.VarSequence;
import icy.painter.Painter;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.image.colormap.FireColorMap;
import icy.image.colormap.JETColorMap;
import plugins.tlecomte.flowdisplay.FlowAngle;
import plugins.tlecomte.flowdisplay.FlowNorm;
import plugins.tlecomte.flowdisplay.VectorFlowPainter;
import plugins.tlecomte.middleburyColorCoder.FlowMiddlebury;

public class OpticalFlowHornSchunck extends EzPlug implements Block
{
	public EzGroup inputGroup = new EzGroup("Input");
	public EzVarSequence sequenceSelector = new EzVarSequence("Sequence");
	public EzVarInteger	channelSelector	= new EzVarInteger("Channel");
	
	public EzGroup modelGroup = new EzGroup("Flow parameters");
	public EzVarDouble alphaSelector = new EzVarDouble("Regularization parameter", 100000., 0., Double.MAX_VALUE, 10000.);
	
	public EzGroup computationGroup = new EzGroup("Computation parameters");
	Double[] values = {1E-2, 1E-3, 1E-4, 1E-5, 1E-6, 1E-7, 1E-8, 1E-9, 1E-10};
	public EzVarDouble epsilonSelector = new EzVarDouble("Tolerance before termination", values, 2 /* default index */, false /*allowUserInput*/);
	public EzVarInteger	iterSelector = new EzVarInteger("Maximum number of iterations", 10000, 1, Integer.MAX_VALUE, 1000);
	
	public EzGroup displayGroup = new EzGroup("Display options for the vector flow overlay");
	public EzVarBoolean hideZeroVelocitiesSelector = new EzVarBoolean("Hide zero velocities", true);
	public EzVarInteger	resolutionSelector = new EzVarInteger("Pixels between neighbour flow arrows", 10, 1, Integer.MAX_VALUE, 1);
	
	public EzGroup outputGroup = new EzGroup("Output options");
	public EzVarBoolean flowMapSelector = new EzVarBoolean("Horizontal and vertical flows", false);
	public EzVarBoolean flowNormSelector = new EzVarBoolean("Flow norm", true);
	public EzVarBoolean colorFlowSelector = new EzVarBoolean("Color-coded flow", true);
	public EzButton axisButton;
	
	VarSequence uSequenceVar = new VarSequence("Horizontal flow", null);
	VarSequence vSequenceVar = new VarSequence("Vertical flow", null);
	
	@Override
	protected void initialize()
	{
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
		      	Sequence axisSequence = FlowMiddlebury.coloredAxes();
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
		
		// output
		addEzComponent(flowMapSelector);
		flowMapSelector.setToolTipText(  "<html>Will output two sequences with the horizontal and<br>"
									   + "vertical displacements, respectively.</html>");
		addEzComponent(flowNormSelector);
		flowNormSelector.setToolTipText(  "<html>Will output a sequence with the norm of the flow.</html>");
		addEzComponent(colorFlowSelector);
		colorFlowSelector.setToolTipText(  "<html>Will output a sequence where the flow is displayed<br>"
				   						 + "with the Middlebury color-code.</html>");
		addEzComponent(axisButton);	
		axisButton.setToolTipText(    "Click here to display a reference image of the color code used"
									+ " to display the 2D flow.");
		outputGroup.addEzComponent(flowMapSelector, flowNormSelector, colorFlowSelector, axisButton);

		addEzComponent(outputGroup);
	}
	
	// declare ourself to Blocks
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(sequenceSelector.getVariable());
		inputMap.add(channelSelector.getVariable());
		inputMap.add(alphaSelector.getVariable());
		inputMap.add(epsilonSelector.getVariable());
		inputMap.add(iterSelector.getVariable());
	}

	// declare ourself to Blocks
	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add(uSequenceVar);
		outputMap.add(vSequenceVar);
	}
	
	@Override
	protected void execute()
	{
		// main plugin code goes here, and runs in a separate thread
	
		// first clean previous executions (will detach the painters)
		clean();
		
		if (getUI() != null) {
			getUI().setProgressBarMessage("Waiting...");
		}
		
        Sequence inputSequence = sequenceSelector.getValue();
        int channel = channelSelector.getValue();
        
        // Check if sequence exists.
        if ( inputSequence == null )
        {
        	if (getUI() != null) {
        		MessageDialog.showDialog("Please open a sequence to use this plugin.", MessageDialog.ERROR_MESSAGE );
    			return;	
        	} else {
        		// FIXME figure out what to do when headless or in BLocks
        		return;
        	}
        }

        int z = 0;
    	
    	int numT = inputSequence.getSizeT();
    	
    	if ( numT < 2 ) {
    		if (getUI() != null) {
    			MessageDialog.showDialog("The input sequence should have at least two successive images\n"
    					+ "for the optical flow computation.\n\n"
    					+ "Note: If you want to compute the optical flow on a z-stack,\n"
    					+ "first convert it to a time sequence\n"
    					+ "(\"Sequence / Image operation\" -> \"Convert to time\").", MessageDialog.ERROR_MESSAGE );
    			return;
    		} else {
        		// FIXME figure out what to do when headless or in BLocks
        		return;
        	}
    	}
    	
    	// define empty sequences for the velocities maps.
    	Sequence uSequence = new Sequence();
    	Sequence vSequence = new Sequence();
        uSequence.setName("Horizontal flow");
        vSequence.setName("Vertical flow");
        
        VectorFlowPainter flowPainter = new VectorFlowPainter();
        if (getUI() != null) {
        	// remove previous painters
    		List<Painter> painters = inputSequence.getPainters(VectorFlowPainter.class);
    		for (Painter painter : painters) {
    			inputSequence.removePainter(painter);
        		inputSequence.painterChanged(painter);
    		}
    		
	    	flowPainter.hideZeroVelocities(hideZeroVelocitiesSelector.getValue());
	    	flowPainter.setResolution(resolutionSelector.getValue());
        }
        
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
        	     	
        	if (getUI() != null) {
            	flowPainter.update_flow_arrows(u1, u2, w , h);
        		getUI().setProgressBarValue((double) (t) / (double) (numT));
        	}
    	}
    	
		// set a JET colormap and make the bounds symmetric
		// so that the zero corresponds to the green at the center of the colormap
		double[] bounds = uSequence.getChannelBounds(0);
		double ext = Math.max(Math.abs(bounds[0]), Math.abs(bounds[1]));
		uSequence.setAutoUpdateChannelBounds(false);
		uSequence.getColorModel().setComponentUserBounds(0, -ext, ext);
		
		bounds = vSequence.getChannelBounds(0);
		ext = Math.max(Math.abs(bounds[0]), Math.abs(bounds[1]));
		vSequence.setAutoUpdateChannelBounds(false);
		vSequence.getColorModel().setComponentUserBounds(0, -ext, ext);
		
		uSequence.getColorModel().setColormap(0, new JETColorMap());
		vSequence.getColorModel().setColormap(0, new JETColorMap());
    	
    	
    	// assign the vars that are used by the Blocks interface
    	uSequenceVar.setValue(uSequence);
    	vSequenceVar.setValue(vSequence);
    	
    	// here goes the non-Block outputs
      	if (getUI() != null) {
	    	// compute a map of the velocity norm
	      	FlowNorm uvNormSequence = new FlowNorm(uSequence, vSequence, inputSequence.getName());
	      	      	     	
            // Create viewers to watch the velocities sequences.
      		if (flowMapSelector.getValue()) {     			
      			addSequence(uSequence);
                addSequence(vSequence);
      		}
            if (flowNormSelector.getValue()) {
    	      	uvNormSequence.getColorModel().setColormap(0, new FireColorMap());
            	addSequence(uvNormSequence);	
            }
            if (colorFlowSelector.getValue()) {
    	      	// compute a map of the velocity angle
    	      	FlowAngle uvAngleSequence = new FlowAngle(uSequence, vSequence, inputSequence.getName());

            	// compute a colored map of the velocity norm+angle, coded with hue and saturation
            	FlowMiddlebury uvColoredSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence, inputSequence.getName());
            	addSequence(uvColoredSequence);
            }

            // add a painter to the sequence to draw the arrows
          	flowPainter.normalize();
    		inputSequence.addPainter(flowPainter);
      	}
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

   	public double getDistance( double x1 , double y1 , double z1 , double x2 , double y2 , double z2 )
   	{
   		double distance = Math.sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );		
   		return distance;
   	}

   	@Override
	public void clean()
	{
		// use this method to clean local variables or input streams (if any) to avoid memory leaks
	}
}
