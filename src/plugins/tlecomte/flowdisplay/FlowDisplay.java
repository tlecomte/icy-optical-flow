package plugins.tlecomte.flowdisplay;

import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;

public class FlowDisplay extends EzPlug {
	public EzVarSequence uxSequenceSelector = new EzVarSequence("ux Sequence");
	public EzVarSequence uySequenceSelector = new EzVarSequence("uy Sequence");
	public EzVarSequence coverSequenceSelector = new EzVarSequence("Cover Sequence");
	public EzVarInteger	resolutionSelector = new EzVarInteger("Pixels between neighbour flow arrows");
	public EzVarBoolean hideZeroVelocitiesSelector = new EzVarBoolean("Hide zero velocities", false);
	
	Sequence uxSequence = null;
	Sequence uySequence = null;
	Sequence coverSequence = null;
	
	protected void initialize() {
		resolutionSelector.setValue(10);
		hideZeroVelocitiesSelector.setValue(true);
		
		// sequence selection
		addEzComponent(uxSequenceSelector);
		uxSequenceSelector.setToolTipText("<html>Choose a sequence for the u_x flow.</html>");
		addEzComponent(uySequenceSelector);
		uySequenceSelector.setToolTipText("<html>Choose a sequence for the u_y flow.</html>");
		addEzComponent(coverSequenceSelector);
		coverSequenceSelector.setToolTipText("<html>Choose a sequence where the flow will be displayed onto.</html>");
		
		// display
		addEzComponent(resolutionSelector);
		addEzComponent(hideZeroVelocitiesSelector);
		hideZeroVelocitiesSelector.setToolTipText(  "<html>If checked, the very smaller flow vectors will not be<br>"
												  + "displayed on top of the sequence, so that the visualization is clearer.</html>");
		
		// display, additional
		//addEzComponent(axisButton);
		//axisButton.setToolTipText(    "Click here to display a reference image of the color code used"
		//							+ " to display the 2D flow.");	
	}
	
	protected void execute() {
		VectorFlowPainter flowPainter = new VectorFlowPainter();
		
		uxSequence = uxSequenceSelector.getValue();
		uySequence = uySequenceSelector.getValue();
		coverSequence = coverSequenceSelector.getValue();
		
		int w = uxSequence.getSizeX();
    	int h = uxSequence.getSizeY();
		int numT = uxSequence.getSizeT();
		int z = 0;
		int channel = 0;
		
		flowPainter.clear();
		flowPainter.hideZeroVelocities(hideZeroVelocitiesSelector.getValue());
		flowPainter.setResolution(resolutionSelector.getValue());
		
		for (int t = 0; t<numT-1; t++) {
			// Get a a direct reference to the data as doubles.
			double[] ux = Array1DUtil.arrayToDoubleArray(uxSequence.getDataXY(t, z, channel), uxSequence.isSignedDataType());
			double[] uy = Array1DUtil.arrayToDoubleArray(uySequence.getDataXY(t, z, channel), uySequence.isSignedDataType());
			
			flowPainter.update_flow_arrows(ux, uy, w, h);
		}
		
		flowPainter.normalize();
			
        // add a painter to the sequence to draw the arrows
		coverSequence.addPainter(flowPainter);
	}

	public void clean() {
	}
}
