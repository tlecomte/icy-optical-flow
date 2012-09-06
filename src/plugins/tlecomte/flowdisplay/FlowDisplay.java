package plugins.tlecomte.flowdisplay;

import java.util.List;

import icy.gui.dialog.ConfirmDialog;
import icy.gui.dialog.MessageDialog;
import icy.painter.Painter;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;

public class FlowDisplay extends EzPlug implements Block {
	public EzVarSequence coverSequenceSelector = new EzVarSequence("Cover Sequence");
	public EzVarSequence uxSequenceSelector = new EzVarSequence("ux Sequence");
	public EzVarSequence uySequenceSelector = new EzVarSequence("uy Sequence");
	public EzVarInteger	resolutionSelector = new EzVarInteger("Pixels between neighbour flow arrows", 10, 1, Integer.MAX_VALUE, 1);
	public EzVarBoolean hideZeroVelocitiesSelector = new EzVarBoolean("Hide zero velocities", true);
	public EzVarBoolean removePreviousPaintersSelector = new EzVarBoolean("Remove previous flow painters", true);
	
	Sequence uxSequence = null;
	Sequence uySequence = null;
	Sequence coverSequence = null;
	
	protected void initialize() {
		// sequence selection
		addEzComponent(coverSequenceSelector);
		coverSequenceSelector.setToolTipText("<html>Choose a sequence where the flow will be displayed onto.</html>");
		addEzComponent(uxSequenceSelector);
		uxSequenceSelector.setToolTipText("<html>Choose a sequence for the u_x flow.</html>");
		addEzComponent(uySequenceSelector);
		uySequenceSelector.setToolTipText("<html>Choose a sequence for the u_y flow.</html>");
		
		// display
		addEzComponent(resolutionSelector);
		addEzComponent(hideZeroVelocitiesSelector);
		hideZeroVelocitiesSelector.setToolTipText(  "<html>If checked, the very smaller flow vectors will not be<br>"
												  + "displayed on top of the sequence, so that the visualization is clearer.</html>");
		addEzComponent(removePreviousPaintersSelector);
		removePreviousPaintersSelector.setToolTipText(  "<html>If checked, the previous flow painters that are present on the sequence<br>"
													  + "will be removed. Uncheck if you want to preserve them.</html>");

		// display, additional
		//addEzComponent(axisButton);
		//axisButton.setToolTipText(    "Click here to display a reference image of the color code used"
		//							+ " to display the 2D flow.");	
	}
	
	// declare ourself to Blocks
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(coverSequenceSelector.getVariable());
		inputMap.add(uxSequenceSelector.getVariable());
		inputMap.add(uySequenceSelector.getVariable());
		inputMap.add(resolutionSelector.getVariable());
		inputMap.add(hideZeroVelocitiesSelector.getVariable());
		inputMap.add(removePreviousPaintersSelector.getVariable());
	}

	// declare ourself to Blocks
	@Override
	public void declareOutput(VarList outputMap) {
		// we have no output, the block adds a painter to the cover sequence
	}
	
	protected void execute() {
		VectorFlowPainter flowPainter = new VectorFlowPainter();

		coverSequence = coverSequenceSelector.getValue();
		uxSequence = uxSequenceSelector.getValue();
		uySequence = uySequenceSelector.getValue();
		
        // Check if sequence exists.
        if ( (uxSequence == null)
        		|| (uySequence == null)
        		|| (coverSequence == null))
        {
    		MessageDialog.showDialog("Please open three sequences to use this plugin.", MessageDialog.ERROR_MESSAGE );
    		return;
        }
		
		int w = uxSequence.getSizeX();
    	int h = uxSequence.getSizeY();
		int numT = uxSequence.getSizeT();
		int coverNumT = coverSequence.getSizeT();
		int z = 0;
		int channel = 0;
		
		if ((w != uySequence.getSizeX())
				|| (h != uySequence.getSizeY())
				|| (numT != uySequence.getSizeT())
				|| ((numT != coverNumT) && (numT + 1 != coverNumT))) {
			MessageDialog.showDialog("Sequences sizes are not compatible.", MessageDialog.ERROR_MESSAGE );
			return;
		}
		
		if ((w != coverSequence.getSizeX())
				|| (h != coverSequence.getSizeY())) {
			boolean ret = ConfirmDialog.confirm("Sizes of flow and cover sequences are different, are you sure you want to continue ?");
			if (!ret) {
				return;
			}
		}
		
		flowPainter.clear();
		flowPainter.hideZeroVelocities(hideZeroVelocitiesSelector.getValue());
		flowPainter.setResolution(resolutionSelector.getValue());
		
		for (int t = 0; t<numT; t++) {
			// Get a a direct reference to the data as doubles.
			double[] ux = Array1DUtil.arrayToDoubleArray(uxSequence.getDataXY(t, z, channel), uxSequence.isSignedDataType());
			double[] uy = Array1DUtil.arrayToDoubleArray(uySequence.getDataXY(t, z, channel), uySequence.isSignedDataType());
			
			flowPainter.update_flow_arrows(ux, uy, w, h);
		}
		
		// usually the flow sequences have one less time step than the input sequence
		// if it is the case, we replicate the last flow image at the end
		if (numT == coverNumT - 1) {
			int t = numT - 1;
			// Get a a direct reference to the data as doubles.
			double[] ux = Array1DUtil.arrayToDoubleArray(uxSequence.getDataXY(t, z, channel), uxSequence.isSignedDataType());
			double[] uy = Array1DUtil.arrayToDoubleArray(uySequence.getDataXY(t, z, channel), uySequence.isSignedDataType());
			
			flowPainter.update_flow_arrows(ux, uy, w, h);			
		}
		
		flowPainter.normalize();
		
		if (removePreviousPaintersSelector.getValue()) {
	    	// remove previous painters
			List<Painter> painters = coverSequence.getPainters(VectorFlowPainter.class);
			for (Painter painter : painters) {
				coverSequence.removePainter(painter);
				coverSequence.painterChanged(painter);
			}			
		}
			
        // add a painter to the sequence to draw the arrows
		coverSequence.addPainter(flowPainter);
	}

	public void clean() {
	}
}
