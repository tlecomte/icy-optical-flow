package plugins.tlecomte.middleburyColorCoder;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import icy.gui.dialog.MessageDialog;
import icy.sequence.Sequence;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.vars.lang.VarSequence;
import plugins.tlecomte.flowdisplay.FlowAngle;
import plugins.tlecomte.flowdisplay.FlowMiddlebury;
import plugins.tlecomte.flowdisplay.FlowNorm;

public class MiddleburyColorCoder extends EzPlug implements Block {
	public EzVarSequence uxSequenceSelector = new EzVarSequence("ux Sequence");
	public EzVarSequence uySequenceSelector = new EzVarSequence("uy Sequence");
	public EzButton axisButton;
	
	Sequence uxSequence = null;
	Sequence uySequence = null;
	
	VarSequence colorSequenceVar = new VarSequence("Color-coded flow", null);
	
	protected void initialize() {
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
		addEzComponent(uxSequenceSelector);
		uxSequenceSelector.setToolTipText("<html>Choose a sequence for the u_x flow.</html>");
		addEzComponent(uySequenceSelector);
		uySequenceSelector.setToolTipText("<html>Choose a sequence for the u_y flow.</html>");
		
		addEzComponent(axisButton);	
		axisButton.setToolTipText(    "Click here to display a reference image of the color code used"
									+ " to display the 2D flow.");
	}
	
	// declare ourself to Blocks
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(uxSequenceSelector.getVariable());
		inputMap.add(uySequenceSelector.getVariable());
	}

	// declare ourself to Blocks
	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add(colorSequenceVar);
	}
	
	protected void execute() {
		uxSequence = uxSequenceSelector.getValue();
		uySequence = uySequenceSelector.getValue();
		
        // Check if sequence exists.
        if ( (uxSequence == null)
        		|| (uySequence == null))
        {
        	if (getUI() != null) {
        		MessageDialog.showDialog("Please choose two displacement sequences to use this plugin.", MessageDialog.ERROR_MESSAGE );
    			return;	
        	} else {
        		// FIXME figure out what to do when headless or in BLocks
        		return;
        	}
        }
		
    	// compute a map of the velocity norm
      	FlowNorm uvNormSequence = new FlowNorm(uxSequence, uySequence, "");
      	
      	// compute a map of the velocity angle
      	FlowAngle uvAngleSequence = new FlowAngle(uxSequence, uySequence, "");
      	
      	// compute a colored map of the velocity norm+angle, coded with hue and saturation
      	FlowMiddlebury colorSequence = new FlowMiddlebury(uvNormSequence, uvAngleSequence, "");
        
    	// assign the vars that are used by the Blocks interface
    	colorSequenceVar.setValue(colorSequence);
    	
      	if (getUI() != null) {
            addSequence(colorSequence);
      	}
	}

	public void clean() {
	}
}
