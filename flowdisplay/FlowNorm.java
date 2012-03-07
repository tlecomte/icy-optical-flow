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

package plugins.tlecomte.flowdisplay;

import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;

public class FlowNorm extends Sequence {
  	public FlowNorm(Sequence uSequence, Sequence vSequence) {
  		int numT = uSequence.getSizeT();
  		
      	for (int i = 0; i<numT; i++) {
    		// get frames
        	Object uImageData = uSequence.getDataXY(i, 0 /*z*/, 0 /* channel */);
        	Object vImageData = vSequence.getDataXY(i, 0 /*z*/, 0 /* channel */);
        	
    		// get frames data
        	double[] uDataBuffer = Array1DUtil.arrayToDoubleArray(uImageData, uSequence.isSignedDataType());
        	double[] vDataBuffer = Array1DUtil.arrayToDoubleArray(vImageData, vSequence.isSignedDataType());
        	
        	int w = uSequence.getSizeX();
        	int h = uSequence.getSizeY();
        	
        	double[] uvNormDataBuffer = new double[w*h];
        	
        	for (int j = 0; j<w*h; j++) {
        		uvNormDataBuffer[j] = Math.sqrt(uDataBuffer[j]*uDataBuffer[j] + vDataBuffer[j]*vDataBuffer[j]);
        	}
        	
            // create the image object
            IcyBufferedImage uvNormImage = new IcyBufferedImage(w, h, 1, DataType.getDataType("double"));
            Object uvNormImageData = uvNormImage.getDataXY(0);
    		// Put the velocities data in output images.
    		Array1DUtil.doubleArrayToArray(uvNormDataBuffer, uvNormImageData);
        	
            // notify to icy that data has changed to refresh internal state and display
            uvNormImage.dataChanged();

            // add the new images to the sequences at a new time point
            setImage(getSizeT(), 0 /*z*/, uvNormImage);
            
            setName("Flow norm");
      	}      		
  	}	
}
