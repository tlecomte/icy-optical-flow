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

public class FlowMiddlebury extends Sequence {
  	public FlowMiddlebury(Sequence uvNormSequence, Sequence uvAngleSequence, String namePrefix) {
		int numT = uvNormSequence.getSizeT();
  		
      	for (int i = 0; i<numT; i++) {
    		// get frames
        	Object uvNormImageData = uvNormSequence.getDataXY(i, 0 /*z*/, 0 /* channel */);
        	Object uvAngleImageData = uvAngleSequence.getDataXY(i, 0 /*z*/, 0 /* channel */);
        	
    		// get frames data
        	double[] uvNormDataBuffer = Array1DUtil.arrayToDoubleArray(uvNormImageData, uvNormSequence.isSignedDataType());
        	double[] uvAngleDataBuffer = Array1DUtil.arrayToDoubleArray(uvAngleImageData, uvAngleSequence.isSignedDataType());
        	
        	int w = uvNormSequence.getSizeX();
        	int h = uvNormSequence.getSizeY();
        	
        	double[] uvRDataBuffer = new double[w*h];
        	double[] uvGDataBuffer = new double[w*h];
        	double[] uvBDataBuffer = new double[w*h];
        	
        	// find max norm
        	double maxNorm = 0;
        	for (int j = 0; j<w*h; j++) {
        		if (maxNorm < uvNormDataBuffer[j]) maxNorm = uvNormDataBuffer[j];
        	}
        	
        	for (int j = 0; j<w*h; j++) {
        		// define HSV colors from velocities
        		double Hue = uvAngleDataBuffer[j]; // Hue in radians, from 0 to 2*pi
        		double Saturation = uvNormDataBuffer[j]/maxNorm; // from 0 to 1
        		double Value = 1.; // from 0 to 1
        		
        		// convert to RGB
        		// chroma
        		double C = Value*Saturation;
        		
        		double H_prime = Hue / (Math.PI/3.);
        		double X = C * (1. - Math.abs((H_prime % 2.) - 1.));
        		double m = Value - C;
        		
    			if (Hue < Math.PI/3.) {
            		uvRDataBuffer[j] = C + m;
            		uvGDataBuffer[j] = X + m;
            		uvBDataBuffer[j] = m;
    			} else if (Hue < 2*Math.PI/3.) {
            		uvRDataBuffer[j] = X + m;
            		uvGDataBuffer[j] = C + m;
            		uvBDataBuffer[j] = m;
    			} else if (Hue < Math.PI) {
            		uvRDataBuffer[j] = m;
            		uvGDataBuffer[j] = C + m;
            		uvBDataBuffer[j] = X + m;
    			} else if (Hue < 4*Math.PI/3.) { 
            		uvRDataBuffer[j] = m;
            		uvGDataBuffer[j] = X + m;
            		uvBDataBuffer[j] = C + m;
    			} else if (Hue < 5*Math.PI/3.) {
            		uvRDataBuffer[j] = X + m;
            		uvGDataBuffer[j] = m;
            		uvBDataBuffer[j] = C + m;
    			} else { // Hue < 2*Math.PI
            		uvRDataBuffer[j] = C + m;
            		uvGDataBuffer[j] = m;
            		uvBDataBuffer[j] = X + m;
    			}
        	}
        	
            // create the image object
            IcyBufferedImage uvRGBImage = new IcyBufferedImage(w, h, 3 /* color! */, DataType.getDataType("double"));
            Object uvRImageData = uvRGBImage.getDataXY(0);
            Object uvGImageData = uvRGBImage.getDataXY(1);
            Object uvBImageData = uvRGBImage.getDataXY(2);
    		// Put the velocities data in output images.
    		Array1DUtil.doubleArrayToArray(uvRDataBuffer, uvRImageData);
    		Array1DUtil.doubleArrayToArray(uvGDataBuffer, uvGImageData);
    		Array1DUtil.doubleArrayToArray(uvBDataBuffer, uvBImageData);
        	
            // notify to icy that data has changed to refresh internal state and display
            uvRGBImage.dataChanged();

            // add the new images to the sequences at a new time point
            setImage(getSizeT(), 0 /*z*/, uvRGBImage);
      	}
      	
      	// ask Icy core not to distort the colors
      	for (int i=0; i<3; i++) {
      		setComponentUserBoundsAutoUpdate(false);
      		getColorModel().setComponentUserMinValue(i, 0.);
      	}
      	
      	setName(namePrefix + " Color-coded flow");
	}
}
