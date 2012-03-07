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
package plugins.tlecomte.simpleopticalflow;

public class fem_functions {
	/* 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double Nbar(double xi, double x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;
	    return (x - xi2)*(x - xi3)/((xi - xi2)*(xi - xi3));
	}
	
	static double[] Nbar_1D(double xi, double[] x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;

        double[] out = new double[x.length]; 
        for (int i = 0; i < x.length; i++)
            out[i] = (x[i] - xi2)*(x[i] - xi3)/((xi - xi2)*(xi - xi3));
	    
        return out;
	}
	
//	static double[][] Nbar_2D(double xi, double[][] x) {
//	    assert (xi == -1. || xi == 0. || xi == 1.);   
//	    double xi2 = (xi + 1 + 1)%3 - 1;
//	    double xi3 = (xi + 2 + 1)%3 - 1;
//	    return (x - xi2)*(x - xi3)/((xi - xi2)*(xi - xi3));
//	}
	
	/* derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double Nbar_prime(double xi, double x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;
	    return (2*x - xi2 - xi3)/((xi - xi2)*(xi - xi3));
	}
	
	/* derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double[] Nbar_prime_1D(double xi, double[] x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;
	    
        double[] out = new double[x.length]; 
        for (int i = 0; i < x.length; i++)
            out[i] = (2*x[i] - xi2 - xi3)/((xi - xi2)*(xi - xi3));
	    
        return out;
	}
	
	/* derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
//	static double[][] Nbar_prime_2D(double xi, double[][] x) {
//	    assert (xi == -1. || xi == 0. || xi == 1.);   
//	    double xi2 = (xi + 1 + 1)%3 - 1;
//	    double xi3 = (xi + 2 + 1)%3 - 1;
//	    return (2*x - xi2 - xi3)/((xi - xi2)*(xi - xi3));
//	}
	
	/* second derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double Nbar_prime_prime(double xi, double x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;
	    return 2./((xi - xi2)*(xi - xi3));
	}
	
	/* second derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double[] Nbar_prime_prime_1D(double xi, double[] x) {
	    assert (xi == -1. || xi == 0. || xi == 1.);   
	    double xi2 = (xi + 1 + 1)%3 - 1;
	    double xi3 = (xi + 2 + 1)%3 - 1;
	    
        double[] out = new double[x.length]; 
        for (int i = 0; i < x.length; i++)
            out[i] = 2./((xi - xi2)*(xi - xi3));
	    
        return out;
	}
	   
	/* second derivative of the 1D quadratic shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
//	static double[][] Nbar_prime_prime_2D(double xi, double[][] x) {
//	    assert (xi == -1. || xi == 0. || xi == 1.);   
//	    double xi2 = (xi + 1 + 1)%3 - 1;
//	    double xi3 = (xi + 2 + 1)%3 - 1;
//	    return 2./((xi - xi2)*(xi - xi3));
//	}
	
	/* 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double Mbar(double xi, double x) {
		assert(xi == 1. || xi == -1.);
		double xi2 = 0.;
	    if (xi == 1.)
	        xi2 = -1.;
	    else if (xi == -1.)
	        xi2 = 1.;
	    return (x - xi2)/(xi - xi2);
	}
	
	/* 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double[] Mbar_1D(double xi, double[] x) {
		assert(xi == 1. || xi == -1.);
		double xi2 = 0.;
	    if (xi == 1.)
	        xi2 = -1.;
	    else if (xi == -1.)
	        xi2 = 1.;
	    
        double[] out = new double[x.length]; 
        for (int i = 0; i < x.length; i++)
            out[i] = (x[i] - xi2)/(xi - xi2);
	    
        return out;
	}
	
	/* 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
//	static double[][] Mbar_2D(double xi, double[][] x) {
//		assert(xi == 1. || xi == -1.);
//		double xi2;
//	    if (xi == 1.)
//	        xi2 = -1.;
//	    else if (xi == -1.)
//	        xi2 = 1.;
//	    return (x - xi2)/(xi - xi2);
//	}
//	
	/* derivative of the 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double Mbar_prime(double xi, double x) {
		assert(xi == 1. || xi == -1.);
		double xi2 = 0.;
	    if (xi == 1.)
	        xi2 = -1.;
	    else if (xi == -1.)
	        xi2 = 1.;
	    return 1./(xi - xi2);
	}
	
	/* derivative of the 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
	static double[] Mbar_prime_1D(double xi, double[] x) {
		assert(xi == 1. || xi == -1.);
		double xi2 = 0.;
	    if (xi == 1.)
	        xi2 = -1.;
	    else if (xi == -1.)
	        xi2 = 1.;
	    
        double[] out = new double[x.length]; 
        for (int i = 0; i < x.length; i++)
            out[i] = 1./(xi - xi2);
	    
        return out;
	}
	
	/* derivative of the 1D linear shape function in the 1D canonical element [-1,1]
	 * xi is the coordinate where the function equals 1
	 * x is the coordinate where it is evaluated
	 */
//	static double[][] Mbar_prime_2D(double xi, double[][] x) {
//		assert(xi == 1. || xi == -1.);
//		double xi2;
//	    if (xi == 1.)
//	        xi2 = -1.;
//	    else if (xi == -1.)
//	        xi2 = 1.;
//	    return 1./(xi - xi2);
//	}
}    