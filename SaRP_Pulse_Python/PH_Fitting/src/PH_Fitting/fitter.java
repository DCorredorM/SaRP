package PH_Fitting;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.fit.*;

import py4j.GatewayServer;
public class fitter {
	public static void main(String[] args) {
		System.out.println("Hola");
		
		//double[] data= {1,2,3,4,5,6};
		//EMPhaseFit v1=EMPhaseFit(data);
		//v1.fit();
	}
	public static EMHyperErlangFit EMHyperErlangFit(List<Double> data) {
		 double[] target = new double[data.size()];
		 for (int i = 0; i < target.length; i++) {
		    target[i] = data.get(i).doubleValue();  // java 1.4 style
		 }

		return new EMHyperErlangFit(target);
	}
	
	public static String hola(String a) {
		return a;
	}
}
