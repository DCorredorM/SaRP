package PH_Fitting;

import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.fit.*;

import py4j.GatewayServer;

public class EMPhaseFitEntryPoint {
	private fitter fitter;
	
	public EMPhaseFitEntryPoint() {
		fitter = new fitter();
	}

	public fitter getEMHyperErlangFit() {
		return fitter;
	}

	public static void main(String[] args) {
		GatewayServer gatewayServer = new GatewayServer(new EMPhaseFitEntryPoint());
		gatewayServer.start();
		System.out.println("Gateway Server Started");
	}
}