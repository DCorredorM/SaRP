package PH_Fitting;

import py4j.GatewayServer;

public class MomentsACPHFitEntryPoint {
	private Mfitter Mfitter;
	
	public MomentsACPHFitEntryPoint() {
		Mfitter = new Mfitter();
	}

	public Mfitter getACPHFit() {
		return Mfitter;
	}

	public static void main(String[] args) {
		GatewayServer gatewayServer = new GatewayServer(new MomentsACPHFitEntryPoint());
		gatewayServer.start();
		System.out.println("Gateway Server Started");
	}

}
