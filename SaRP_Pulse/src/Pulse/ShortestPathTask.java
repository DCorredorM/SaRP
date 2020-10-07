/**
 * This class is used for parallelizing the SP algorithm.
 * 
 * Ref.: Lozano, L. and Medaglia, A. L. (2013). 
 * On an exact method for the constrained shortest path problem. Computers & Operations Research. 40 (1):378-384.
 * DOI: http://dx.doi.org/10.1016/j.cor.2012.07.008 
 * 
 * 
 * @author L. Lozano & D. Duque
 * @affiliation Universidad de los Andes - Centro para la Optimizaci�n y Probabilidad Aplicada (COPA)
 * @url http://copa.uniandes.edu.co/
 * 
 */
package Pulse;

public class ShortestPathTask implements Runnable {

	/**
	 * The dijkstra for distance
	 */
	private DukqstraDist spDist; 
	/**
	 * The dijkstra for time
	 */
	private DukqstraTime spTime;
	/**
	 * The dijkstra for Exptime
	 */
	private DukqstraExpTime spExpTime;
	/**
	 * A boolean indicator to know if it is for time or for distance
	 */
	private int algoRuning;
	/**
	 * The main method. Creates a task
	 * @param quienEs if it is for time or for distance
	 * @param dD
	 * @param dT
	 */
	public ShortestPathTask(int quienEs, DukqstraDist dD, DukqstraTime dT, DukqstraExpTime dET) {
		//quienES 1 Dist, 0 Time;
		algoRuning = quienEs;
		if(quienEs==1){
			spDist = dD;
		}else if (quienEs==0){
			spTime = dT;
		}else {
			spExpTime = dET;
		}
	}
	/**
	 * This method runs the desired dijkstra
	 */
	public void run() {
		// TODO Auto-generated method stub
		if(algoRuning==1){
			spDist.runAlgDist();
		}else if (algoRuning==0){
			spTime.runAlgTime();
		}else {
			spExpTime.runAlgExpTime();
		}
	}

	
	
}
