/**
 * Data structure for an edge used for the SP implementation (not for the pulse!)
 * 
 * Ref.: Lozano, L. and Medaglia, A. L. (2013). 
 * On an exact method for the constrained shortest path problem. Computers & Operations Research. 40 (1):378-384.
 * DOI: http://dx.doi.org/10.1016/j.cor.2012.07.008 
 * 
 * 
 * @author D. Duque
 * @affiliation Universidad de los Andes - Centro para la Optimizaci�n y Probabilidad Aplicada (COPA)
 * @url http://copa.uniandes.edu.co/
 * 
 */

package Pulse;

import jphase.ContPhaseVar;



public class EdgePulse {
	
	/**
	 * The edge distance
	 */
	private int eDist;
	/**
	 * The edge time random variable
	 */
	private ContPhaseVar eTimeRV;
	
	/**
	 * The edge free folw speed time
	 */
	private int eMTime;
	
	
	/**
	 * The edge free folw speed time
	 */
	private int eMeanTime;
	
	/**
	 * 
	 */
	private EdgePulse nextE;
	/**
	 * The edge id
	 */
	private int id;
	/**
	 * The source node
	 */
	private VertexPulse source;
	/**
	 * The target node
	 */
	private VertexPulse target;
	
	/**
	 * This method creates an edge
	 * @param d the edge distance
	 * @param timeRV the edge time
	 * @param nT the edge tail
	 * @param nH the edge head
	 * @param nid the edge id
	 */
	public EdgePulse(int d ,int pMinTime, ContPhaseVar timeRV,  VertexPulse nT, VertexPulse nH, int nid) {
		// TODO Auto-generated constructor stub
		eDist = d;
		eTimeRV = timeRV;
		eMTime=pMinTime;
		eMeanTime=eMTime + (int) eTimeRV.moment(1);
		//System.out.println("min: "+eMTime);
		//System.out.println("meam: "+eMeanTime);
		this.source = nT;
		this.target = nH;
		this.id = nid;
	}
	
	/**
	 * This method
	 * @param e
	 */
	public void addNextCommonTailEdge(EdgePulse e){
		if(nextE!= null){
			nextE.addNextCommonTailEdge(e);
		}
		else{
			nextE = e;
		}
	}

	/**
	 * This method returns the next edge
	 * @return next edge
	 */
	public EdgePulse getNext()
	{
		return nextE;
	}
	
	/**
	 * This method changes the next edge
	 * @param e
	 */
	public void setNextE(EdgePulse e ){
		nextE = e;
	}
	/**
	 * This method returns the edge distance
	 * @return edge distance
	 */
	public int getWeightDist(){
		return eDist;
	}
	/**
	 * This method returns the edge free flow time
	 * @return edge time
	 */
	public int getWeightTime(){
		return eMTime;
	}
	
	/**
	 * This method returns the edge mean time
	 * @return edge time
	 */
	public int getWeightMeanTime(){
		return eMeanTime;
	}
	
	
	/**
	 * This method returns the edge free flow time
	 * @return edge time
	 */
	public double getWeightTimeRV(){
		return eMTime;
	}
	
	/**
	 * This method returns the tail
	 * @return edge tail
	 */
	public VertexPulse getSource(){
		return source;
	}
	
	/**
	 * This method returns the head
	 * @return edge head
	 */
	public VertexPulse getTarget(){
		return target;
	}
	/**
	 * This method returns the edge id
	 * @return
	 */
	public int getID()
	{
		return id;
	}
	/**
	 * This method finds the edge associated with a head node
	 * @param targetT
	 * @return
	 */
	public EdgePulse findEdgebyTarget( VertexPulse targetT)
	{
		if(targetT.getID() == this.target.getID())
		{
			return this;
		}else{
			if(nextE!= null)
			{
				return nextE.findEdgebyTarget(targetT);
			}
		}
		return null;
	}
	
	/**
	 * This method is the criteria to sort arcs
	 * @return the sum of the min dist and the min time
	 */
	public double getCompareCriteria(){
		return target.getMinDist() + target.getMinTime();
	}

	
}
