/**
 * This class represents the final node in a network. The pulse procedure is overrided!
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

import java.util.ArrayList;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;



public class FinalVertexPulse extends VertexPulse {
	
	/**
	 * Node id
	 */
	public int id;
	
	/**
	 * SP stuff
	 */
	private VertexPulse bLeft;
	private VertexPulse bRigth;
	private EdgePulse reverseEdges;
	
	/**
	 * Bounds
	 */
	int minDist;
	int maxTime;
	int minTime;
	int maxDist;
	int minExpTime;
	int maxExpTime;
	
	/**
	 * 
	 */
	private boolean inserted;
	/**
	 * 
	 */
	int c=0;
	int d=0;
	/**
	 * Creates the final node
	 * @param i
	 */
	public FinalVertexPulse(int i) {
		super(i);
		id = i;
		inserted = false;
		minDist = infinity;
		minTime = infinity;
		maxTime = 0;
		maxDist = 0;
		bLeft = this;
		bRigth = this;
	}

	/**
	 * Returns the node id
	 */
	public int  getID()
	{
		return id;
	}
	
	
	/**
	 * Adds an edge
	 */
	public void addReversedEdge(EdgePulse e)
	{
		if(reverseEdges!=null){
			reverseEdges.addNextCommonTailEdge(e);
		}else
			reverseEdges = e;
	}
	
	
	public EdgePulse findEdgeByTarget(VertexPulse target){
		if(reverseEdges!=null){
			reverseEdges.findEdgebyTarget(target);
		}
		return null;
	}

	public EdgePulse getReversedEdges() {
		if(reverseEdges!= null){
			return reverseEdges;
		}
		double [ ] alpha = new double[] {};
		double [ ] [ ] A = new double[][] {{}};
		ContPhaseVar PH=new DenseContPhaseVar(alpha, A );
		return new EdgePulse(1,1,PH, this,this , -1);
	}
	
	public void setMinDist(int c){
		minDist = c;
	}
	

	public int getMinDist(){
		return minDist;
	}
	
	public void setMaxTime(int mt){
		maxTime = mt;
	}
	
	public int getMaxTime(){
		return maxTime;
	}
	public void setMaxExpTime(int mt){
		maxTime = mt;
	}
	
	public int getMaxExpTime(){
		return maxTime;
	}
	
	public void setMinTime(int t){
		minTime = t;
	}
	
	public int getMinTime(){
		return minTime;
	}
	public void setMinExpTime(int t){
		minTime = t;
	}
	
	public int getMinExpTime(){
		return minTime;
	}
	
	public void setMaxDist(int md){
		maxDist = md;
	}
	
	public int getMaxDist(){
		return maxDist;
	}
	
	
	/**
	 * Unlink a vertex from the bucket
	 * @return true, if the buckets gets empty
	 */
	public boolean unLinkVertexDist(){
		if(bRigth.getID() == id){
			bLeft=this;
			bRigth=this;
			return true;
		}else{
			bLeft.setRigthDist(bRigth);
			bRigth.setLeftDist(bLeft);
			bLeft = this;
			bRigth = this;
			return false;
		}
	}
	/**
	 * Insert a vertex in a bucket. 
	 * New vertex is inserted on the left of the bucket entrance 
	 * @param v vertex in progress to be inserted
	 */
	public void insertVertexDist(VertexPulse v) {
		v.setLeftDist(bLeft);
		v.setRigthDist(this);
		bLeft.setRigthDist(v);
		bLeft = v;
	}
	
	
	public void setLeftDist(VertexPulse v){
		bLeft= v;
	}
	public void setRigthDist(VertexPulse v){
		bRigth= v;
	}
	public VertexPulse getBLeftDist(){
		return bLeft;
	}
	public VertexPulse getBRigthDist(){
		return bRigth;
	}
	public void setInsertedDist(){
		inserted = true;
	}
	public boolean isInserteDist(){
		return inserted;
	}
	/**
	 * Unlink a vertex from the bucket
	 * @return true, if the buckets gets empty
	 */
	public boolean unLinkVertexTime(){
		if(bRigth.getID() == id){
			bLeft=this;
			bRigth=this;
			return true;
		}else{
			bLeft.setRigthTime(bRigth);
			bRigth.setLeftTime(bLeft);
			bLeft = this;
			bRigth = this;
			return false;
		}
	}
	/**
	 * Insert a vertex in a bucket. 
	 * New vertex is inserted on the left of the bucket entrance 
	 * @param v vertex in progress to be inserted
	 */
	public void insertVertexTime(VertexPulse v) {
		v.setLeftTime(bLeft);
		v.setRigthTime(this);
		bLeft.setRigthTime(v);
		bLeft = v;
	}
	
	
	public void setLeftTime(VertexPulse v){
		bLeft= v;
	}
	public void setRigthTime(VertexPulse v){
		bRigth= v;
	}
	public VertexPulse getBLeftTime(){
		return bLeft;
	}
	public VertexPulse getBRigthTime(){
		return bRigth;
	}
	public void setInsertedTime(){
		inserted = true;
	}
	public boolean isInserteTime(){
		return inserted;
	}

	/**
	 * Unlink a vertex from the bucket
	 * @return true, if the buckets gets empty
	 */
	public boolean unLinkVertexExpTime(){
		if(bRigth.getID() == id){
			bLeft=this;
			bRigth=this;
			return true;
		}else{
			bLeft.setRigthExpTime(bRigth);
			bRigth.setLeftExpTime(bLeft);
			bLeft = this;
			bRigth = this;
			return false;
		}
	}
	/**
	 * Insert a vertex in a bucket. 
	 * New vertex is inserted on the left of the bucket entrance 
	 * @param v vertex in progress to be inserted
	 */
	public void insertVertexExpTime(VertexPulse v) {
		v.setLeftExpTime(bLeft);
		v.setRigthExpTime(this);
		bLeft.setRigthExpTime(v);
		bLeft = v;
	}
	
	
	public void setLeftExpTime(VertexPulse v){
		bLeft= v;
	}
	public void setRigthExpTime(VertexPulse v){
		
		bRigth= v;
	}
	public VertexPulse getBLeftExpTime(){
		return bLeft;
	}
	public VertexPulse getBRigthExpTime(){
		return bRigth;
	}
	public void setInsertedExpTime(){
		inserted = true;
	}
	public boolean isInserteExpTime(){
		return inserted;
	}


	public void reset(){
		inserted = false;
	}
	public void setBounds(int MT, int MD){
		maxDist = MD- minDist;
		maxTime = MT - minTime;
	}
	/**
	 * This is the pulse for the last node !!!
	 * If the path is feasible and updates the primal bound the information on the best solution found is updated
	 */
	public void pulse(int pCost, ContPhaseVar ptRV,double pProb, double ptmin,double pMean,  ArrayList<Integer> path,double[] pData){		
		// Add node id to path
		path.add(id);
//		System.out.println("Llegue aca (final) "+id+" - "+pCost+" - "+pMean+" - "+PulseGraph.PrimalBound + " - "+pProb);
		// If the path is feasible and updates the primal bound the information on the best solution found is updated
		if(pCost<PulseGraph.PrimalBound && pProb >= PulseGraph.alpha){
			// Update the best solution known
			PulseGraph.Path.clear();
			PulseGraph.Path.addAll(path);
			PulseGraph.PrimalBound=pCost;
			PulseGraph.Mean=pMean;
			PulseGraph.minTime=ptmin;
			
			//The best distance
			
			PulseGraph.PH = ptRV;
			//PulseGraph.finalProb = ptRV.cdf(PulseGraph.TimeC-ptmin);
			PulseGraph.finalProb=pProb;
		}
		// Remove the node id to backtrack
		path.remove((path.size()-1));
	}	
}
