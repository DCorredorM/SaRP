/**
 * This class represents a node, contains the pulse main logic.
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

import org.graalvm.compiler.debug.CSVUtil.Escape;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.fit.EMHyperErlangFit;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;


public class VertexPulse {
	/**
	 * This array contains the indexes for all the outgoing arcs from this node
	 */
	ArrayList<Integer> magicIndex;

	/**
	 * Labels for dominance prunning
	 */

	ArrayList<Label> labels;

	/**
	 * Boolean that tells if the node is visited for first time
	 */
	boolean firstTime = true;
	/**
	 * Bounds to reach the end node
	 */
	int minDist;
	int maxTime;
	int minTime;
	int maxDist;



	/**
	 * SP stuff
	 */
	public static final int infinity = (int)Double.POSITIVE_INFINITY;
	/**
	 * The edge that is coming to the node
	 */
	private EdgePulse reverseEdges;
	/**
	 * The vertex id
	 */
	public int id;
	private VertexPulse leftDist;
	private VertexPulse rigthDist;
	private VertexPulse leftTime;
	private VertexPulse rigthTime;
	private boolean insertedDist;
	private boolean insertedTime;

	/**
	 * Creates a node
	 * @param i the id
	 */
	public VertexPulse(int i) {
		id = i;
		insertedDist = false;
		minDist = infinity;
		minTime = infinity;
		maxTime = 0;
		maxDist = 0;

		leftDist = this;
		rigthDist = this;
		leftTime = this;
		rigthTime = this;

		labels = new ArrayList<Label>();

		magicIndex = new ArrayList<Integer>();
	}

	/**
	 * Returns the node id
	 * @return
	 */
	public int  getID()
	{
		return id;
	}

	/**
	 * Adds an edge to the coming arcs list
	 * @param e the edge
	 */
	public void addReversedEdge(EdgePulse e)
	{
		if(reverseEdges!=null){
			reverseEdges.addNextCommonTailEdge(e);
		}else
			reverseEdges = e;
	}

	/**
	 * Returns the list of reversed edges
	 * @return
	 */
	public EdgePulse getReversedEdges() {
		if(reverseEdges!= null){
			return reverseEdges;
		}
		double [ ] alpha = new double[] {};
		double [ ] [ ] A = new double[][] {{}};
		ContPhaseVar PH=new DenseContPhaseVar(alpha, A );
		return new EdgePulse(1,1,PH, this,this , -1);
	}

	/**
	 * Sets the minimum distance
	 * @param c
	 */
	public void setMinDist(int c){
		minDist = c;
	}

	/**
	 * Returns the minimum distance
	 * @return
	 */
	public int getMinDist(){
		return minDist;
	}

	/**
	 * Sets the maximum time
	 * @param d
	 */
	public void setMaxTime(int d){
		maxTime = d;
	}
	/**
	 * Returns the maximum time
	 * @return
	 */
	public int getMaxTime(){
		return maxTime;
	}
	/**
	 * Sets the minimum time
	 * @param ntj
	 */
	public void setMinTime(int ntj){
		minTime = ntj;
	}
	/**
	 * Returns the minimum time
	 * @return
	 */
	public int getMinTime(){
		return minTime;
	}
	/**
	 * Sets the maximum distance
	 * @param md
	 */
	public void setMaxDist(int md){
		maxDist = md;
	}
	/**
	 * Returns the maximum distance
	 * @return
	 */
	public int getMaxDist(){
		return maxDist;
	}


	/**
	 * Unlink a vertex from the bucket
	 * @return true, if the buckets gets empty
	 */
	public boolean unLinkVertexDist(){
		if(rigthDist.getID() == id){
			leftDist=this;
			rigthDist=this;
			return true;
		}else{
			leftDist.setRigthDist(rigthDist);
			rigthDist.setLeftDist(leftDist);
			leftDist = this;
			rigthDist = this;
			return false;
		}
	}
	/**
	 * 
	 * @return
	 */
	public boolean unLinkVertexTime(){
		if(rigthTime.getID() == id){
			leftTime=this;
			rigthTime=this;
			return true;
		}else{
			leftTime.setRigthTime(rigthTime);
			rigthTime.setLeftTime(leftTime);
			leftTime = this;
			rigthTime = this;
			return false;
		}
	}

	public void fastUnlinkDist(){
		leftDist=this;
		rigthDist=this;
	}
	public void fastUnlinkTime(){
		leftTime = this;
		rigthTime = this;
	}
	public void unlinkRighBoundDist()
	{
		rigthDist = null;
	}
	public void unlinkRighBoundTime()
	{
		rigthTime = null;
	}
	/**
	 * Insert a vertex in a bucket. 
	 * New vertex is inserted on the left of the bucket entrance 
	 * @param v vertex in progress to be inserted
	 */
	public void insertVertexDist(VertexPulse v) {
		v.setLeftDist(leftDist);
		v.setRigthDist(this);
		leftDist.setRigthDist(v);
		leftDist = v;
	}

	/**
	 * 
	 * @param v
	 */
	public void insertVertexTime(VertexPulse v) {
		v.setLeftTime(leftTime);
		v.setRigthTime(this);
		leftTime.setRigthTime(v);
		leftTime = v;
	}

	/**
	 * Distance basic methods
	 */
	public void setLeftDist(VertexPulse v){
		leftDist= v;
	}
	public void setRigthDist(VertexPulse v){
		rigthDist= v;
	}
	public VertexPulse getBLeftDist(){
		return leftDist;
	}
	public VertexPulse getBRigthDist(){
		return rigthDist;
	}
	public void setInsertedDist(){
		insertedDist = true;
	}
	public boolean isInserteDist(){
		return insertedDist;
	}
	/**
	 * Time basic methods
	 */
	public void setLeftTime(VertexPulse v){
		leftTime= v;
	}
	public void setRigthTime(VertexPulse v){
		rigthTime= v;
	}
	public VertexPulse getBLeftTime(){
		return leftTime;
	}
	public VertexPulse getBRigthTime(){
		return rigthTime;
	}
	public void setInsertedTime(){
		insertedTime = true;
	}
	public boolean isInsertedTime(){
		return insertedTime;
	}




	/**
	 * 
	 */
	public void reset(){
		insertedDist = false;
	}

	// This is the pulse procedure
	/**
	 * This is the pulse procedure
	 * @param PTime
	 * @param PDist
	 * @param path
	 */
	public void pulse1(int pCost, ContPhaseVar ptRV,double pProb, double ptmin,double pMean,  ArrayList<Integer> path ) 
	{
		//System.out.println("Llegue con: "+pCost+" - "+pMean+" - "+path);
		// if a node is visited for first time, sort the arcs
		if (this.firstTime) {
			this.firstTime = false;
			this.Sort(this.magicIndex); 
			leftDist = null;
			rigthDist = null;
			reverseEdges = null;
		}

		// Label update
		changeLabels(pCost,pProb,path); 

		//Check the path completion:

		// Check for cycles
		if (PulseGraph.Visited[id]==0) {  

			// Add the node to the path
			path.add(id);

			// Update the visit indicator
			PulseGraph.Visited[id]=1;
			//System.out.println(path.toString());
			//System.out.println(DataHandler.Arcs[magicIndex.get(0)][0] + " - "+DataHandler.Arcs[magicIndex.get(0)][1] );
			//System.out.println(DataHandler.Arcs[magicIndex.get(1)][0] + " - "+DataHandler.Arcs[magicIndex.get(1)][1] );
			//System.out.println(DataHandler.Arcs[magicIndex.get(2)][0] + " - "+DataHandler.Arcs[magicIndex.get(2)][1] );
			//System.out.println(DataHandler.Arcs[magicIndex.get(3)][0] + " - "+DataHandler.Arcs[magicIndex.get(3)][1] );		

			// Pulse all the head nodes for the outgoing arcs only if getting to this node passes all the tests
			for (int i = 0; i < magicIndex.size(); i++) { 
				// Update distance and time

				int e = (Integer) magicIndex.get(i);
				int a = Fitter.Arcs[e][1];
				int newCost = (pCost + Fitter.Distance[e]);
				double newMean = pMean + Fitter.Mean[e];

				//System.out.println("viene");
				//System.out.println(ptRV.description());
				//System.out.println("se le suma "+e);
				//System.out.println(Fitter.TimeRV[e].description());



				ContPhaseVar newtRV = ptRV.sum(Fitter.TimeRV[e]);

				//System.out.println("queda");
				//System.out.println(newtRV.description());

				double newTmin=ptmin +Fitter.MinTime[e];					
				//double t=PulseGraph.TimeC-ptmin-PulseGraph.vertexes[a].getMinTime();
				//double newProb=newtRV.cdf(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())); //Coputes the probability of arriving on time to this node
				double t_bar=Math.max(0,PulseGraph.TimeC);
				double cota=newMean/Math.pow(t_bar,1);					
				double newProb=0;


				/*
					double iniTime = System.nanoTime();
					double newProb=newtRV.cdf(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())); //Coputes the probability of arriving on time to this node
					double finalTime = (System.nanoTime() - iniTime)/1000000000;
					System.out.println("No se usa la cota de "+cota+" y da "+ newProb+" y me deomro "+finalTime+" para una T de "+newtRV.getMatrixArray().length );
					System.out.println(1-newProb+"<="+cota);
				 */
				//
				//System.out.println("Llegue con2: "+newCost+" - "+newProb+ " - "+path + " - "+a);					

				// Pruning strategies: infeasibility, bounds and labels
				//System.out.println(PulseGraph.vertexes[a].CheckLabels(newCost,newMean,pRates,path,a));
				//System.out.println(checkInfeasibility(pRates,path));
				//System.out.println((newCost + PulseGraph.vertexes[a].getMinDist() < PulseGraph.PrimalBound));


				/*
				if (checkInfeasibility(newProb,path) && (newCost + PulseGraph.vertexes[a].getMinDist() < PulseGraph.PrimalBound) && !PulseGraph.vertexes[a].CheckLabels(newCost,newProb,path,a)){
					// If not pruned the pulse travels to the next head node
					PulseGraph.vertexes[a].pulse(newCost,newtRV, newProb,newTmin,newMean,path);										
				}
				 */
				if( (newCost + PulseGraph.vertexes[a].getMinDist() >= PulseGraph.PrimalBound)) {
					//System.out.println("Pode por cota");
					PulseGraph.Bound+=1;
				}

				if ((newCost + PulseGraph.vertexes[a].getMinDist() < PulseGraph.PrimalBound) ){
					PulseGraph.calc_needed+=1;
					if (cota<1-PulseGraph.alpha) {
						newProb=1-cota;
						//System.out.println("Se usa la cota de "+ newProb);
					}else if(t_bar<=0 ){
						//System.out.println("t_bar<0" + t_bar);
					}else {
						double iniTime = System.nanoTime();
						newProb=1-cota;
						//System.out.println("Empiezo a calcular la prob");
						try {
							newProb=newtRV.cdf(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())); //Coputes the probability of arriving on time to this node
						} catch (Exception e2) {
							// TODO: handle exception
							System.out.println("Este es el T: " + Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime()));
						}
						System.out.println("Este es el T: " + Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())+"\n"+"la prob es: "+newProb);


						/*
						try {
							System.out.println("Empiezo a calcular la prob");							
							newProb=newtRV.cdf(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())); //Coputes the probability of arriving on time to this node
						} catch (Exception e2) {

							System.out.println("No se pudo calcular la probabilidad");							
						}finally {

							System.out.println("No se pudo calcular la probabilidad");
						}
						 */

						double finalTime = (System.nanoTime() - iniTime)/1000000000;
						//System.out.println("No se usa la cota de "+cota+" y da "+ newProb+" y me deomro "+finalTime+" para una T de "+newtRV.getMatrixArray().length );
						//System.out.println(1-newProb+"<="+cota);
						PulseGraph.calc_made+=1;
					}

					try {
						newProb=newtRV.cdf(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime())); //Coputes the probability of arriving on time to this node
					} catch (Exception e2) {
						// TODO: handle exception
						System.out.println(Math.max(0,PulseGraph.TimeC-newTmin-PulseGraph.vertexes[a].getMinTime()));
					}


					if(!checkInfeasibility(newProb,path)) {
						//System.out.println("Pode por infact");
						PulseGraph.Infeasibility+=1;

					}
					if(!PulseGraph.vertexes[a].CheckLabels(newCost,newProb,path,a) ) {
						//System.out.println("Pode por dominancia");
						PulseGraph.Dominance+=1;
					}

					if (checkInfeasibility(newProb,path) && !PulseGraph.vertexes[a].CheckLabels(newCost,newProb,path,a) ){
						// If not pruned the pulse travels to the next head node
						if ((System.nanoTime()-PulseGraph.pulse_time)/1000000000<5000){
							PulseGraph.vertexes[a].pulse1(newCost,newtRV, newProb,newTmin,newMean,path);
						}

					}



				}


			}
			// Updates path and visit indicator for backtrack
			path.remove((path.size() - 1));
			PulseGraph.Visited[id]=0;
		}

	}


	/**
	 * 
	 * @param pCost
	 * @param ptRV
	 * @param pProb
	 * @param ptmin
	 * @param pMean
	 * @param path
	 */
	public void pulse2(int pCost, ContPhaseVar ptRV,double pProb, double ptmin,double pMean,  ArrayList<Integer> path ) {
		// if a node is visited for first time, sort the arcs
		if (this.firstTime) {
			this.firstTime = false;
			this.Sort(this.magicIndex); 
			leftDist = null;
			rigthDist = null;
			reverseEdges = null;
		}

		// Label update
		changeLabels(pCost,pProb,path); 

		// Check for cycles
		if (PulseGraph.Visited[id]==0) {

			// Add the node to the path
			path.add(id);			

			// Update the visit indicator
			PulseGraph.Visited[id]=1;

			// Pulse all the head nodes for the outgoing arcs only if getting to this node passes all the tests

			for (int i = 0; i < magicIndex.size(); i++) {

				// Update distance and time
				int e = (Integer) magicIndex.get(i); //Arc index
				int a = Fitter.Arcs[e][1];			// Head Node
				int newCost = (pCost + Fitter.Distance[e]);// New cost
				double newMean = pMean + Fitter.Mean[e];// New mean
				ContPhaseVar newtRV = ptRV.sum(Fitter.TimeRV[e]);// New time random variable.

				double newTmin=ptmin +Fitter.MinTime[e];// New free flow time					
				double newProb=0;				

				//Bound pruning
				if ((newCost + PulseGraph.vertexes[a].getMinDist() < PulseGraph.PrimalBound) ){
					newProb=calcProb(newtRV,newTmin,a);
					//Dominance pruning
					if(!PulseGraph.vertexes[a].CheckLabels(newCost,newProb,path,a) ) {
						//Infeasibility pruning
						if (checkInfeasibility(newProb,path)) {
							//If not pruned the pulse travels to the next head node
							if ((System.nanoTime()-PulseGraph.pulse_time)/1000000000<PulseGraph.pulseTimeLimit){
								PulseGraph.vertexes[a].pulse2(newCost,newtRV, newProb,newTmin,newMean,path);
							}
						}else {
							//System.out.println("Pode por dominancia");
							PulseGraph.Infeasibility+=1;
						}

					}else {
						//System.out.println("Pode por dominancia");
						PulseGraph.Dominance+=1;
					}
				}else {
					//System.out.println("Pode por cota");
					PulseGraph.Bound+=1;
				}
			}
			//Updates path and visit indicator for backtrack
			//System.out.println(path);
			path.remove((path.size() - 1));
			PulseGraph.Visited[id]=0;
		}				
	}



	public void pulse(int pCost, ContPhaseVar ptRV,double pProb, double ptmin,double pMean,  ArrayList<Integer> path,double[] pData) {
//		System.out.println("Llegue aca: "+id+" - "+pCost+" - "+pMean+" - "+PulseGraph.PrimalBound + " - "+pProb);
		// if a node is visited for first time, sort the arcs
		if (this.firstTime) {
			this.firstTime = false;
			this.Sort(this.magicIndex); 
			leftDist = null;
			rigthDist = null;
			reverseEdges = null;
		}

		// Label update
		changeLabels(pCost,pProb,path); 

		// Check for cycles
		if (PulseGraph.Visited[id]==0) {

			// Add the node to the path
			path.add(id);			

			// Update the visit indicator
			PulseGraph.Visited[id]=1;

			// Pulse all the head nodes for the outgoing arcs only if getting to this node passes all the tests

			for (int i = 0; i < magicIndex.size(); i++) {

				// Update distance and time
				int e = (Integer) magicIndex.get(i); //Arc index
				int a = Fitter.Arcs[e][1];			// Head Node
				int newCost = (pCost + Fitter.Distance[e]);// New cost
				double newMean = pMean + Fitter.Mean[e];// New mean
				
				
				
				
				double[] newData=Fitter.sum(pData, Fitter.Data[e]);							
				ContPhaseVar newtRV=newPHVar( ptRV, path,e,newData);								
				
				double newTmin=ptmin +Fitter.MinTime[e];// New free flow time					
				double newProb=0;
				
//				System.out.println("Las means son: "+newMean+"\t"+(newtRV.expectedValue())+" - "+Fitter.mean(newData));
				//Bound pruning
				if ((newCost + PulseGraph.vertexes[a].getMinDist() < PulseGraph.PrimalBound) ){
					newProb=calcProb(newtRV,newTmin,a);
					//Dominance pruning
					if(!PulseGraph.vertexes[a].CheckLabels(newCost,newProb,path,a) ) {
						//Infeasibility pruning
						if (checkInfeasibility(newProb,path)) {
							//If not pruned the pulse travels to the next head node
							if ((System.nanoTime()-PulseGraph.pulse_time)/1000000000<PulseGraph.pulseTimeLimit){								
								PulseGraph.vertexes[a].pulse(newCost,newtRV, newProb,newTmin,newMean,path,newData);
							}
						}else {
							//System.out.println("Pode por dominancia");
							PulseGraph.Infeasibility+=1;
						}

					}else {
						//System.out.println("Pode por dominancia");
						PulseGraph.Dominance+=1;
					}
				}else {
					//System.out.println("Pode por cota");
					PulseGraph.Bound+=1;
				}
			}
			//Updates path and visit indicator for backtrack
			//System.out.println(path);
			path.remove((path.size() - 1));
			PulseGraph.Visited[id]=0;
		}				
	}



	public static ContPhaseVar newPHVar(ContPhaseVar pPH,ArrayList<Integer> pPath,int arc,double [] pData) {


		if (pPath.size()%PulseGraph.refit==0) {	
			pPath.add(Fitter.Arcs[arc][1]);
			ContPhaseVar ph=Fitter.fitPath(pPath);
//			System.out.println("El pat es: "+pPath);
//			EMHyperErlangFit EMfit = new EMHyperErlangFit(pData);
//			ContPhaseVar v1=EMfit.fit(Fitter.N_Phase);
//			ContPhaseVar ph=new DenseContPhaseVar(v1.getVectorArray(),v1.getMatrixArray());
//			System.out.println("Pase por aca!!");
			pPath.remove(pPath.size()-1);
			return ph; 

		} else {
//			System.out.println("Pase por aca...");
			return pPH.sum(Fitter.TimeRV[arc]);// New time random variable.
		}

	}


	/**
	 * Computes the probability of arriving on time to node pHeadNode
	 * @param pTimeRV
	 * @param pTMin
	 * @param pHeadNode
	 * @return
	 */
	public static double calcProb(ContPhaseVar pTimeRV ,Double pTMin,int pHeadNode) {
		double prob=0;
//		System.out.println("Esre es el T con el que calculo: "+(PulseGraph.TimeC-pTMin-PulseGraph.vertexes[pHeadNode].getMinTime())+"\t\t\n"+(PulseGraph.TimeC)+"\t"+(pTMin)+"\t"+(PulseGraph.vertexes[pHeadNode].getMinTime()));
		int n=pTimeRV.getNumPhases();		
		double trace=0;
		for (int i=0;i<n;i++) {
			trace+=pTimeRV.getMatrix().get(i, i);
		}
//		System.out.println(pTimeRV.toString()+"\n"+trace);
		try {			
			if(-1*trace<3000) {
			prob=pTimeRV.cdf(Math.max(0,PulseGraph.TimeC-pTMin-PulseGraph.vertexes[pHeadNode].getMinTime())); //Coputes the probability of arriving on time to this node			
			}
//			System.out.println(prob);
		} catch (Exception e) {
			// TODO: handle exception
		}
		return prob;
		//return 0.95;
	}

	/**
	 * This method sorts the arcs
	 * @param set
	 */
	private void Sort(ArrayList<Integer> set) 
	{
		QS(magicIndex, 0, magicIndex.size()-1);
	}

	/**
	 * 
	 * @param e
	 * @param b
	 * @param t
	 * @return
	 */
	public int colocar(ArrayList<Integer> e, int b, int t)
	{
		int i;
		int pivote, valor_pivote;
		int temp;

		pivote = b;
		valor_pivote = PulseGraph.vertexes[Fitter.Arcs[e.get(pivote)][1]].getCompareCriteria();

		for (i=b+1; i<=t; i++){
			if (PulseGraph.vertexes[Fitter.Arcs[e.get(i)][1]].getCompareCriteria() < valor_pivote){
				pivote++;    
				temp= e.get(i);
				e.set(i, e.get(pivote));
				e.set(pivote, temp);
			}
		}
		temp=e.get(b);
		e.set(b, e.get(pivote));
		e.set(pivote, temp);
		return pivote;
	} 


	/**
	 * 
	 * @param e
	 * @param b
	 * @param t
	 */
	public void QS(ArrayList<Integer> e, int b, int t)
	{
		int pivote;
		if(b < t){
			pivote=colocar(e,b,t);
			QS(e,b,pivote-1);
			QS(e,pivote+1,t);
		}  
	}
	/**
	 * This method verifies dominance
	 * @param PTime
	 * @param PDist
	 * @return true if the label (pProb,pCost) is dominated
	 */
	public boolean CheckLabels(int pCost,double pProb ,ArrayList<Integer> pPath,int node)
	// Label pruning strategy
	{
		for (int i = 0; i < PulseGraph.vertexes[node].labels.size(); i++) {
			if (PulseGraph.vertexes[node].labels.get(i).dominateLabel(pCost,pProb)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * 
	 * @param pProb
	 * @param pPath
	 * @return True if pPath is feasible i.e. pProb<alpha
	 */
	public boolean checkInfeasibility(double pProb,ArrayList<Integer> pPath) {

		if(pProb < PulseGraph.alpha) {
			return false;
		}
		return true;
	}

	/**
	 * This method modifies the labels.
	 * @param PTime the current pulse time
	 * @param PDist the current pulse distance
	 */
	private void changeLabels(int pCost, double pProb,ArrayList<Integer> pPath) {
		labels.add(new Label(pCost,pProb));

		//To - Do
		/**
		if (labels.size() == 0) {
			labels.add(new Label(pCost,attributes));
		} else if (labels.size() == 1) {
			labels.add( < labels.get(0).attributes[0] ? 0 : 1, new Label(pCost,attributes));
		} else if (labels.size() == 2) {
			labels.add(objs[1] < labels.get(1).attributes[1] ? 0 : 1, new Label(pCost,attributes));
		}else {
			if (labels.size() < DataHandler.numLabels) {
				insertLabel1(objs);
			} else {
				if(DataHandler.numLabels > 2) {
					int luck = DataHandler.r.nextInt(DataHandler.numLabels-2)+2;
					labels.remove(luck);
					insertLabel1(objs);
				}
			}
		}
		 */

	}


	/**
	 * This method inserts a label based on their time
	 * @param objs
	 */

	/**
	private void insertLabel1(double[] objs) {
		Label np = new Label(objs);
		double cScore = np.attributes[0];

		boolean cond = true;
		int l = 0;
		int r = labels.size();
		int m = (int) ((l + r) / 2);
		double mVal = 0;
		if (labels.size() == 0) {
			labels.add(np);
			return;
		} else if (labels.size() == 1) {
			mVal = labels.get(m).attributes[0];;
			labels.add(cScore < mVal ? 0 : 1, np);
			return;
		} else {
			mVal = labels.get(m).attributes[0];
		}

		while (cond) {
			if (r - l > 1) {
				if (cScore < mVal) {
					r = m;
					m = (int) ((l + r) / 2);
				} else if (cScore > mVal) {
					l = m;
					m = (int) ((l + r) / 2);
				} else {
					labels.set(m, np);
					return;
				}
				mVal = labels.get(m).attributes[0];;
			} else {
				cond = false;
				if (l == m) {
					labels.add(cScore < mVal ? l : l + 1, np);
				} else if (r == m) {
					System.out.println("esto no pasa ");
					labels.add(cScore < mVal ? r : Math.min(r + 1, labels.size()), np);
				} else {
					System.err.println(VertexPulse.class +  " insert label, error");
				}
				return;

			}
		}

	}
	 */

	/**
	 * 
	 * @return
	 */
	public int getCompareCriteria(){
		if (PulseGraph.PrimalBound<9999999) {
			return getMinDist();
		}else {
			System.out.println("Entre aca, minTime es "+getMinTime());
			return getMinTime();
		}
	}


}

