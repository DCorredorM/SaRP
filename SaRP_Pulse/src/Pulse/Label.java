package Pulse;

public class Label {
	
public double prob;

public int cost;

	public Label(int pCost, double pProb) {
		prob = pProb;
		cost = pCost;
	}
	/**
	 * True if this dominates param
	 * @param pProb
	 * @return
	 */
	public boolean dominateLabel(int pCost, double pProb) {
		//System.out.println("At0: "+ attributes[0]);
		//System.out.println("At1: "+ attributes[1]);
		
		if(pProb < prob && pCost > cost) {
			return true;
		}else {
			return false;
		}
		
		
	}
	
	@Override
	public String toString() {
		return "";
	}
}
