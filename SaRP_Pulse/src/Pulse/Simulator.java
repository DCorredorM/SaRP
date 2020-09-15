package Pulse;
import org.jblas.*;
import static org.jblas.DoubleMatrix.*;
import static org.jblas.MatrixFunctions.*;
import java.util.ArrayList;

public class Simulator {

	//Exponential variables
	/**
	 * 
	 * @param pRates
	 * @param path
	 * @return
	 */
	public static double calculateProb(ArrayList<Double> pRates,ArrayList<Double> path ) {
		double[] p = new double[pRates.size()];
		
		//System.out.println(pRates);
		//ToDO - iguales
		
		for(int i = 0;i<pRates.size();i++) {
			double rateAct = pRates.get(i);
			p[i] = 1;
			for(int j = 0;j<pRates.size();j++) {
				if(i!=j) {
					//System.out.println(i+" - "+rateAct + " - "+pRates.get(j) + " - "+(double)(rateAct/pRates.get(j)));
					p[i] = p[i]*(1-(double)(rateAct/pRates.get(j)));
				}
			}
			//System.out.println("p  "+i+"  " +p[i]);	
		}
		double suma = 0;
		
		for(int i = 0;i<pRates.size();i++) {
			suma += Math.exp(-pRates.get(i)*PulseGraph.TimeC)/p[i];
		}
		
		//System.out.println("suma: "+suma);
		
		/**System.out.println("Prob a: "+suma);
		System.out.println("Prob b: "+(calculateProbMatrix(pRates,path)));
		*/
		System.out.println("ME llamaron");
		return suma;
		
	}
	
	/**
	 * Devuelve la probabilidad de llegar tarde es decir P(T(p)>Tmax))
	 * @param pRates
	 * @param path
	 * @return
	 */
	public static double calculateProbMatrix(ArrayList<Double> pRates,ArrayList<Integer> path) {
		
		if (pRates.size()>=1) {
			DoubleMatrix t= DoubleMatrix.zeros(1,pRates.size());
			DoubleMatrix T=DoubleMatrix.zeros(pRates.size(),pRates.size());
			
			t.put(0,0,1);
			
			/**
			if(pRates.size()==1) {
				System.out.println("Antes");
				System.out.println(T);
				System.out.println(pRates.size());
				System.out.println(pRates.get(0));
				System.out.println("Path"+path);
				
			}
			*/
			
			//System.out.println(T);
			for(int i =0;i<=T.rows-1;i++) {
				if (i<T.rows-1) {
					T.put(i, i, -pRates.get(i));
					T.put(i, i+1, pRates.get(i));					
				}else {
					T.put(i, i, -pRates.get(i));
				}
			}
			
			/**System.out.println(T.rows+"x"+T.columns);
			
			for(int i =0;i<=T.rows-1;i++) {
				for(int j =0;j<=T.columns-1;j++) {
					System.out.print(T.get(i,j)+"\t");
				}
				System.out.println("");
			}
			System.out.println("");
			*/
			
			DoubleMatrix e_T=expm(T.mul(PulseGraph.TimeC));
			double prob= t.mmul(e_T.mmul(DoubleMatrix.ones(pRates.size(),1))).get(0,0);
			return prob;
			
		}else {
			return (0);
		}
	}
	
}




