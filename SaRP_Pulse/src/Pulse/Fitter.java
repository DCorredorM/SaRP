package Pulse;

import java.io.BufferedReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;
import java.util.PrimitiveIterator.OfDouble;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.function.Min;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.stat.descriptive.summary.Sum;
import org.graalvm.compiler.core.common.type.ArithmeticOpTable.UnaryOp.Not;

import com.sun.org.apache.xalan.internal.xsltc.compiler.sym;
import com.sun.org.apache.xml.internal.security.keys.content.RetrievalMethod;
import com.sun.org.apache.xpath.internal.operations.And;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections; 
import java.util.List; 

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.MatrixUtils;
import jphase.fit.EMHyperErlangFit;
import no.uib.cipr.matrix.Vector;
import sun.security.util.Length;

//Reads DIMACS, simulates data and saves file with fitted Phase Type (alpha, T)
public class Fitter {

	/**
	 *Maximum number of phases
	 */
	static int N_Phase;

	/**
	 *Size of the sample
	 */
	static int NumSample;

	/**
	 * Number of arcs
	 */
	int NumArcs;
	/**
	 * Number of nodes
	 */
	static int NumNodes;
	/**
	 * Destination node
	 */
	int LastNode;
	/**
	 * Source node
	 */
	int Source;
	/**
	 * All the arcs in the network stored in a vector where Arcs[i][0]= Tail for arc i and Arcs[i][1]= Head for arc i 
	 */
	static int[][] Arcs;
	/**
	 * The distance attribute for any arc i
	 */
	static int[] Distance;
	/**
	 * The time attribute for any arc i
	 */
	static double[] Time;

	/**
	 * The time attribute for any arc i
	 */
	static double[] Mean;

	/**
	 * The free flow time attribute for any arc i
	 */
	static int[] MinTime;

	/**
	 * The data for each arc
	 */
	static double [][] Data;

	/**
	 * The time attribute for any arc i
	 */
	static ContPhaseVar[] TimeRV;

	/**
	 * Data structure for storing the graph
	 */
	private PulseGraph Gd;

	private int networkId;

	private String acro;
	private String acroData;

	static int numLabels;


	private ContPhaseVar ph;

	static Random r = new Random(0);

	//static double alpha = 0.90; //Confidence level of the chance constraint
	/**
	 * Read data from an instance
	 * @param numNodes
	 * @param numArcs
	 * @param sourceNode
	 * @param lastNode
	 * @param netId
	 */
	public Fitter(int numNodes, int numArcs, int sourceNode, int lastNode, int netId, String acronym,String acronymD,int n) {

		//Retrieves the information of the instance
		NumSample=200;
		N_Phase=n;
		NumArcs = numArcs;
		NumNodes = numNodes;
		LastNode = lastNode;
		Source = sourceNode;
		networkId = netId;
		acro = acronym;
		acroData=acronymD;

		//Creates the list of arcs. A list of distances and a list of times   --- Serian independientes del sentido de la red ! 
		Arcs = new int[numArcs][2]; // Array of array of lenght 2 of lenght numArcs 
		Distance = new int[numArcs]; //Array of int containing distance or cost
		Time = new double[numArcs]; //Array of double containing mean of travel time
		TimeRV = new ContPhaseVar [numArcs]; //Array of DenseContPhaseVar containig the RV of each time of lenght numArcs
		MinTime = new int[numArcs]; //Array of double containing free flow travel time of the arcs
		Mean = new double[numArcs]; //Array of double containing the mean of the travel time of the distribution 
		Data= new double[numArcs][NumSample];


		double [ ] [ ] A = new double[][] {{0}};
		double [ ]  tau = new double[] {0};
		ContPhaseVar ptRV=new DenseContPhaseVar(tau, A);
		ph=ptRV;
		//Creates the graph
		Gd = new PulseGraph(NumNodes);
	}


	/**
	 * This procedure creates the nodes for the graph
	 */
	public void upLoadNodes(){
		// All nodes are VertexPulse except the final node
		for (int i = 0; i < NumNodes; i++) {
			if(i!=(LastNode)){
				Gd.addVertex(new VertexPulse(i) ); //Primero lo creo, y luego lo meto. El id corresponde al n�mero del nodo
			}
		}
		// The final node is a FinalVertexPulse 
		FinalVertexPulse vv = new FinalVertexPulse(LastNode);
		Gd.addFinalVertex(vv);		
	}

	/**
	 * This procedure returns a graph
	 * @return the graph
	 */
	public PulseGraph getGd()
	{
		return Gd;
	}

	public static double mean(double[] m) {
		double sum = 0;
		for (int i = 0; i < m.length; i++) {
			sum += m[i];
		}
		return sum / m.length;
	}

	/**
	 * Checks whether a vector is sub-stochastic or not 
	 * @param a Vector to check
	 * @return true if the vector is sub-stochastic, false otherwise
	 */
	public static boolean checkSubStochasticVector(double[] a){
		double Epsilon=-10^(-10);
		boolean res = false; 
		boolean nonneg = true; 
		double cumsum = a[0]; 
		for (int i = 1; i < a.length; i++){
			if (a[i] < Epsilon){
				nonneg = false;
				break;
			}
			cumsum += a[i];    		
		}    	
		if (cumsum <= 1 + Epsilon && nonneg) res = true;  
		return res;
	}

	/**
	 * This procedure reads data from a data file in DIMACS format
	 * @throws NumberFormatException
	 * @throws IOException
	 */

	public void ReadDimacsF() throws NumberFormatException, IOException {
		File file2 = null;

		file2 = new File(acro);
		BufferedReader bufRdr2 = new BufferedReader(new FileReader(file2));
		String line2 = null;
		int row2 = 0;
		while ((line2 = bufRdr2.readLine()) != null && row2 < NumArcs) {			
			String[] Actual = line2.split("\t");	

			Arcs[row2][0] = Integer.parseInt(Actual[0])-1;
			Arcs[row2][1] =  Integer.parseInt(Actual[1])-1;
			Distance[row2] = (int)Math.floor(Double.parseDouble(Actual[2]));
			MinTime[row2]= (int)Math.floor(Double.parseDouble(Actual[3]));

			double[] tau;
			double [][] T;
			tau=string_to_list(Actual[4].split(", "));		


			T=new double[N_Phase][N_Phase];
			String[] TT= Actual[5].split("],");
			for(int i = 0; i<N_Phase;i++) {				
				T[i]=string_to_list(TT[i].split(", "));
			}			
			if (checkSubStochasticVector(tau)) {
				TimeRV[row2] = new DenseContPhaseVar(tau, T);
			}else {								
				double cumsum=0;
				for (double p : tau) {
					cumsum+=p;
				}

				for (int i = 0; i < tau.length; i++) {
					tau[i]=tau[i]/cumsum;
				}
				TimeRV[row2] = new DenseContPhaseVar(tau, T);
			}

			Mean[row2]=TimeRV[row2].moment(1);
			if (Mean[row2]<0) {
				System.out.println("problema "+ Mean[row2]);

			}

			//System.out.println("prob"+TimeRV[row2].cdf(PulseGraph.TimeC));
			row2++;


			/*
			for (int i=0;i<N_Phase;i++) {
				for (int j=0;j<N_Phase;j++) {
					System.out.print(T[i][j]+"\t");
				}
				System.out.println("");
			}
			System.out.println("");
			 */
			/*
			double var = Math.random()*10;		    
		    double phi=Math.sqrt(var+Time[row2]*Time[row2]);
		    double mu = Math.log(Time[row2]*Time[row2]/phi);
		    double sigma=Math.sqrt(Math.log(phi*phi/(Time[row2]*Time[row2])));
		    LogNormalDistribution logNormal = new LogNormalDistribution(mu,sigma);
		    double[] data=new double[NumSample];

		    for (int i = 0; i < data.length; i++) {
		        data[i]=logNormal.sample();
		    } 

		    Arrays.sort(data);
		    MinTime[row2]= data[0];
		    for (int i = 0; i < data.length; i++) {
		        data[i]=data[i]-MinTime[row2];
		    } 

		    EMHyperErlangFit fitter = new EMHyperErlangFit(data);
			ContPhaseVar v1 = fitter.fit(N_Phase);
			TimeRV[row2] = new DenseContPhaseVar(v1.getVectorArray(),v1.getMatrixArray());
			 */

			//System.out.println("min time "+MinTime[row2]);

			/*
		    System.out.println("mean "+logNormal.getNumericalMean());
		    System.out.println("var "+logNormal.getNumericalVariance());
		    System.out.println("min "+data[0]);
		    System.out.println("max "+data[data.length-1]);
		    System.out.println("Time "+Time[row2]);
		    System.out.println("mean "+mean(data));

		    for (int i = 0; i < data.length; i++) {
		        System.out.print(data[i]+",");;
		    }

		    System.out.println("Termine fitting");
			System.out.println(v1.description());
			v1.getVectorArray();
			System.out.println(fitter.getLogLikelihood());
			 */


		}		
	}
	public void ReadDimacsF_fit() throws NumberFormatException, IOException {
		System.out.println("entre");
		File file2 = null;
		file2 = new File("./networks/"+acro);
		BufferedReader bufRdr2 = new BufferedReader(new FileReader(file2));
		String line2 = null;
		int row2 = 0;
		System.out.println(NumArcs);
		while ((line2 = bufRdr2.readLine()) != null && row2 < NumArcs) {			
			String[] Actual = line2.split(";");
			//System.out.println(line2);
			Arcs[row2][0] = Integer.parseInt(Actual[0])-1;
			Arcs[row2][1] =  Integer.parseInt(Actual[1])-1;
			Distance[row2] = Integer.parseInt(Actual[2]);
			Time[row2]=Double.parseDouble(Actual[3]);

			double var = Math.random()*10;		    
			double phi=Math.sqrt(var+Time[row2]*Time[row2]);
			double mu = Math.log(Time[row2]*Time[row2]/phi);
			double sigma=Math.sqrt(Math.log(phi*phi/(Time[row2]*Time[row2])));
			LogNormalDistribution logNormal = new LogNormalDistribution(mu,sigma);
			double[] data=new double[NumSample];

			for (int i = 0; i < data.length; i++) {
				data[i]=logNormal.sample();
			} 

			Arrays.sort(data);
			MinTime[row2]=(int) Math.floor( data[0]);
			for (int i = 0; i < data.length; i++) {
				data[i]=data[i]-MinTime[row2];
			} 

			EMHyperErlangFit EMfit = new EMHyperErlangFit(data);
			ContPhaseVar v1 = EMfit.fit(N_Phase);
			//System.out.println(v1.description());
			TimeRV[row2] = new DenseContPhaseVar(v1.getVectorArray(),v1.getMatrixArray());

			Mean[row2]=Time[row2];

			//System.out.println(v1.description());
			//System.out.println("prob: ");
			//System.out.println(v1.cdf(85-MinTime[row2]));
			row2++;

			/*
			for (int i=0;i<N_Phase;i++) {
				for (int j=0;j<N_Phase;j++) {
					System.out.print(T[i][j]+"\t");
				}
				System.out.println("");
			}
			System.out.println("");
			 */
			/*

			 */

			//System.out.println("min time "+MinTime[row2]);

			/*
		    System.out.println("mean "+logNormal.getNumericalMean());
		    System.out.println("var "+logNormal.getNumericalVariance());
		    System.out.println("min "+data[0]);
		    System.out.println("max "+data[data.length-1]);
		    System.out.println("Time "+Time[row2]);
		    System.out.println("mean "+mean(data));

		    for (int i = 0; i < data.length; i++) {
		        System.out.print(data[i]+",");;
		    }

		    System.out.println("Termine fitting");
			System.out.println(v1.description());
			v1.getVectorArray();
			System.out.println(fitter.getLogLikelihood());
			 */


		}
	}

	/**
	 * Loads the data
	 * @throws IOException 
	 * @throws  
	 */
	public void loadData() throws IOException {		
		File file2 = new File(acroData);
		BufferedReader bufRdr2 = new BufferedReader(new FileReader(file2));

		String line2 = null;
		int row2 = 0;
		double[] tData=null;
		while ((line2 = bufRdr2.readLine()) != null && row2 < NumArcs) {
			String[] Actual = line2.split("\t");
			//System.out.println(line2);
			int t = Integer.parseInt(Actual[0])-1;
			int h =  Integer.parseInt(Actual[1])-1;
			Data[row2]=string_to_list(Actual[2].replace("[", "").replace("]", "").split(", "));
			row2++;
		}
	}

	public ContPhaseVar fitPath1(ArrayList<Integer> pPath) throws IOException {

		File file2 = new File(acroData);
		BufferedReader bufRdr2 = new BufferedReader(new FileReader(file2));

		String line2 = null;
		int row2 = 0;
		double[] tData=null;
		while ((line2 = bufRdr2.readLine()) != null && row2 < NumArcs) {
			String[] Actual = line2.split("\t");
			//System.out.println(line2);
			int t = Integer.parseInt(Actual[0])-1;
			int h =  Integer.parseInt(Actual[1])-1;

			if (inPath(t, h, pPath)) {
				double [] data = string_to_list(Actual[2].replace("[", "").replace("]", "").split(", "));
				if (tData!=null) {
					tData=sum(tData,data);
				}else {
					tData=data;
				}
			}
		}

		EMHyperErlangFit EMfit = new EMHyperErlangFit(tData);
		ContPhaseVar ph =  EMfit.fit(N_Phase);
		//System.out.println(v1.description());
		//TimeRV[row2] = new DenseContPhaseVar(v1.getVectorArray(),v1.getMatrixArray());


		return ph;		
	}


	public static ContPhaseVar fitPath(ArrayList<Integer> pPath)  {
		double[] tData=null;
		for (int i = 0; i < pPath.size()-1; i++) {
			if (tData!=null) {
				
				try {
					tData=sum(tData,Data[getArcId(pPath.get(i),pPath.get(i+1))]);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}else {
				try {
					tData=Data[getArcId(pPath.get(i),pPath.get(i+1))];
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		EMHyperErlangFit EMfit = new EMHyperErlangFit(tData);
		ContPhaseVar pha =  EMfit.fit(N_Phase);
		double [] tau=pha.getVectorArray();
		
		if (checkSubStochasticVector(tau)) {
			DenseContPhaseVar ph= new DenseContPhaseVar(tau, pha.getMatrixArray());
			return ph;
		}else {								
			double cumsum=0;
			for (double p : tau) {
				cumsum+=p;
			}

			for (int i = 0; i < tau.length; i++) {
				tau[i]=tau[i]/cumsum;
			}
			DenseContPhaseVar ph= new DenseContPhaseVar(tau, pha.getMatrixArray());
			return ph;
		}
		
		
				
	}
	
	public static int getArcId(int i,int j) throws Exception {
		int cont=0;
		boolean f=false;
		for (int[] js : Arcs) {
			if(i==js[0]&j==js[1]) {
				f=true;
				return cont;
			}
			cont++;
		}
		throw new Exception("Arc ("+i+", "+j+") Not found");
	}
	
	public static double[] sum(double[] d1, double[]d2) {
		double[] s=new double[d1.length];
		for (int i = 0; i < d1.length; i++) {
			s[i]=d1[i]+d2[i];			
		}
		return s;		
	}
	public boolean inPath(int i, int j, ArrayList<Integer>pPath) {

		for (int k = 0; k < pPath.size()-1; k++) {
			if (pPath.get(k)==i & pPath.get(k+1)==j) {
				return true;
			} 
		}
		return false;
	}

	public double[] string_to_list(String[] s) {
		double[] l=new double[s.length];
		for(int i=0;i < s.length;i++) {

			l[i]= Double.parseDouble(s[i].replace("[", "").replace("]", ""));
		}
		return l;
	}

	/**
	 * This procedure reads data from a data file in DIMACS format
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public void ReadDimacsB() throws NumberFormatException, IOException {
		File file2 = null;
		file2 = new File("./networks/"+acro);
		BufferedReader bufRdr2 = new BufferedReader(new FileReader(file2));
		String line2 = null;
		int row2 = 0;
		while ((line2 = bufRdr2.readLine()) != null && row2 < NumArcs) {
			String[] Actual = line2.split("\t");
			Arcs[row2][1] = Integer.parseInt(Actual[0]);
			Arcs[row2][0] =  Integer.parseInt(Actual[1]);
			Distance[row2] = Integer.parseInt(Actual[2]);
			Time[row2] = Double.parseDouble(Actual[3]);
			row2++;
		}
	}

}
