/**
 * This is the main class for the pulse algorithm. 
 * 
 * Ref.: Lozano, L. and Medaglia, A. L. (2013). 
 * On an exact method for the constrained shortest path problem. Computers & Operations Research. 40 (1):378-384.
 * DOI: http://dx.doi.org/10.1016/j.cor.2012.07.008 
 * 
 * 
 * @author L. Lozano & D. Duque & N. Cabrera
 * @affiliation Universidad de los Andes - Centro para la Optimizaci�n y Probabilidad Aplicada (COPA)
 * @url http://copa.uniandes.edu.co/
 * 
 */

package Pulse;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.nio.file.Path;
import java.nio.file.Paths;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import sun.java2d.pipe.AlphaColorPipe;

public class PulseMain {

	
	public static void main(String[] args) throws IOException, InterruptedException {
		
		
		try {
			String[] arg= {args[0],args[1]};		
			runPulse(arg);			
		} catch (Exception e) {
			String city="Chicago-Sketch";
			String city_file=Paths.get("../data/Networks"+"/"+city).normalize().toString();
//			System.out.println(Paths.get(city_file).toAbsolutePath().normalize());
			String[] arg= {city_file,city};		
			runPulse(arg);
		}
		
		
	}
	
	/**
	 * This method runs the pulse
	 * @param args
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void runPulse(String[] args) throws IOException, InterruptedException {
		
		String city_file=args[0];
		String city=args[1];
		
		File config_File = new File(city_file+"/"+city+"_config.txt");
		BufferedReader bufRedr = new BufferedReader(new FileReader(config_File));
				
		//Reads the file
		String actLine = null;
		String [] information = new String [9];
		int rowA = 0;
		int colA = 0;

		while((actLine = bufRedr.readLine()) != null && rowA < 9){	
			String [] info1 = actLine.split(":");			
			information[rowA] = info1[1];
			//System.out.println(rowA+" "+info1[0]+" "+info1[1]);
			rowA++;
		}
		double T_max=Double.parseDouble(information[3]);
		int source =Integer.parseInt(information[4]);
		int target =Integer.parseInt(information[5]);
		int num_nodes=Integer.parseInt(information[2]);
		int num_arcs=Integer.parseInt(information[1]);
		int N_Phases=Integer.parseInt(information[6]);
		double alpha=Double.parseDouble(information[7]);
		double timeLimit=Double.parseDouble(information[8]);


		String net_file = city_file +"/"+information[0];
		//Choose pulse direction: 1: Forward -  2: Backward

		int direction = 1;

		//Initialize the graph

		PulseGraph network = null;

		Fitter dataA = null;
		//Forward direction: 
		if(direction == 1) {
			dataA = new Fitter(num_nodes,num_arcs,source,target,0,net_file,N_Phases);
			dataA.ReadDimacsF();
			network = S_createGraph(dataA);

		}
		//Backward direction: 
		if(direction == 2) {
			dataA = new Fitter(num_nodes,num_arcs,source,target,0,net_file,N_Phases );
			dataA.ReadDimacsF();
			network = S_createGraph(dataA);
		}

		//Initial bounds

		SP(dataA,network);

		//Time limit tightness


		network.SetConstraint(T_max);
		network.Setaplpha(alpha);


		//Pulse algorithm parameters:

		//First primal bound

		//int MD=network.getVertexByID(dataA.Source-1).getMaxDist();
		int MD = 9999999;
		double prob=0;
		network.setPrimalBound(MD);
		network.setFinalProb(0.0);
		//network.TimeStar = network.getVertexByID(dataA.Source-1).getMinTime();

		//Size of Q

		dataA.numLabels = 10;

		//Initial weights

		int[] initialWeights = new int[2];
		initialWeights[0] = 0;
		initialWeights[1] = 0;

		//Initial path

		ArrayList<Integer> finalPath = new ArrayList();
		//Run the pulse and recover the path:

		//Starts the clock

		double iniTime = System.nanoTime(); 
		network.pulseTimeLimit=timeLimit;
		network.pulse_time= System.nanoTime();		
		double [ ] [ ] A = new double[][] {{0}};
		double [ ]  tau = new double[] {0};
		ContPhaseVar ptRV=new DenseContPhaseVar(tau, A);
		network.getVertexByID(dataA.Source).pulse(0, ptRV, 1, 0,0, finalPath);

		/*
		try {
			network.getVertexByID(dataA.Source).pulse(0, ptRV, 1, 0,0, finalPath);
		}
		catch (Exception e) {
			System.out.println(e);
		}
		 */
		//finalPath = completeThePath(dataA,network);
		finalPath = network.Path;

		//Ends the clock

		double finalTime = (System.nanoTime() - iniTime)/1000000000;


		//Print results

		String text =source+"\t"+target+"\t"+ finalTime+"\t"+network.PrimalBound+"\t"+network.minTime+"\t"+ network.finalProb+"\t"+network.Bound+"\t"+network.Infeasibility+"\t"+network.Dominance+"\t"+finalPath;
		System.out.println(text);
		String results_path=city_file+"/Results/"+information[0];

//		System.out.println(results_path);
		try(FileWriter fw = new FileWriter(results_path, true);
				BufferedWriter bw = new BufferedWriter(fw);
				PrintWriter out = new PrintWriter(bw))
		{
			out.println(text);
			//more code

		} catch (IOException e) {
			//exception handling left as an exercise for the reader
		}

	}


	/**
	 * This method runs the pulse
	 * @param args
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void main1(String[] args) throws IOException, InterruptedException {

		//The main file of the networks
		//String Mfile= "C:/Users/d.corredor/OneDrive - Universidad de los andes/Thesis/3. Case_study/Chicago/Networks";
		//String Mfile= "C:/Users/d.corredor/OneDrive - Universidad de los andes/Thesis/3. Case_study/Transportation_networks/Networks";
		String Mfile="/Users/davidcorredor/Documents/PulseAlgorithm/SaRP/data/Networks";

		//Networks to run.
		String [] cities=new String[] {"Chicago-Sketch"};//"SiouxFalls","Chicago-Sketch","GoldCoast","Philadelphia","Sydney"

		//Number of instances to run:
		int[] instances_run =new int[] {1,2,3,4,5,6,7};//,5,6,12,13,14,16,19};
		//Creates the networks and runs the insatnces for every city

		for (int i = 0; i < cities.length; i++) {
			String city= cities[i];
			System.out.println("Empece con "+city);		
			String city_file = Mfile+"/"+city;//The cities folder path

			//Opens the file
			File inst_File = new File(city_file+"/"+city+"_Instances_mod.txt");
			BufferedReader bufRedrr = new BufferedReader(new FileReader(inst_File));

			// Reads and runs all instances:
			String actLine1 = null;
			String [] instance = new String [10];
			int row1 = 0;

			int source;
			int target;
			double T_max;
			int N_Phases;

			while((actLine1 = bufRedrr.readLine()) != null && row1 < 42 ){	
				//in(row1,instances_run)
				if (row1>-1 && in(row1,instances_run) ) {
					String [] info = actLine1.split("\t");
					source=(int) Double.parseDouble(info[0])-1;
					target=(int) Double.parseDouble(info[1])-1;
					T_max= Double.parseDouble(info[2]);

					System.out.println(source+"\t"+target+"\t"+T_max);
					//Opens the file
					File config_File = new File(city_file+"/"+city+"_config.txt");
					BufferedReader bufRedr = new BufferedReader(new FileReader(config_File));

					//Reads the file

					String actLine = null;
					String [] information = new String [7];
					int rowA = 0;
					int colA = 0;

					while((actLine = bufRedr.readLine()) != null && rowA < 7){	
						String [] info1 = actLine.split(":");
						information[rowA] = info1[1];
						//System.out.println(rowA+" "+info1[0]+" "+info1[1]);
						rowA++;
					}
					int num_nodes=Integer.parseInt(information[2]);
					int num_arcs=Integer.parseInt(information[1]);
					N_Phases=Integer.parseInt(information[6]);



					String net_file = city_file +"/"+information[0];
					//Choose pulse direction: 1: Forward -  2: Backward

					int direction = 1;

					//Initialize the graph

					PulseGraph network = null;

					Fitter dataA = null;
					//Forward direction: 
					if(direction == 1) {
						dataA = new Fitter(num_nodes,num_arcs,source,target,i,net_file,N_Phases);
						dataA.ReadDimacsF();
						network = S_createGraph(dataA);

					}
					//Backward direction: 
					if(direction == 2) {
						dataA = new Fitter(num_nodes,num_arcs,source,target,i,net_file,N_Phases );
						dataA.ReadDimacsF();
						network = S_createGraph(dataA);
					}

					//Initial bounds

					SP(dataA,network);

					//Time limit tightness


					network.SetConstraint(T_max);
					network.Setaplpha(0.9);


					//Pulse algorithm parameters:

					//First primal bound

					//int MD=network.getVertexByID(dataA.Source-1).getMaxDist();
					int MD = 9999999;
					double prob=0;
					network.setPrimalBound(MD);
					network.setFinalProb(0.0);
					//network.TimeStar = network.getVertexByID(dataA.Source-1).getMinTime();

					//Size of Q

					dataA.numLabels = 10;

					//Initial weights

					int[] initialWeights = new int[2];
					initialWeights[0] = 0;
					initialWeights[1] = 0;

					//Initial path

					ArrayList<Integer> finalPath = new ArrayList();
					//Run the pulse and recover the path:

					//Starts the clock

					double iniTime = System.nanoTime(); 
					network.pulse_time= System.nanoTime(); 
					double [ ] [ ] A = new double[][] {{0}};
					double [ ]  tau = new double[] {0};
					ContPhaseVar ptRV=new DenseContPhaseVar(tau, A);
					network.getVertexByID(dataA.Source).pulse(0, ptRV, 1, 0,0, finalPath);

					/*
					try {
						network.getVertexByID(dataA.Source).pulse(0, ptRV, 1, 0,0, finalPath);
					}
					catch (Exception e) {
						System.out.println(e);
					}
					 */
					//finalPath = completeThePath(dataA,network);
					finalPath = network.Path;

					//Ends the clock

					double finalTime = (System.nanoTime() - iniTime)/1000000000;


					//Print results

					String text =source+"\t"+target+"\t"+ finalTime+"\t"+network.PrimalBound+"\t"+network.minTime+"\t"+ network.finalProb+"\t"+network.Bound+"\t"+network.Infeasibility+"\t"+network.Dominance+"\t"+finalPath;
					System.out.println(text);
					String results_path=city_file+"/Results/"+information[0];


					try(FileWriter fw = new FileWriter(results_path, true);
							BufferedWriter bw = new BufferedWriter(fw);
							PrintWriter out = new PrintWriter(bw))
					{
						out.println(text);
						//more code

					} catch (IOException e) {
						//exception handling left as an exercise for the reader
					}

				}				
				row1++;
			}
		}



		/*


		//Create the test instance:

			File testFile = new File("./networks/Config2.txt");

			BufferedReader bufRedr = new BufferedReader(new FileReader(testFile));

			String actLine = null;

			String [] information = new String [6];

			int rowA = 0;
			int colA = 0;

			while((actLine = bufRedr.readLine()) != null && rowA < 6){	
				String [] info = actLine.split(":");
				information[rowA] = info[1];
				rowA++;
			}


			//Choose pulse direction: 1: Forward -  2: Backward

			int direction = 1;

			//Initialize the graph

			PulseGraph network = null;
			//DataHandler dataA = null;
			Fitter dataA = null;


		 */
		/*
		//Forward direction: 

			if(direction == 1) {

				/**
				dataA = new DataHandler(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[4]),Integer.parseInt(information[5]),1,information[0]);
				dataA.ReadDimacsF();
				network = createGraph(dataA);	
				//coment!
				dataA = new Fitter(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[4]),Integer.parseInt(information[5]),1,information[0]);
				dataA.ReadDimacsF();
				network = S_createGraph(dataA);

			}


		//Backward direction: 

			if(direction == 2) {
				/*
				dataA = new DataHandler(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[5]),Integer.parseInt(information[4]),1,information[0]);	
				dataA.ReadDimacsB();
				network = createGraph(dataA);
				coment

				dataA = new Fitter(Integer.parseInt(information[2]),Integer.parseInt(information[1]),Integer.parseInt(information[4]),Integer.parseInt(information[5]),1,information[0]);
				dataA.ReadDimacsF_fit();
				network = S_createGraph(dataA);
			}

		//Initial bounds

			SP(dataA,network);
		//Time limit tightness

			double k = 5.5;


			//Time limit
			System.out.println(network.getVertexByID(dataA.Source).getMaxTime());
			System.out.println(network.getVertexByID(dataA.Source).getMinTime());
			//double timeLimit = (double) ((network.getVertexByID(dataA.Source).getMaxTime() - network.getVertexByID(dataA.Source).getMinTime())*k + network.getVertexByID(dataA.Source).getMinTime());
			double timeLimit = 12271;
			System.out.println("tl: "+timeLimit);
			//System.out.println("tl'"+((13015-9310)*k+9310));
			//double timeLimit = 150;
			network.SetConstraint(timeLimit);

		//Pulse algorithm parameters:

			//First primal bound

			//int MD=network.getVertexByID(dataA.Source-1).getMaxDist();
			int MD = 9999999;
			network.setPrimalBound(MD);
			//network.TimeStar = network.getVertexByID(dataA.Source-1).getMinTime();

			//Size of Q

			dataA.numLabels = 10;

			//Initial weights

			int[] initialWeights = new int[2];
			initialWeights[0] = 0;
			initialWeights[1] = 0;

			//Initial path

			ArrayList<Integer> finalPath = new ArrayList();

			//for (int i = 0; i < network.getNumNodes(); i++) {
			//	System.out.println("t "+i+": "+network.vertexes[i].getMinTime());
			//	System.out.println("d "+i+": "+network.vertexes[i].getMinDist());
			//} 	

		//Run the pulse and recover the path:

			//Starts the clock

			double iniTime = System.nanoTime(); 
			double [ ] [ ] A = new double[][] {{0}};
			double [ ]  tau = new double[] {0};
			ContPhaseVar ptRV=new DenseContPhaseVar(tau, A);

			network.getVertexByID(dataA.Source).pulse(0, ptRV, 1, 0,0, finalPath);
			//finalPath = completeThePath(dataA,network);
			finalPath = network.Path;

			//Ends the clock

			double finalTime = (System.nanoTime() - iniTime)/1000000000;

		//Print results

		System.out.println("The final path is: "+finalPath);
		System.out.println("The total cost is: "+network.PrimalBound);
		System.out.println("The time limit is: "+network.TimeC+"sec");
		System.out.println("The time limit is: "+network.TimeC/60+"min");
		System.out.println("The probability of arriving in time is: "+(network.finalProb));
		System.out.println("The computational time is: "+finalTime+" s");
		System.out.println(network.getVertexByID(dataA.Source).getMaxTime());
		System.out.println(network.getVertexByID(dataA.Source).getMinTime());
		//double timeLimit = (double) ((network.getVertexByID(dataA.Source).getMaxTime() - network.getVertexByID(dataA.Source).getMinTime())*k + network.getVertexByID(dataA.Source).getMinTime());
		//double timeLimit = (13015-9310)*k+9310;
		System.out.println("tl: "+timeLimit);
		 */
	}



	/**
	 * This method creates a network
	 * @param data
	 * @return the final graph
	 */
	/**
	private static PulseGraph createGraph(DataHandler data) {
		int numNodes = data.NumNodes;
		PulseGraph Gd = new PulseGraph(numNodes);
		for (int i = 0; i < numNodes; i++) {
			if(i!=(data.LastNode-1)){
				Gd.addVertex(new VertexPulse(i) ); //Primero lo creo, y luego lo meto. El id corresponde al n�mero del nodo
			}
		}
		FinalVertexPulse vv = new FinalVertexPulse(data.LastNode-1);
		Gd.addFinalVertex(vv);
		for(int i = 0; i <data.NumArcs; i ++){
			Gd.addWeightedEdge( Gd.getVertexByID(data.Arcs[i][0]), Gd.getVertexByID(data.Arcs[i][1]),data.Distance[i],data.Time[i], i);			
		}
		return Gd;
	}
	 */

	/**
	 * This method creates a Stochastic_network
	 * @param data fitter
	 * @return the final graph
	 */

	private static PulseGraph S_createGraph(Fitter  data) {
		int numNodes = data.NumNodes;
		PulseGraph Gd = new PulseGraph(numNodes);
		for (int i = 0; i < numNodes; i++) {
			if(i!=(data.LastNode)){

				Gd.addVertex(new VertexPulse(i) ); //Primero lo creo, y luego lo meto. El id corresponde al n�mero del nodo
			}
		}

		FinalVertexPulse vv = new FinalVertexPulse(data.LastNode);
		Gd.addFinalVertex(vv);


		for(int i = 0; i <data.NumArcs; i ++){

			Gd.addWeightedEdge( Gd.getVertexByID(data.Arcs[i][0]), Gd.getVertexByID(data.Arcs[i][1]),data.Distance[i],data.MinTime[i],data.TimeRV[i], i);
		}
		return Gd;
	}

	/**
	 * This method runs a shortest path for cost and times
	 * @param data
	 * @param network
	 * @throws InterruptedException
	 */
	/**
	private static void SP(DataHandler data, PulseGraph network) throws InterruptedException {
		// Create two threads and run parallel SP for the initialization		
		//Thread tTime = new Thread();
		Thread tDist = new Thread();
		// Reverse the network and run SP for distance and time 
		DukqstraDist spDist = new DukqstraDist(network, data.LastNode-1);
		//DukqstraTime spTime = new DukqstraTime(network, data.LastNode-1);
		tDist = new Thread(new ShortestPathTask(1, spDist, null));
		//tTime = new Thread(new ShortestPathTask(0, null,  spTime));
		tDist.start();
		//tTime.start();
		tDist.join();
		//tTime.join();
	}
	 */

	/**
	 * This method runs a shortest path for cost
	 * @param data
	 * @param network
	 * @throws InterruptedException
	 */
	private static void SP(Fitter data, PulseGraph network) throws InterruptedException {
		// Create two threads and run parallel SP for the initialization		
		Thread tTime = new Thread();
		Thread tDist = new Thread();

		// Reverse the network and run SP for distance and time 
		DukqstraDist spDist = new DukqstraDist(network, data.LastNode);
		DukqstraTime spTime = new DukqstraTime(network, data.LastNode);
		tDist = new Thread(new ShortestPathTask(1, spDist, null));

		tTime = new Thread(new ShortestPathTask(0, null,  spTime));
		tDist.start();

		tTime.start();
		tDist.join();

		tTime.join();
	}

	/**
	 * This method completes the path using the minimum cost or minimum time path
	 * @param data
	 * @param network
	 * @return the final path
	 */
	public static ArrayList<Integer> completeThePath(DataHandler data, PulseGraph network){
		ArrayList<Integer> path = new ArrayList<Integer>();
		for(int i = 0;i<network.Path.size();i++){
			path.add(network.Path.get(i));
		}


		int nodoInicial = network.getFinalNode();
		boolean termine = false;
		double costoAcumulado = network.getFinalCost();
		double tiempoAcumulado = network.getFinalTime();


		while(termine == false) {
			int nodoActual = data.LastNode;
			for(int i = 0; i < network.getVertexByID(nodoInicial).magicIndex.size(); i++) {
				int e = (Integer) network.getVertexByID(nodoInicial).magicIndex.get(i);
				int a = DataHandler.Arcs[e][1];

				if(network.best == 1){
					if(costoAcumulado + DataHandler.Distance[e] + network.getVertexByID(a).minDist == network.PrimalBound ) {
						if(tiempoAcumulado+ DataHandler.Time[e] + network.getVertexByID(a).maxTime == network.ProbStar) {
							costoAcumulado+=DataHandler.Distance[e];
							tiempoAcumulado+=DataHandler.Time[e];
							nodoActual = a;	
						}
					}
				}else{
					if(costoAcumulado + DataHandler.Distance[e] + network.getVertexByID(a).maxDist == network.PrimalBound ) {
						if(tiempoAcumulado+ DataHandler.Time[e] + network.getVertexByID(a).minTime == network.ProbStar) {
							costoAcumulado+=DataHandler.Distance[e];
							tiempoAcumulado+=DataHandler.Time[e];
							nodoActual = a;	
						}
					}
				}
			}

			path.add(nodoActual);
			if(nodoActual == data.LastNode) {
				termine = true;
			}else {
				nodoInicial = nodoActual;
			}
		}

		return path;
	}

	/**
	 * This method tells if the final path was found by a minimum time or a minimum cost path completion
	 * @param network
	 * @return 1: cost 2: time
	 */
	public static String whoFindThePath(PulseGraph network) {
		String rta = "There is no point";
		if(network.best == 1) {
			return "Forward cost path completion";
		}
		else if (network.best == 2) {
			return "Forward time path completion";
		}
		return rta;
	}
	public static boolean in(int i, int[] list) {
		boolean is =false;

		for (int j = 0; j < list.length; j++) {
			if (i==list[j]) {
				is=true;
				break;
			}
		}

		return is;

	}




}
