import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.InvocationTargetException;
import java.sql.Connection;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

public class CalcPart1 {
	
	private static long time_start = 0;
	private static long time_stop = 0;
	private static long time2_start = 0;
	private static long time2_stop = 0;
	
  public CalcPart1() {}
  
  
  /***** Principal Component Analyze Calculating 
 * @throws IOException 
 * @throws NumberFormatException ******/
  public double[][] buildPComponents(Functions functions, List<List<Object>> dynamicND,
		  double[] arrayEigenValue,double[][] eigneVectors,
		  SignalFrimework signalFrimework, Connection conn, String sql,List<Object> nodeObjects) throws NumberFormatException, IOException{
	  
	  Functions functions11 = new Functions(); 
	  Functions functions1 = new Functions(); 
	  
	  double[][] data1 = functions.createMatrixArray(dynamicND);
	  DrawData(data1, nodeObjects);
	  //double[] signal1 = functions11.DrawSinDiagram(data1, 38, 2);
	  //functions11.go(signal1);
	
      /**************** Calculating PCA *******************/
	  time_start = System.currentTimeMillis();
	  
		arrayEigenValue = functions.CovarianceMatrix(data1, arrayEigenValue);
		eigneVectors = functions.ReturnEigneVector(data1, eigneVectors, arrayEigenValue);
		double[][] buildPC = functions.buildPrincipalComponents1(2, eigneVectors,arrayEigenValue);
        
		
		time_stop = System.currentTimeMillis();
    	System.out.println(" Cas celk...:[ms]" + (getTimePeriod()));
    	//System.in.read();
    	
		int n = buildPC[0].length;
		double[] y = new double[n + 1];
		double[] y1 = new double[n + 1];
		/********** Build Realtime Graph ***********/

		// the function y = sin(4x) + sin(20x), sampled at n+1 points
		// between x = 0 and x = pi
		y = functions.DrawSinDiagram(buildPC, 0,1);
		y1 = functions1.DrawSinDiagram(buildPC, 1,1);
		double[] yTest = functions.DrawSinDiagram(buildPC, 0,300);
		double[] y1Test = functions1.DrawSinDiagram(buildPC, 1,300);
		double[] yTest1 = functions.DrawSinDiagram(buildPC, 0,30000);
		double[] y1Test1 = functions1.DrawSinDiagram(buildPC, 1,30000);
		double[][] frekwSign = {yTest1,y1Test1};
		signalFrimework.ShowScarredDoubleDiagram( buildPC);
		functions.go(yTest);

		functions1.go(y1Test);
		
		signalFrimework.ShowTimeDiagram(y,  conn, sql);
		signalFrimework.ShowTimeDiagram(y1,  conn, sql);
		
		return buildPC;
  }
  
  /****** END 
 * @throws ParseException 
 * @throws InterruptedException 
 * @throws InvocationTargetException *******/
  public double[][] PrincipalComponentsClustering( double[][] buildPC, SignalFrimework signalFrimework,
		  CalculateMetrix metrix, Clustering calculateClustering, Connection conn, String sql) throws ParseException, InvocationTargetException, InterruptedException{

	//double[][] metr = {y,y1};
	  Functions functions = new Functions();
	  Functions functions1 = new Functions();
	  Functions functions1111 = new Functions();
	  Functions functions2222 = new Functions();
		double[][] metrixFinal = metrix.getAllMatrix(buildPC, signalFrimework);
		double[] arrayEigenValueFinal = null;
		double[][] eigneVectorsFinal = null;
		Object[][] res = functions.calculatingEigneValAndVect(metrixFinal);
		arrayEigenValueFinal = functions.CovarianceMatrix(metrixFinal, arrayEigenValueFinal);
		eigneVectorsFinal = functions.ReturnEigneVector(metrixFinal, eigneVectorsFinal, arrayEigenValueFinal);
		//double[] arrayEigenValueFinal = new double[buildPC[0].length];
		//double[][] eigneVectorsFinal = new double[buildPC.length][buildPC[0].length];
		
		/*for(int i =0; i < res[0].length; i++) {
			arrayEigenValueFinal[i] = Double.parseDouble(res[2][i].toString());
			 eigneVectorsFinal[0][i]  = Double.parseDouble(res[0][i].toString());
			 eigneVectorsFinal[1][i]  = Double.parseDouble(res[1][i].toString());
		}*/
		
		
		double[][] calculatingDimensionalPC = functions1111.buildPrincipalComponents1(2, eigneVectorsFinal,arrayEigenValueFinal);
		
		
		double[] y11 = new double[buildPC[0].length + 1];
		double[] y22 = new double[buildPC[0].length + 1];
		y11 = functions.DrawSinDiagram(calculatingDimensionalPC, 0,1);
		y22 = functions1.DrawSinDiagram(calculatingDimensionalPC, 1,1);
		
		signalFrimework.ShowScarredDoubleDiagram( calculatingDimensionalPC);
		
		int minPoints = 10;
		double eps = 6;
		
		List<Clustering.Point> points = new ArrayList<Clustering.Point>();
		System.out.println("Original Points: ");
		for ( int j = 0; j < calculatingDimensionalPC[0].length; ++j) {
			Clustering.Point p = calculateClustering.new Point(calculatingDimensionalPC[0][j],calculatingDimensionalPC[1][j]);
			points.add(p);
		}
		//Collections.sort(points);
		List<Clustering.Cluster> clusters =calculateClustering.getClusters(points, minPoints, eps);
		Object[] finallCluster = null;
		for ( Clustering.Cluster c : clusters ) {
			finallCluster = c.points.toArray();
		}
		
		double[][] clustersR = calculateClustering.CorrectClusterRes(finallCluster);
		
		signalFrimework.ShowScarredDoubleDiagram( clustersR);
	
		functions1111.go(clustersR[0]);

		functions2222.go(clustersR[1]);
		
		signalFrimework.ShowTimeDiagram(clustersR[0],  conn, sql);
		signalFrimework.ShowTimeDiagram(clustersR[1],  conn, sql);
		
        double[][] disqFrekw = signalFrimework.FreqDuisq( clustersR[0]);
        signalFrimework.showDiagram(disqFrekw);
        
        return clustersR;
  }
  
  public void DrawData(double[][] data1,List<Object> nodeObjects) throws NumberFormatException, IOException {
	  
	  System.out.println(" Jaky merici hodnoty ze stroje budeme sledovat: ");
	  while(true) {
		  System.out.println(" Vyberte hodnotu ze seznamu: ");
		  for(int i = 0; i <nodeObjects.size(); i++) {
			  System.out.println(" Node Number: " +i+" "+nodeObjects.get(i));
		  }
		  BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
		  //System.out.println(" Please write niumber from 1 -10: ");
		  try {
		  int number = Integer.parseInt(bufferRead.readLine().toString());
		  Functions func = new Functions();
			double[] signal1 = func.DrawSinDiagram(data1, number, 2);
			func.go(signal1);
			if(number == 0) {
				  break;
			  }
		  }
	  catch(NumberFormatException e) {}
	  
  }
	  
	  

}
  
	public static int getTimePeriod() {
		return (int) (time_stop - time_start);
	}
	
	/**
	 * @return
	 */
	public int getTimePeriodMeasure() {
		return (int) (time2_stop - time2_start);
	}
	
}
