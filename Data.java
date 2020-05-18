//import static org.junit.Assert.assertEquals;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
//import de.lmu.ifi.dbs.elki.math.MeanVariance;
import java.awt.Color;
import java.lang.Object;
import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.CategoryDataset;


//import org.junit.Test;

class Data extends JFrame {
	/******* Variables *******/
	static SignalFrimework signalFrimework = new SignalFrimework();
	static String[] parseL44;
	static String[] parseL4;
	static Clustering calculateClustering = new Clustering();
	static Functions functions1111 = new Functions();
	static Functions functions2222 = new Functions();
	static Functions functions11 = new Functions();
	static Functions functions12 = new Functions();
	static Functions functions12223 = new Functions();
	static JSONObject obj = null;
	static JSONArray mToolData = null;
	JSONParser parser = new JSONParser();
	static double[] arrayEigenValue = null;
	static double[][] eigneVectors = null;
	static Functions functions = new Functions();
	Functions functions1 = new Functions();
	static CalcPart1 calc1 = new CalcPart1();
	static CalculateMetrix metrix = new CalculateMetrix();
	//static String timestampPoint1 = null;
	//static String timestampPoint2 = null;
	static List<List<Object>> dynamicND = new ArrayList<List<Object>>();
	static List<List<Object>> dynamicND1 = new ArrayList<List<Object>>();
	
	private static long time_start = 0;
	private static long time_stop = 0;
	private static long time2_start = 0;
	private static long time2_stop = 0;
    /******** End ********/
	
	
	/******* Constructors **********/
	  public Data(double[] arr1, double[] arr2) {
	    XYDataset dataset = signalFrimework.createDataset(arr1, arr2);
	    // Create chart
	    JFreeChart chart = ChartFactory.createScatterPlot(
	        "PCA MAtrix", 
	        "PCA1-Axis", "PCA2-Axis", dataset);
	    //Changes background color
	    XYPlot plot = (XYPlot)chart.getPlot();
	    plot.setBackgroundPaint(new Color(255,228,196));
	    // Create Panel
	    ChartPanel panel = new ChartPanel(chart);
	    setContentPane(panel);
	  }
	 
	  public Data(double[] arr, Connection conn, String sql) {
		// TODO Auto-generated constructor stub
		//super(title);
	        JFreeChart chart = signalFrimework.createChart(signalFrimework.createDataset(arr, conn,sql));
	        ChartPanel chartPanel = new ChartPanel(chart) {
	            @Override
	            public Dimension getPreferredSize() {
	                return new Dimension(650, 400);
	            }
	        };
	        chartPanel.setMouseZoomable(true, false);
	        add(chartPanel);
	}
	  

	  public Data(String appTitle, double[][] arr) {
	    super(appTitle);
	    // Create Dataset
	    CategoryDataset dataset = signalFrimework.createDataset( arr);
	    //Create chart
	    JFreeChart chart=ChartFactory.createBarChart(
	        "Principal Components", //Chart Title
	        "PC", // Category axis
	        "Freq", // Value axis
	        dataset,
	        PlotOrientation.VERTICAL,
	        true,true,false
	       );
	    ChartPanel panel=new ChartPanel(chart);
	    setContentPane(panel);
	  }
	  /********** End ***********/

	/***** Execute Function 
	 * @throws IOException 
	 * @throws java.text.ParseException *****/
	public static  void execute() throws SQLException, JSONException, ParseException, IOException, java.text.ParseException {
		
		//Connection conn = DriverManager.getConnection("jdbc:postgresql://localhost:5432/opcuadb1", "postgres",
				//"iaremko89q");
		 Connection conn =
		 DriverManager.getConnection("jdbc:postgresql://vm28.os.zps:5432/opcuadb",
		 "dbguest", "xxxxxx");
		Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
		String machineToolNumber = "2";
		//String timestampPoint1 = "'2019-12-11 10:40:07.914'";
   		//String timestampPoint2 = "'2019-12-11 11:10:47.809'";
		String timestampPoint1 = "'2019-12-11 09:42:20.638'";
   		String timestampPoint2 = "'2019-12-11 09:55:10.829'";
		//String timestampPoint1 = "'2019-12-11 10:02:25.494'";
   		//String timestampPoint2 = "'2019-12-11 10:33:22.557'";
		
		//String timestampPoint1 = "'2019-12-11 10:40:07.914'";
   		//String timestampPoint2 = "'2019-12-11 11:10:47.809'";
		//String timestampPoint1 = "'2019-12-11 10:45:00.279'";
   		//String timestampPoint2 = "'2019-12-11 10:50:00.073'";
   		
   		/*String timestampPoint1 = "'2019-10-08 10:50:10.939'";
   		String timestampPoint2 = "'2019-10-08 11:20:05.414'";*/
        List<Object> timestampDB = new ArrayList<Object>();
        List<Object> timestampDB22 = new ArrayList<Object>();
        List<Object> nodeObjects = new ArrayList<Object>();
        List<Object> nodeObjects1 = new ArrayList<Object>();
        String sql1 = "select msr_json, msr_txdt from opcmsr where" + " dev_id0 ="+machineToolNumber ;
        ResultSet rs = st.executeQuery(sql1);
        parseL44 = functions12223.takeDataFromDB(rs, parseL44, dynamicND1,nodeObjects1);
        functions12223.fillDataFromDB(rs, parseL44, dynamicND1, obj, mToolData,timestampDB22);
		/*String timestampPoint1 = "'2020-01-08 08:32:50.807'";
		String timestampPoint2 = "'2020-01-08 09:22:55.701'";*/
		
		/*String timestampPoint1 = "'2019-12-11 09:41:15.166'";
		String timestampPoint2 = "'2019-12-11 10:00:29.747'";*/
        /* Comparing Timestamps Values Here */
        /* Calculating comparation *****/
        
       
      		
        /*String timestampPoint1 = "'2019-12-11 10:43:43.605'";
		String timestampPoint2 = "'2019-12-11 11:13:30.436'";*/
		
		/*String timestampPoint1 = "'2019-12-11 10:43:43.605'";
		String timestampPoint2 = "'2019-12-11 11:13:30.436'";*/
		
		/*SimpleDateFormat format = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss",Locale.US);
		 Date d1 = format.parse(timestampPoint1.toString());
	     Date   d2 = format.parse(timestampPoint2);
         //in milliseconds
	     long diff = d2.getMinutes() - d1.getMinutes();
	     System.out.println(" Timestamp Diff: "+diff);*/
        /*System.out.println(" We are here "+timestampDB.size());
        for(int i = 0; i <timestampDB22.size(); i++) {
        	System.out.println(" Timestamp [i ]: "+i+" : "+timestampDB22.get(i));
        }
		//for (int i = 0; i< timestampDB.size(); i ++) {
			System.out.println(" Select Time Data Period from: "+0+" to "+timestampDB22.size());
			System.out.println(" Select Time Stamp from: "+timestampDB22.get(0)+" to "+timestampDB22.get(timestampDB22.size()-1));
			BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
			  //System.out.println(" Please write niumber from 1 -10: ");
			  try {
			  int number = Integer.parseInt(bufferRead.readLine().toString());
			  int index = number+3000;
			  System.out.println("Index: "+index);
			  System.out.println(timestampDB22.get(number));
			  timestampPoint1 = "'"+timestampDB22.get(number).toString()+"'";
			  timestampPoint2 = "'"+timestampDB22.get(index).toString()+"'";
				
			  }
		  catch(NumberFormatException e) {}
		//}
			  System.out.println(" Timesttamp Point1: "+timestampPoint1);
			  System.out.println(" Timesttamp Point2: "+timestampPoint2);*/
		
		String sql = "select msr_json, msr_txdt from opcmsr where msr_txdt between "+timestampPoint1+"and"+timestampPoint2+ "AND dev_id0 ="+machineToolNumber ;
		 rs = st.executeQuery(sql);
		// System.out.println(rs.getString(1));
		System.out.println(" Start: ");
		
		/* Dynamically Create N-Dimensional Matrix */
		parseL4 = functions.takeDataFromDB(rs, parseL4, dynamicND,nodeObjects);
		System.out.println(" Dynamic 2D size: " + dynamicND.size());
		functions.fillDataFromDB(rs, parseL4, dynamicND, obj, mToolData,timestampDB);
		//System.out.println(" Dynamic 2D size: " + dynamicND.size());
		
		try {
			//System.out.println(" Timestamps: "+timestampDB.size());
			//System.out.println(" nodeObjects: "+nodeObjects.size());
			/*for(int i =0; i < nodeObjects.size(); i++) {
				System.out.println(" nodeObjects: "+nodeObjects.get(i));
			}*/
            
			//double[][] data1 = functions.createMatrixArray(dynamicND);
			time_start = System.currentTimeMillis();
			double[][] buildPC = calc1.buildPComponents(functions, dynamicND, arrayEigenValue,eigneVectors,
					   signalFrimework,  conn,  sql, nodeObjects);
			time_stop = System.currentTimeMillis();
	    	System.out.println(" Cas celk...:[ms]" + (getTimePeriod()));
	    	//System.in.read();
			/**** End ****/
           
			double[][] clusteringPC = calc1.PrincipalComponentsClustering(buildPC, signalFrimework, metrix,calculateClustering,conn,  sql);
	        

			
		} catch (Exception e) {
			System.out.println(e);
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

class Complex {
	public final double re;
	public final double im;

	public Complex() {
		this(0, 0);
	}

	public Complex(double r, double i) {
		re = r;
		im = i;
	}

	public Complex add(Complex b) {
		return new Complex(this.re + b.re, this.im + b.im);
	}

	public Complex sub(Complex b) {
		return new Complex(this.re - b.re, this.im - b.im);
	}

	public Complex mult(Complex b) {
		return new Complex(this.re * b.re - this.im * b.im, this.re * b.im + this.im * b.re);
	}

	@Override
	public String toString() {
		String result = re+" , " + im;
		return result;
	}
}




