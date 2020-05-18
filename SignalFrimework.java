import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.log;
import static java.lang.Math.sin;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.lang.reflect.InvocationTargetException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.ItemLabelAnchor;
import org.jfree.chart.labels.ItemLabelPosition;
import org.jfree.chart.labels.StandardCategoryItemLabelGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.time.Millisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

import Jama.Matrix;
import Jama.QRDecomposition;

public class SignalFrimework {
	
	public SignalFrimework() {
		
	}
	/******* Drawing column statistical diagram 
	 * @throws InterruptedException 
	 * @throws InvocationTargetException *******/
	public void showDiagram(double[][] disqFrekw) throws InvocationTargetException, InterruptedException {
	     SwingUtilities.invokeAndWait(()->{
	            Data example=new Data("Bar Chart Window",disqFrekw);
	            example.setSize(800, 400);
	            example.setLocationRelativeTo(null);
	            example.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	            example.setVisible(true);
	          });
	        
	}
	  public CategoryDataset createDataset(double[][] arr) {
		    DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		    
		    for(int i =0; i < arr[0].length; i++) {
		    	dataset.addValue(arr[0][i]," Freq",String.valueOf(arr[1][i]));
		    }

		    return dataset;
		  }
	/******* Calculating Freaquencies and distinquish *******/
	public double[][] FreqDuisq(double[] buildPC){
		
		final double[] distinct = Arrays.stream(buildPC)
	            .distinct()
	            .toArray();
	        
	        double[] frekwnts = new double[distinct.length];
	        double frekw = 0.0;
	        for(int i =0; i < distinct.length; i ++) {
	        	 
	        	for(int j = 0; j < buildPC.length; j++ ) {
	        		if(distinct[i] == buildPC[j]) {
	        			frekw = frekw+1;
	        			frekwnts[i] = frekw;
	        		}
	        	}
	        	System.out.println(" Number: "+distinct[i]+" : "+" Frekw: "+frekwnts[i]);
	        	frekw = 0.0;
	        	
	        	
	        }
	        double[][] result = {distinct, frekwnts};
		return result;
	}
	/******* End *******/
	
	/******* Calculate FFT ******/
	public static int bitReverse(int n, int bits) {
		int reversedN = n;
		int count = bits - 1;

		n >>= 1;
		while (n > 0) {
			reversedN = (reversedN << 1) | (n & 1);
			count--;
			n >>= 1;
		}

		return ((reversedN << count) & ((1 << bits) - 1));
	}

	static void fft(Complex[] buffer) {

		int bits = (int) (log(buffer.length) / log(2));
		for (int j = 1; j < buffer.length / 2; j++) {
			int swapPos = bitReverse(j, bits);
			Complex temp = buffer[j];
			buffer[j] = buffer[swapPos];
			buffer[swapPos] = temp;
		}

		for (int N = 2; N <= buffer.length; N <<= 1) {
			for (int i = 0; i < buffer.length; i += N) {
				for (int k = 0; k < N / 2; k++) {
					try {
						int evenIndex = i + k;
						int oddIndex = i + k + (N / 2);
						Complex even = buffer[evenIndex];
						Complex odd = buffer[oddIndex];
						double term = (-2 * PI * k) / (double) N;
						Complex exp = (new Complex(cos(term), sin(term)).mult(odd));
						buffer[evenIndex] = even.add(exp);
						buffer[oddIndex] = even.sub(exp);
					} catch (Exception e) {
					}
				}
			}
		}
	}

	/****** END   *****/
	
	/***** Fill FFT arrays *****/
	
	public double[] fftArray(SignalFrimework signalFrimework, double[] y1, double[] fft, double[] fft1) {
		Complex[] cinput = new Complex[y1.length];
		for (int i = 0; i < y1.length - 1; i++)
			cinput[i] = new Complex(y1[i], 0.0);
		signalFrimework.fft(cinput);
		fft = new double[cinput.length];
		fft1 = new double[cinput.length];
		System.out.println("Results:");
		int counter = 0;
		for (Complex c : cinput) {
			counter++;
			System.out.println(c);
			try {
				String[] a = c.toString().split(",");
				fft[counter] = Double.parseDouble(a[0]);
				fft1[counter] = Double.parseDouble(a[1]);
			} catch (Exception e) {}
		}
		return fft;
	}

	/******* END *********/
	/****** Calculate Shanonn Entropy ******/ 
	 @SuppressWarnings("boxing")
	  public static double getShannonEntropy(String s) {
	    int n = 0;
	    Map<Character, Integer> occ = new HashMap<>();
	 
	    for (int c_ = 0; c_ < s.length(); ++c_) {
	      char cx = s.charAt(c_);
	      if (occ.containsKey(cx)) {
	        occ.put(cx, occ.get(cx) + 1);
	      } else {
	        occ.put(cx, 1);
	      }
	      ++n;
	    }
	 
	    double e = 0.0;
	    for (Map.Entry<Character, Integer> entry : occ.entrySet()) {
	      char cx = entry.getKey();
	      double p = (double) entry.getValue() / n;
	      e += p * log2(p);
	    }
	    return -e;
	  }
	 
	  private static double log2(double a) {
	    return Math.log(a) / Math.log(2);
	  }
	  
	  /******* Get Entropy Array ******/
	  public double[] getEntropyRes(double[] fft, SignalFrimework signalFrimework,
			  double[] entropyGraph) {
		  
		  for (int i = 0; i < fft.length; i++) {
				 
				double entropy = getShannonEntropy(String.valueOf(fft[i]));
				entropyGraph[i] = entropy;
				//System.out.printf("Shannon entropy of %40s: %.12f%n", "\"" + fft[i] + "\"", entropy);
			}
		  return entropyGraph;
	  
	  }
	  /****** END *****/
	  /**** Calculate Variation and Kovariation Koeficient Value ****/
	  public double[] Variation(double[] y) {
			 //double x_avg = 6.55; /* average of the data*/
			 double sum = 0;
			 for(int i = 0; i<y.length;i++) {
				sum +=y[i]; 
			 }
			 double x_avg = sum / y.length;
			 double[] l = new double [y.length]; //temporary data
			 double temp1=0;
			 double var=0;
			 double SD=0; //standard deviation
			 for (int h=0; h<y.length; h++){
			 l[h] = ((Double.parseDouble (String.valueOf(y[h]))-x_avg)*(y[h])-x_avg);
			 temp1 = temp1+l[h];
			 }
			 var = temp1/y.length;
			 SD = Math.pow((var),(0.5));
			 double[] m = new double [y.length]; //sum of y(h)
			 double[] p = new double [y.length]; //temporary data
			 double[] y11 = new double [y.length]; //ACF
			 for (int h=0; h<y.length; h++){
			 for (int t=1; t<y.length-h; t++){
			 p [h] = p[h]+(y[t+h]-x_avg)*(y[t-1]-x_avg);
			 m [h] = (1.0/y.length)*p[h];
			 }
			 y11[h] = m[h]/var;
			 }
			return y11;
	  }
	  /**** END ****/
	  /******* Draw Column Diagram-Graph **********/
	  /******* Creating a Dataset *********/
	  public  CategoryDataset createDataset(double[] array) {
	        String row = "Row";
	        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
	        for(int i =0; i < array.length; i =i+100) {
	        	dataset.addValue(array[i], row, String.valueOf(i));
	        }
	        return dataset;
	    }
	  /***** End *********/
	  /**** Plotting Chart ********/
	  public JFreeChart createChart(CategoryDataset dataset) {
	        CategoryAxis categoryAxis = new CategoryAxis("");
	        ValueAxis valueAxis = new NumberAxis("");
	        valueAxis.setVisible(false);
	        BarRenderer renderer = new BarRenderer() {

	            @Override
	            public Paint getItemPaint(int row, int column) {
	                switch (column) {
	                    case 0:
	                        return Color.red;
	                    case 1:
	                        return Color.yellow;
	                    case 2:
	                        return Color.blue;
	                    case 3:
	                        return Color.orange;
	                    case 4:
	                        return Color.gray;
	                    case 5:
	                        return Color.green.darker();
	                    default:
	                        return Color.red;
	                }
	            }
	        };
	        renderer.setDrawBarOutline(false);
	        renderer.setBaseItemLabelGenerator(new StandardCategoryItemLabelGenerator());
	        renderer.setBasePositiveItemLabelPosition(new ItemLabelPosition(
	            ItemLabelAnchor.OUTSIDE12, TextAnchor.BOTTOM_CENTER));
	        renderer.setBaseItemLabelsVisible(Boolean.TRUE);
	        renderer.setBarPainter(new StandardBarPainter());
	        CategoryPlot plot = new CategoryPlot(dataset, categoryAxis, valueAxis, renderer);
	        JFreeChart chart = new JFreeChart("", JFreeChart.DEFAULT_TITLE_FONT, plot, false);
	        chart.setBackgroundPaint(Color.white);
	        return chart;
	    }
	  /****** End *****/
	  /**** Display Diagram *********/
	  public void display(double[] array) {
	        JFrame f = new JFrame("BarChart");
	        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	        f.add(new ChartPanel(createChart(createDataset(array))));
	        f.pack();
	        f.setLocationRelativeTo(null);
	        f.setVisible(true);
	    }
	  /**** End ****/
	  
	  /******** Showing All ******/
	  
	  public void ShowAll(double[] kovariation) {
		  EventQueue.invokeLater(() -> {
	            new SignalFrimework().display(kovariation);
	        });
	  }
	  /************ End All ***/
	  
	  /***** Level-Shift Metric *****/
	  public static Double[] maxsInEveryWindows(double[] arr, int k) {
		    Deque<Integer> deque = new ArrayDeque<Integer>();
		    /* Process first k (or first window) elements of array */
		    for (int i = 0; i < k; i++) {
		        // For very element, the previous smaller elements are useless so
		        // remove them from deque
		        while (!deque.isEmpty() && arr[i] >= arr[deque.peekLast()]) {
		            deque.removeLast(); // Remove from rear
		        }
		        // Add new element at rear of queue
		        deque.addLast(i);
		    }
		    List<Double> result = new ArrayList<Double>();
		    // Process rest of the elements, i.e., from arr[k] to arr[n-1]
		    for (int i = k; i < arr.length; i++) {
		        // The element at the front of the queue is the largest element of
		        // previous window, so add to result.
		        result.add(arr[deque.getFirst()]);
		        // Remove all elements smaller than the currently
		        // being added element (remove useless elements)
		        while (!deque.isEmpty() && arr[i] >= arr[deque.peekLast()]) {
		            deque.removeLast();
		        }
		        // Remove the elements which are out of this window
		        while (!deque.isEmpty() && deque.getFirst() <= i - k) {
		            deque.removeFirst();
		        }
		        // Add current element at the rear of deque
		        deque.addLast(i);
		    }
		    // Print the maximum element of last window
		    result.add(arr[deque.getFirst()]);

		    return result.toArray(new Double[0]);
		}
	  
	  
	  /********** Level Shift Metric Test ***********/
	// finds the max and its index
	    static double[] findMaxByIteration(double arr[], int start, int end)
	    {
	        double max, max_ndx; 

	        max = arr[(int)start];
	        max_ndx = start;
	        for (int i=start; i<end; i++)
	        {
	            if (arr[i] > max)
	            {
	                max = arr[i];
	                max_ndx = i;
	            }    
	        }

	        double result[] = {max, max_ndx};

	        return result;
	    }

	    // optimized to skip iteration, when previous windows max element 
	    // is present in current window
	    static double[] optimizedPrintKMax(double arr[], int n, int k)
	    {
	        int i, j;
	        double[] lShiftMatric = new double[arr.length];
	        double max, max_ndx;

	        // for first window - find by iteration.    
	        double result[] = findMaxByIteration(arr, 0, k);

	        //System.out.println(" Level Shift Matric: "+ result[0]);

	        max = result[0];
	        max_ndx = result[1];   

	         for (j=1; j <= (n-k); j++)
	         {
	            // if previous max has fallen out of current window, iterate and find
	            if (max_ndx < j)  
	            {
	                result = findMaxByIteration(arr, j, j+k);
	                max = result[0];
	                max_ndx = result[1];   
	            } 
	            // optimized path, just compare max with new_elem that has come into the window 
	            else 
	            {
	                int new_elem_ndx = j + (k-1);
	                if (arr[new_elem_ndx] > max)
	                {
	                    max = arr[new_elem_ndx];
	                    max_ndx = new_elem_ndx;
	                }      
	            }
	            lShiftMatric[j] = max;
	            //System.out.println(" Level Shift Finall :  "+ max);
	         }
	         return lShiftMatric;
	    } 
	    
	    /******* Test Level Shift *******/
	    /***** End *****/
	    
	  /***** End *****/
	    /****** Computing Variance Change ******/
		static double variance(double a[],  
	            int n) 
	{ 
	// Compute mean (average  
	// of elements) 
	double sum = 0; 

	for (int i = 0; i < n; i++) 
	sum += a[i]; 
	double mean = (double)sum /  
	       (double)n; 

	// Compute sum squared  
	// differences with mean. 
	double sqDiff = 0; 
	for (int i = 0; i < n; i++)  
	sqDiff += (a[i] - mean) *  
	       (a[i] - mean); 

	return (double)sqDiff / n; 
	} 

	static double standardDeviation(double arr[],  
	                     int n) 
	{ 
	return Math.sqrt(variance(arr, n)); 
	} 
	/***** End *****/
	/***** Variance Algorithm 2 Test  Rolling Variance Algorithm *****/
	public double[] VarianceTest(double[] fft, double[] accVarRes, double mean,double accVar, int n) {
		//var queue = new Queue(fft.length);
	    Queue<Double> queue = new PriorityQueue<>();
	    for(int i =1; i < fft.length; i++) {
	    	queue.add(fft[i]);
	        //queue.Enqueue(fft[i]);
	        if (n < fft.length)
	        {
	            // Calculating first variance
	            n++;
	            double delta = fft[i] - mean;
	            mean += delta / n;
	            accVar += delta * (fft[i] - mean);
	        }
	        else
	        {
	            // Adjusting variance
	            //double then = queue.Dequeue();
	        	double then = queue.peek();
	            double prevMean = mean;
	            mean += (fft[i] - then) / fft.length;
	            accVar += (fft[i] - prevMean) * (fft[i] - mean) - (then - prevMean) * (then - mean);
	        }

	        if (n == fft.length)
	        	accVarRes[i] = accVar / (fft.length - 1);
	            /*yield return*/// System.out.println(accVar / (fft.length - 1));
	    }
	    return accVarRes;
	}
	/***** END *****/
	/***** Calculating Curvature *****/
	 public static double[] createAndShowGUI(double[][] PCA, double[] curveature)
	    {
		 CurvatureFromThreePointsPanel cpf = new CurvatureFromThreePointsPanel(PCA, curveature);
		 return cpf.CurvatureFromThreePoints(PCA,  curveature);
	    }
 
	 public void CurvtureCalc(double[][] PCA, double[] curveature) {
		 SwingUtilities.invokeLater(new Runnable()
	        {
	            @Override
	            public void run()
	            {
	                createAndShowGUI(PCA, curveature);
	                
	            }
	        });
	 }
	 
/***** END *****/
	  
	  
/************** Build Scarred PCA Graph ********************/
	public XYDataset createDataset(double[] arr1, double[] arr2) {
		    XYSeriesCollection dataset = new XYSeriesCollection();
		    XYSeries series1 = new XYSeries("PCA1");
		    for(int i =0; i < arr1.length; i++) {
		    series1.add(i, arr1[i]);
		    }

		    dataset.addSeries(series1);
		    
		    XYSeries series2 = new XYSeries("PCA2");
		    for(int i =0; i < arr2.length; i++) {
		    series2.add(i, arr2[i]);
		    }

		    dataset.addSeries(series2);

		    return dataset;
		  }
/************* END **********************/

/*********** Calculating Spikiness **********/
	/******* Calculating spikiness Test Version *******/
    // Implementing output data like double array
	public double[] calcSpikiness(double[] buildPC){
	    double[] median = new double[buildPC.length];
	    double[] finalSignal = new double[buildPC.length];
	    //Calculating mean
	    double w = 0.0;
	    for (int i = 0; i < buildPC.length-1; i++)
	    	for (int j = 0; j < buildPC.length; j++) {
	    		if(buildPC[i] < buildPC[j]) {
	    			w = buildPC[i];
	    			median[i] = median[j];
	    			median[j] = w;
	    			//System.out.println(" Median: "+ w);
	    		}
	    		
	    	}
	    /***** Final Signa Calculating *****/
	    for(int i = 0; i <median.length; i++) {
	      finalSignal[i] = buildPC[i] - median[i];
	      //System.out.println(" finalSignal[i]: "+finalSignal[i]);
	    }
	   //Skipes = difference between basical signal and median signal
	    // For example h[n] = |f[n]- g[n]|
	    return finalSignal;
	}
/********** End ***********************/
/****** Caqlculating fspots step1 *******/
	/* divide cuting intervals */
	

	static double[][] countIntervals(double[][] arr, int V, int N, int multiplicateVariable) 
    { 
        // Variable to store the count of intervals 
        int count = 0; 
        double[][] ranges = new double[arr.length][arr[0].length];
        // Variables to store start and end of an interval 
        double li, ri; 
  
        for (int i = 0; i < arr.length; i++) { 
        	for(int j =0; j <arr[0].length; j++) {
            li = arr[i][0]*multiplicateVariable; 
            ri = arr[i][j]*multiplicateVariable; 
  
            // Implies V lies in the interval 
            // so increase count 
            if (V >= li && V <= ri) {
            	ranges[0][j] = li;
                ranges[1][j] = ri;
                count++; 
            }
        } }
        return ranges; 
    } 
	
	/***** Step2 Calculating RLE ******/
	  public double[] calculatingRLE(double[] arr) 
	    { 
		    double[] rle = new double[arr.length];
	        int n = arr.length; 
	        for (int i = 0; i < n; i++) { 
	  
	            // Count occurrences of current character
	        	int count = 0;
	        	if(i < n - 1 && (arr[i] == arr[i+1])) {
	        	  count ++;
	        	}
	            rle[i] = arr[i];
	        } 
	        
	        return arr;
	    } 
	  /****** Calculating Flat Spots *****/
	  public double[] calculatingFlatspots(double[][] arr) {
	    	 double[] flatSpots = new double[arr[0].length];
	    	 for(int i = 0; i <arr[0].length; i++) {
	    		 if(arr[0][i]>arr[1][i]) {
	    		    flatSpots[i] = arr[0][i];
	    		 }
	    		 else {flatSpots[i] = arr[1][i];}
	    	 }
	    	 /*for(int i =0; i< flatSpots.length;i++) {
	    		 if(flatSpots[i] ==0) {
	    			flatSpots[i] =1; 
	    		 }
	    	 }*/
	    	 return flatSpots;
	     
	     }
	  /****** End ******/
	  
	  
	  /********* Drawing Time Diagram ***************/
	  
	  public JFreeChart createChart(XYDataset dataset) {
	        JFreeChart chart = ChartFactory.createTimeSeriesChart(
	            "Sensor Data Values", "Time", "Value", dataset, true, true, false);
	        XYPlot plot = (XYPlot) chart.getPlot();
	        plot.setDomainCrosshairVisible(true);
	        plot.setRangeCrosshairVisible(true);
	        XYItemRenderer r = plot.getRenderer();
	        if (r instanceof XYLineAndShapeRenderer) {
	            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) r;
	            renderer.setBaseShapesVisible(true);
	            renderer.setBaseShapesFilled(true);
	        }
	        DateAxis axis = (DateAxis) plot.getDomainAxis();
	        axis.setDateFormatOverride(new SimpleDateFormat("dd-MMM-yyyy HH:mm"));
	        return chart;
	    }

	    public static XYDataset createDataset(double[] arr,Connection conn, String sql) {
	    	
	        TimeSeries series = new TimeSeries("Value");
	        DateFormat simpleDateFormat=new SimpleDateFormat("yyyy-MM-dd");
	        try {
	            Statement st = conn.createStatement();
	            st = conn.createStatement();
	            ResultSet rs = st.executeQuery(sql);
	            int counter = 0;
	            while (rs.next()) {
	            	if(counter == arr.length-1) {
	            		break;
	            	}
	            	counter++;
	            	Timestamp  time = rs.getTimestamp("msr_txdt");
	                series.add(new Millisecond(time), arr[counter]);
	           }
	                
	        } catch (Exception e) {
	            e.printStackTrace(System.err);
	        }
	        TimeSeriesCollection dataset = new TimeSeriesCollection(series);
	        return dataset;
	    }
		
	    public void ShowTimeDiagram(double[] clustersR, Connection conn, String sql) {
	    	 SwingUtilities.invokeLater(new Runnable() {
		            @Override
		            public void run() {
		                /*Data demo = new Data(buildPC[1],conn,sql);*/
		            	Data demo = new Data(clustersR,conn,sql);
		                demo.setDefaultCloseOperation( javax.swing.WindowConstants.EXIT_ON_CLOSE );
		                demo.pack();
		                demo.setLocationRelativeTo(null);
		                demo.setVisible(true);
		            }
		        });
	    }
	    
	    public void ShowScarredDoubleDiagram(double[][] buildPC) {
	    	SwingUtilities.invokeLater(() -> {
		        Data example = new Data(buildPC[0],buildPC[1]);
		        example.setSize(800, 400);
		        example.setLocationRelativeTo(null);
		        example.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		        example.setVisible(true);
		      });
	    }
	  /********* End ***************************/
	  

/******* !!!!!!!!!! End Class *******************/	
	
}


class CurvatureFromThreePointsPanel extends JPanel 
implements MouseListener, MouseMotionListener
{
private List<Point2D> pointList;
private Point2D draggedPoint;

public CurvatureFromThreePointsPanel(double[][] PCA, double[] curveature)
{
    this.pointList = new ArrayList<Point2D>();
    for(int i = 0; i < PCA[0].length; i++ ) {
      pointList.add(new Point2D.Double(PCA[0][i],PCA[1][i]));
    }
    Point2D p0 = pointList.get(0);
    Point2D p1 = pointList.get(1);
    Point2D p2 = pointList.get(2); 
    for(int i = 0;  i < pointList.size(); i++) {
    	if(i+2 == pointList.size()-1) {
    		break;
    	}
    	curveature[i] = computeCurvature(pointList.get(i), pointList.get(i+1), pointList.get(i+2));
    }

    addMouseListener(this);
    addMouseMotionListener(this);
}

/***** Calculating Curveatures **********/
public double[] CurvatureFromThreePoints(double[][] PCA, double[] curveature)
{
    this.pointList = new ArrayList<Point2D>();
    for(int i = 0; i < PCA[0].length; i++ ) {
      pointList.add(new Point2D.Double(PCA[0][i],PCA[1][i]));
    }
    Point2D p0 = pointList.get(0);
    Point2D p1 = pointList.get(1);
    Point2D p2 = pointList.get(2); 
    for(int i = 0;  i < pointList.size(); i++) {
    	if(i+2 == pointList.size()-1) {
    		break;
    	}
    	curveature[i] = computeCurvature(pointList.get(i), pointList.get(i+1), pointList.get(i+2));
    }

    return curveature;
}

private static double computeCurvature(Point2D p0, Point2D p1, Point2D p2)
{
    double dx1 = p1.getX() - p0.getX();
    double dy1 = p1.getY() - p0.getY();
    double dx2 = p2.getX() - p0.getX();
    double dy2 = p2.getY() - p0.getY();
    double area = dx1 * dy2 - dy1 * dx2;
    double len0 = p0.distance(p1);
    double len1 = p1.distance(p2);
    double len2 = p2.distance(p0);
    return 4 * area / (len0 * len1 * len2);
}


// Adapted from https://stackoverflow.com/a/4103418
private static Point2D computeCircleCenter(
    Point2D p0, Point2D p1, Point2D p2)
{
    double x0 = p0.getX();
    double y0 = p0.getY();
    double x1 = p1.getX();
    double y1 = p1.getY();
    double x2 = p2.getX();
    double y2 = p2.getY();
    double offset = x1 * x1 + y1 * y1;
    double bc = (x0 * x0 + y0 * y0 - offset) / 2.0;
    double cd = (offset - x2 * x2 - y2 * y2) / 2.0;
    double det = (x0 - x1) * (y1 - y2) - (x1 - x2) * (y0 - y1);
    double invDet = 1 / det;
    double cx = (bc * (y1 - y2) - cd * (y0 - y1)) * invDet;
    double cy = (cd * (x0 - x1) - bc * (x1 - x2)) * invDet;
    return new Point2D.Double(cx, cy);
}

@Override
protected void paintComponent(Graphics gr)
{
    super.paintComponent(gr);
    Graphics2D g = (Graphics2D)gr;

    g.setColor(Color.RED);
    for (Point2D p : pointList)
    {
        double r = 5;
        g.draw(new Ellipse2D.Double(p.getX()-r, p.getY()-r, r+r, r+r));
    }

    g.setColor(Color.BLACK);

    Point2D p0 = pointList.get(0);
    Point2D p1 = pointList.get(1);
    Point2D p2 = pointList.get(2);
    double curvature = computeCurvature(p0, p1, p2);
    g.drawString("Curvature: "+curvature, 10,  20);

    Point2D center = computeCircleCenter(p0, p1, p2);
    double radius = center.distance(p0);
    g.draw(new Ellipse2D.Double(
        center.getX() - radius, center.getY() - radius,
        radius + radius, radius + radius));
}

@Override
public void mouseDragged(MouseEvent e)
{
    if (draggedPoint != null)
    {
        draggedPoint.setLocation(e.getX(), e.getY());
        repaint();

        System.out.println("Points: ");
        for (Point2D p : pointList)
        {
            System.out.println("    "+p);
        }
    }
}



@Override
public void mousePressed(MouseEvent e)
{
    final double thresholdSquared = 10 * 10;
    Point2D p = e.getPoint();
    Point2D closestPoint = null;
    double minDistanceSquared = Double.MAX_VALUE;
    for (Point2D point : pointList)
    {
        double dd = point.distanceSq(p);
        if (dd < thresholdSquared && dd < minDistanceSquared)
        {
            minDistanceSquared = dd;
            closestPoint = point;
        }
    }
    draggedPoint = closestPoint;
}

@Override
public void mouseReleased(MouseEvent e)
{
    draggedPoint = null;
}

@Override
public void mouseMoved(MouseEvent e)
{
    // Nothing to do here
}


@Override
public void mouseClicked(MouseEvent e)
{
    // Nothing to do here
}

@Override
public void mouseEntered(MouseEvent e)
{
    // Nothing to do here
}


@Override
public void mouseExited(MouseEvent e)
{
    // Nothing to do here
}


}

