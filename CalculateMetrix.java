import java.util.Arrays;

public class CalculateMetrix {
	
	Functions functions00 = new Functions();
	Functions functions0 = new Functions();
	Functions functions = new Functions();
	Functions functions1 = new Functions();
	Functions functions2 = new Functions();
	Functions functions3 = new Functions();
	Functions functions4 = new Functions();
	Functions functions5 = new Functions();
	Functions functions6 = new Functions();
	Functions functions7 = new Functions();
	Functions functions8 = new Functions();
	Functions functions9 = new Functions();
	Functions functions10 = new Functions();
	Functions functions11 = new Functions();
	Functions functions12 = new Functions();
	
	public double[][] getAllMatrix(double[][] buildPC, SignalFrimework signalFrimework) {
		
		
		/********* Calculating Trand **************/
		
		Forecast forecast = new Forecast();
		double[] trand1 = forecast.fitTrand(buildPC[0], 0.8,0.2);
		double[] trand2 = forecast.fitTrand(buildPC[1], 0.8,0.2);
		functions00.go(trand1);
		/********** End ******************/
		/********* Calculating FFT ****************/
		double [] fft = null;
        double [] fft1 = null;
        double [] fft2 = null;
        double [] fft3 = null;
        fft = signalFrimework.fftArray( signalFrimework, buildPC[0], fft, fft1);
        fft1 = signalFrimework.fftArray( signalFrimework, buildPC[0], fft, fft1);
        fft2 = signalFrimework.fftArray( signalFrimework, buildPC[1], fft2, fft3);
        fft3 = signalFrimework.fftArray( signalFrimework, buildPC[1], fft2, fft3);
		/********* END ****************/
		/********* Calculating Entropy *******************/
		double[][] entropyGraph = new double[buildPC.length][buildPC[0].length];
		entropyGraph[0] = signalFrimework.getEntropyRes(buildPC[0],  signalFrimework, entropyGraph[0]);
		entropyGraph[1] = signalFrimework.getEntropyRes(buildPC[0],  signalFrimework, entropyGraph[1]);
		functions6.go(entropyGraph[1]);
		/******** End *************/
		/******* Kovariation ***********/
		// TODO code application logic here
		double[][] kovariation = new double[buildPC.length][buildPC[0].length];
		kovariation[0] =    signalFrimework.Variation(buildPC[0]);
		kovariation[1] =    signalFrimework.Variation(buildPC[1]);
		signalFrimework.display(kovariation[1]);
		
		/***** Level Shift Matric *****/
		double[][] lShiftMatric = new double[buildPC.length][buildPC[0].length];
        lShiftMatric[0] = signalFrimework.optimizedPrintKMax(buildPC[0], buildPC[0].length, 3);
        lShiftMatric[1] = signalFrimework.optimizedPrintKMax(buildPC[1], buildPC[1].length, 3);
        functions8.go(lShiftMatric[1]);
		
        /***** Variance Change *****/
	    double[][] accVarRes = new double[buildPC.length][buildPC[0].length];
	    double mean = 0;
	    double accVar = 0;
	    accVarRes[0] = signalFrimework.VarianceTest( buildPC[0],  accVarRes[0],  mean, accVar, buildPC[0].length);
	    accVarRes[1] = signalFrimework.VarianceTest( buildPC[1],  accVarRes[1],  mean, accVar, buildPC[1].length);
	    functions9.go(accVarRes[0]);
        /***** End *****/
	    /***** Calculating Curvature *****/
	    /***** Returns just one Array *****/
	    double[] curveatures1 = new double[buildPC[0].length];
	    double[][] curveatures22 = {buildPC[0],buildPC[1]};
	    curveatures1 = signalFrimework.createAndShowGUI(curveatures22, curveatures1);
	    functions10.go(curveatures1);
	    /***** END *****/
	    /****** Calculating Spikiness ******/
	    double[][] spikiness = new double[buildPC.length][buildPC[0].length];
	    spikiness[0] = signalFrimework.calcSpikiness(buildPC[0]);
	    spikiness[1] = signalFrimework.calcSpikiness(buildPC[1]);
	    signalFrimework.display(spikiness[1]);
	     int V = 10; 
	     // length of the array 
	     int N = buildPC.length;
	     double[][]cuttings = signalFrimework.countIntervals(buildPC, V, N, 1000);
	     double[] rle1 = signalFrimework.calculatingRLE(cuttings[0]);
	     double[] rle2 = signalFrimework.calculatingRLE(cuttings[1]);
	     /***** Return just on dimensional double array *****/
	     double[][] finalRle = {rle1, rle2};
	     double[] flatSpots = signalFrimework.calculatingFlatspots(finalRle);
	     signalFrimework.display(flatSpots);
	     
	     double[][] finalMetricsDimension = {trand1,trand2,fft,fft1,fft2,fft3,entropyGraph[0],entropyGraph[1],kovariation[0],kovariation[1],
	            		   lShiftMatric[0],lShiftMatric[1],accVarRes[0],accVarRes[1],
	            		   /*curveatures,*/curveatures1,/*getPolynomCurveatures,*/spikiness[0],spikiness[1],flatSpots};
	               
	    return finalMetricsDimension;
		
		
	}

}
