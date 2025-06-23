package aligner;

import java.util.Arrays;

import shared.Shared;

/**
 *Aligns two sequences to return ANI.
 *Uses only 2 arrays and avoids traceback.
 *Gives an exact answer.
 *Calculates rstart and rstop without traceback.
 *Limited to length 2Mbp with 21 position bits.
 *Restricts alignment to a fixed band around the diagonal.
 *
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date April 24, 2025
 */
public class BandedByteAligner implements IDAligner{

	/** Main() passes the args and class to Test to avoid redundant code */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	public BandedByteAligner() {}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final String name() {return "BandedByte";}
	@Override
	public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {return alignStatic(a, b, pos, rStart, rStop);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Tests for high-identity indel-free alignments needing low bandwidth */
	private static int decideBandwidth(byte[] query, byte[] ref) {
		int bandwidth=Math.min(100, 4+Math.max(query.length, ref.length)/8);
		int subs=0;
		for(int i=0, minlen=Math.min(query.length, ref.length); i<minlen && subs<bandwidth; i++) {
			subs+=(query[i]!=ref[i] ? 1 : 0);
		}
		return Math.min(subs+1, bandwidth);
	}

	/**
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @return Identity (0.0-1.0).
	 */
	public static final float alignStatic(byte[] query0, byte[] ref0, int[] posVector) {
		// Swap to ensure query is not longer than ref
		if(posVector==null && query0.length>ref0.length) {
			byte[] temp=query0;
			query0=ref0;
			ref0=temp;
		}

		final byte[] query=Factory.encodeByte(query0, (byte)15);
		final byte[] ref=Factory.encodeByte(ref0, (byte)31);
		
		final int qLen=query.length;
		final int rLen=ref.length;
		Visualizer viz=(output==null ? null : new Visualizer(output, 0, 0));
		
		// Banding parameters
		final int bandWidth=decideBandwidth(query, ref);
		// Initialize band limits for use outside main loop
		int bandStart=1, bandEnd=rLen;

		// Create arrays for current and previous rows
		byte[] prev=new byte[rLen+1], curr=new byte[rLen+1];
		Arrays.fill(curr, BAD);
		
		// Fill alignment matrix
		for(int i=1; i<=qLen; i++){
			// Calculate band boundaries 
			bandStart=Math.max(1, Math.min(i-bandWidth, rLen-bandWidth));
			bandEnd=Math.min(rLen, i+bandWidth);
			
            if(debug) {
                System.err.println("\nRow "+i+": bw="+bandWidth+"; j is "+bandStart+" to "+bandEnd);
            }
			
			//Clear stale data to the left of the band
			curr[bandStart-1]=BAD;

			// Clear first column score
			curr[0]=(byte)(i*INS);//TODO: This needs to be changed for relative scoring
			
			//Cache the query
			final byte q=query[i-1];
			
			if(Shared.SIMD) {
				//shared.SIMDAlign.alignBandVector(q, ref, bandStart, bandEnd, prev, curr, MATCH, N_SCORE, SUB, INS);
				assert(false) : "TODO";
			}else {

				// Process only cells within the band
				for(int j=bandStart; j<=bandEnd; j++){
					final byte r=ref[j-1];

					// Branchless score calculation
					final boolean isMatch=(q==r);
					final boolean hasN=((q|r)>=15);
					final byte scoreAdd=isMatch ? MATCH : (hasN ? N_SCORE : SUB);

					// Read adjacent scores
					final byte pj1=prev[j-1], pj=prev[j];
					final byte diagScore=(byte)(pj1+scoreAdd);// Match/Sub
					final byte upScore=(byte)(pj+INS);

					// Find max using conditional expressions
					final byte maxDiagUp=(byte)Math.max(diagScore, upScore);//This is fine

					curr[j]=maxDiagUp;
		            
		            // Output debug info if debug flag is set
		            if(debug) {
		                System.err.println("\nCell i="+i+", j="+j+": maxDiagUp="+maxDiagUp);
		                System.err.println("pj1="+pj1+", pj="+pj+", scoreAdd="+scoreAdd+", INS="+INS);
		                System.err.println("match="+isMatch+", diag="+diagScore+", up="+upScore);
		            }
				}
			}
			
			//Tail loop for deletions
			byte leftCell=curr[bandStart-1];
			for(int j=bandStart; j<=bandEnd; j++){
				final byte maxDiagUp=curr[j];
				final byte leftScore=(byte)(leftCell+DEL);
				leftCell=(byte)Math.max(maxDiagUp, leftScore);
				curr[j]=leftCell;

		        
				// Output debug info if debug flag is set
	            if(debug) {
	                System.err.println("\nCell i="+i+", j="+j+": score="+leftCell);
	                System.err.println("maxDiagUp="+maxDiagUp+", leftScore="+leftScore);
	            }
			}
	        
	        // Output debug info on row maxima
	        if(debug) {//TODO
	        	System.err.println("\nRow "+i+" scores: "+Arrays.toString(curr));
//	            int rowAbsScore=maxScore-MIDPOINT+totalAdjustment;
//	            System.err.println("\nRow i="+i+": rowScore="+maxScore+", rowAbsScore="+rowAbsScore);
	        }
			
			if(viz!=null) {viz.print(curr, bandStart, bandEnd, rLen);}
			if(loops>=0) {loops+=(bandEnd-bandStart+1);}

			// Swap rows
			byte[] temp=prev;
			prev=curr;
			curr=temp;
		}
		if(viz!=null) {viz.shutdown();}
		return postprocess(prev, qLen, bandStart, bandEnd, posVector);
	}

	/**
	 * Use alignment information to calculate identity and starting coordinate.
	 * @param prev Most recent score row
	 * @param qLen Query length
	 * @param bandStart Beginning of score band for the previous row
	 * @param bandEnd End of score band for the previous row
	 * @param posVector Optional array for returning reference start/stop coordinates.
	 * @return Identity
	 */
	private static final float postprocess(byte[] prev, int qLen, int bandStart, int bandEnd, int[] posVector) {
        
		// Find best score outside of main loop
		long maxScore=Long.MIN_VALUE;
		int maxPos=bandEnd;
		for(int j=bandStart; j<=bandEnd; j++){
			long score=prev[j];
			if(score>maxScore){
				maxScore=score;
				maxPos=j;
			}
		}
        if(debug) {
            System.err.println("\npostprocess: prev="+Arrays.toString(prev)+", qLen="+qLen+", bandStart="+bandStart+", bandEnd="+bandEnd);
            System.err.println("maxScore="+maxScore+", maxPos="+maxPos);
        }

		// Extract alignment information
		final int originPos=0;//Unknown
		final int endPos=maxPos;
		if(posVector!=null){
			posVector[0]=originPos;
			posVector[1]=endPos-1;
		}

		float identity=(maxScore+qLen)/(2f*qLen);
		
		return identity;
	}

	/**
	 * Lightweight wrapper for aligning to a window of the reference.
	 * @param query Query sequence
	 * @param ref Reference sequence
	 * @param posVector Optional int[2] for returning {rStart, rStop} of the optimal alignment.
	 * If the posVector is null, sequences may be swapped so that the query is shorter.
	 * @param rStart Alignment window start.
	 * @param to Alignment window stop.
	 * @return Identity (0.0-1.0).
	 */
	public static final float alignStatic(final byte[] query, final byte[] ref, 
			final int[] posVector, int refStart, int refEnd) {
		refStart=Math.max(refStart, 0);
		refEnd=Math.min(refEnd, ref.length-1);
		final int rlen=refEnd-refStart+1;
		final byte[] region=(rlen==ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd));
		final float id=alignStatic(query, region, posVector);
		if(posVector!=null) {
			posVector[0]+=refStart;
			posVector[1]+=refStart;
		}
		return id;
	}

	static long loops=-1; //-1 disables.  Be sure to disable this prior to release!
	public long loops() {return loops;}
	public void setLoops(long x) {loops=x;}
	public static String output=null;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	// Scoring constants
	private static final byte MATCH=1;
	private static final byte SUB=-1;
	private static final byte INS=-1;
	private static final byte DEL=-1;
	private static final byte N_SCORE=0;
	private static final byte BAD=Byte.MIN_VALUE/2;

	// Run modes
	private static final boolean PRINT_OPS=false;
	public static boolean debug=true;
//	public static boolean GLOBAL=false; //Cannot handle global alignments

}
