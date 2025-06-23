package ml;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.LineParser1;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.IntList;

public class DataLoader {
	
	private DataLoader(String fname_){
		fnames=new File(fname_).exists() ? 
				new String[] {fname_} : Tools.commaPattern.split(fname_);
	}

	
	public static Matrix[] split(Matrix m, float fraction, int maxLines1, boolean exclusive) {
		if(!exclusive) {fraction=1;}//For convenience...
		
		Matrix[] array=new Matrix[] {new Matrix(), new Matrix()};
		@SuppressWarnings("unchecked")
		ArrayList<float[]>[] inputs=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
		@SuppressWarnings("unchecked")
		ArrayList<float[]>[] outputs=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
		@SuppressWarnings("unchecked")
		ArrayList<float[]>[] weights=new ArrayList[] {new ArrayList<float[]>(), new ArrayList<float[]>()};
		
		for(Matrix n : array) {
			n.dims=m.dims;
			n.numInputs=m.numInputs;
			n.numOutputs=m.numOutputs;
			n.columns=m.columns;
		}
		
		Random randy=new Random(0);
		for(int i=0; i<m.inputs.length; i++) {
			float[] in=m.inputs[i], out=m.outputs[i], wt=m.weights[i];
			int positive=(out[0]>=0.5f ? 1 : 0);
			int negative=positive^1;
			int pick=(randy.nextFloat()>=fraction || inputs[1].size()>=maxLines1) ? 0 : 1;
			inputs[pick].add(in);
			outputs[pick].add(out);
			weights[pick].add(wt);
			array[pick].numNegative+=negative;
			array[pick].numPositive+=positive;
			array[pick].validLines++;
		}
		
		for(int set=0; set<array.length; set++){
			Matrix n=array[set];
			n.inputs=new float[inputs[set].size()][];
			n.outputs=new float[outputs[set].size()][];
			n.weights=new float[weights[set].size()][];
			for(int i=0; i<n.inputs.length; i++) {
				n.inputs[i]=inputs[set].get(i);
				n.outputs[i]=outputs[set].get(i);
				n.weights[i]=weights[set].get(i);
			}
			n.data=new float[][][] {n.inputs, n.outputs, n.weights};
			if(exclusive || set>0) {n.detectRange();}
		}
		
		if(!exclusive) {array[0]=m;}
		
		return array;
	}
	
	public static SampleSet[] load(String fname, int maxLines0, boolean shuffleRaw, 
			float splitFraction, int maxLines1, boolean exclusive, float balance) {
		DataLoader dl=new DataLoader(fname);
		dl.load(maxLines0, shuffleRaw || splitFraction>0, balance);
		
		lastValidLines=dl.validLines;
		lastInvalidLines=dl.invalidLines;
//		assert(false) : dl.matrix.weights[0][0];
		if(splitFraction<=0) {
			return new SampleSet[] {new SampleSet(dl.matrix)};
		}else {
			Matrix[] array=split(dl.matrix, splitFraction, maxLines1, exclusive);
			return new SampleSet[] {new SampleSet(array[0]), new SampleSet(array[1])};
		}
	}
	
	private void load(final int maxLines, final boolean shuffleRaw, final float balance) {
//		bf.reset();
		matrix=new Matrix();
		ArrayList<float[]> inputList=new ArrayList<float[]>();
		ArrayList<float[]> outputList=new ArrayList<float[]>();
		ArrayList<float[]> weightList=new ArrayList<float[]>();
		byte[] s=null;
		final int max=(shuffleRaw ? Shared.MAX_ARRAY_LEN : maxLines);
		
		for(String f : fnames) {
			FileFormat ff=FileFormat.testInput(f, FileFormat.TEXT, null, true, false);
			ByteFile bf=ByteFile.makeByteFile(ff);
			for(s=bf.nextLine(); s!=null && validLines<max; s=bf.nextLine()){
				if(s.length>0) {
					if(s[0]=='#') {//TODO: When processing multiple files, 
						if(Tools.startsWith(s, "#dims")) {
							matrix.dims=parseIntArray(s, delimiter, true);
							matrix.numInputs=matrix.dims[0];
//							matrix.numOutputs=matrix.dims[matrix.dims.length-1];
							matrix.numOutputs=matrix.dims[1];
							weighted=(matrix.dims.length>2 && matrix.dims[2]==1);
							assert(matrix.dims.length>1) : matrix.dims.length+", "+Arrays.toString(matrix.dims)+", '"+new String(s)+"'";
						}else if(Tools.startsWith(s, "#inputs")) {
							matrix.numInputs=parseInt(s);
						}else if(Tools.startsWith(s, "#outputs")) {
							matrix.numOutputs=parseInt(s);
						}else if(Tools.startsWith(s, "##")) {
							matrix.columns=new ArrayList<String>(Arrays.asList(new String(s).split("\t")));
							matrix.columns.set(0, matrix.columns.get(0).substring(2));//Trim ##
						}else {
							//comment
						}
					}else {
						if(matrix.numInputs==0) {
							int terms=Tools.split(s, 0, (byte)'\t').size();
							matrix.numOutputs=1;
							matrix.numInputs=terms-(matrix.numOutputs+(weighted ? 1 : 0));
							System.err.println("Inferring "+matrix.numInputs+" inputs, "+matrix.numOutputs+" output, "+(weighted ? 1 : 0)+" weights.");
						}
						assert(matrix.numInputs>0 & matrix.numOutputs>0) : 
							"Number of inputs and outputs must be in data file header, e.g. '#inputs 5'";
						float[] inputs=new float[matrix.numInputs];
						float[] outputs=new float[matrix.numOutputs];
						float[] weights=new float[] {1};
//						System.err.println("Attempting to parse line; i="+inputs.length+", o="+outputs.length+", w="+weights.length+", weighted="+weighted);
						boolean valid=parseDataLine(s, inputs, outputs, weights);
//						assert(false);
						if(valid) {
							inputList.add(inputs);
							outputList.add(outputs);
							weightList.add(weights);
							validLines++;
						}else {
							invalidLines++;
						}
					}
				}
			}
			bf.close();
			if(validLines>=max) {break;}
		}
		if(shuffleRaw) {shuffle(inputList, outputList, weightList, maxLines);}
		if(balance>0) {balance(inputList, outputList, weightList, balance);}
		matrix.inputs=new float[inputList.size()][];
		matrix.outputs=new float[outputList.size()][];
		matrix.weights=new float[weightList.size()][];
		for(int i=0; i<matrix.inputs.length; i++) {
			matrix.inputs[i]=inputList.get(i);
			matrix.outputs[i]=outputList.get(i);
			matrix.weights[i]=weightList.get(i);
		}
		matrix.data=new float[][][] {matrix.inputs, matrix.outputs, matrix.weights};
		matrix.initializeRange();
//		assert(false) : validLines+", "+invalidLines+", "+matrix.inputs.length+"\n"+matrix.toString();
//		assert(false) : Arrays.toString(matrix.inputs[0])+"\n"+Arrays.toString(matrix.outputs[0])+
//			", "+matrix.targetOutputRangeMax;
	}
	
	private static void shuffle(ArrayList<float[]> inputList, ArrayList<float[]> outputList, ArrayList<float[]> weightList, int maxLines) {
		final int size=inputList.size();
		ArrayList<Triple> list=new ArrayList<Triple>(inputList.size());
		for(int i=0; i<inputList.size(); i++){
			Triple p=new Triple(inputList.get(i), outputList.get(i), weightList.get(i));
			list.add(p);
		}
		Random randy=new Random(SampleSet.shuffleSeed);
		Collections.shuffle(list, randy);
		inputList.clear();
		outputList.clear();
		weightList.clear();
		for(int i=0, lim=Tools.min(maxLines, list.size()); i<lim; i++) {
			Triple p=list.get(i);
			inputList.add(p.in);
			outputList.add(p.out);
			weightList.add(p.w);
		}
		assert(inputList.size()==size);
		assert(outputList.size()==size);
		assert(weightList.size()==size);
	}
	
	private static void balance(ArrayList<float[]> inputList, ArrayList<float[]> outputList, ArrayList<float[]> weightList, float mult) {
		int pos=0, neg=0;
		assert(mult>0 && mult<=1);
		for(float[] out : outputList) {
			if(out[0]>=0.5f) {pos++;} else {neg++;}
		}
		final int target=(int)(Tools.max(pos, neg)*mult);
		if(pos>=target && neg>=target) {return;}
		assert(pos>=target || neg>=target);
		if(pos<1 || neg<1) {
			throw new RuntimeException("Can't balance with zero examples: pos="+pos+", neg="+neg);
		}
		for(int i=0; i<inputList.size() && (pos<target || neg<target); i++) {
			float[] in=inputList.get(i), out=outputList.get(i), weight=weightList.get(i);
			if(out[0]>=0.5 && pos<target) {
				pos++;
				inputList.add(in);
				outputList.add(out);
				weightList.add(weight);
			}else if(out[0]<0.5f && neg<target) {
				pos++;
				inputList.add(in);
				outputList.add(out);
				weightList.add(weight);
			}
		}
		shuffle(inputList, outputList, weightList, inputList.size());
	}
	
	boolean parseDataLine(byte[] line, float[] inputs, float[] outputs, float[] weights) {
		lp.set(line);
		int pos=0;
		for(int i=0; i<inputs.length; i++) {
			inputs[i]=lp.parseFloat(pos);
			pos++;
		}
		if(weighted) {
			weights[0]=lp.parseFloat(pos);
			pos++;
		}else {
			weights[0]=1;
		}
		for(int i=0; i<outputs.length; i++) {
			outputs[i]=lp.parseFloat(pos);
			pos++;
		}
		assert(pos==lp.terms()) : "\nExtra characters for line '"+new String(line)+
			"'; numInputs="+matrix.numInputs+", numOutputs="+matrix.numOutputs+
			", "+pos+", "+line.length+"\n"+Arrays.toString(inputs)+"\n"+Arrays.toString(outputs)+"\n";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int parseInt(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseInt(line, idx+1, line.length);
	}
	
	public static int[] parseIntArray(final byte[] line, final byte delimiter, boolean parseTitle){
		int a=0, b=0;
		IntList list=new IntList(3);
		
		if(parseTitle) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing Title: "+new String(line);
			b++;
			a=b;
		}
		
		while(a<line.length) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing element "+list.size+": '"+new String(line)+"'";
			int x=Parse.parseInt(line, a, b);
//			assert(x>0) : new String(line);
			list.add(x);
			b++;
			a=b;
		}
		return list.toArray();
	}
	
	/*--------------------------------------------------------------*/
	
	private static class Triple{
		Triple(float[] in_, float[] out_, float[] w_){
			in=in_;
			out=out_;
			w=w_;
		}
		final float[] in, out, w;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	LineParser1 lp=new LineParser1(delimiter);
	
	String fnames[];
	Matrix matrix;
//	SampleSet ss;
	int pos=0;
	long validLines=0;
	long invalidLines=0;
	
	static long lastValidLines=0;
	static long lastInvalidLines=0;
	static boolean weighted=false;
	
	public static final byte delimiter='\t';
	public static boolean IGNORE_BAD_LINES=false;
	
	
}
