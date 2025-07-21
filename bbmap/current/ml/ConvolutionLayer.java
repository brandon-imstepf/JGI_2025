package ml;

import java.util.Random;

/**
 * @author Brandon Imstepf
 * @date 7-11-2025
 * 1D Convolution layer for sequence data.
 * Handles multiple input channels from previous layers.
 */
public class ConvolutionLayer extends Layer {
    
    public ConvolutionLayer(String name, int inputChannels, int inputLength, 
                           int numFilters, int filterSize) {
        this.name = name;
        this.inputChannels = inputChannels;
        this.inputLength = inputLength;
        this.numFilters = numFilters;
        this.filterSize = filterSize;
        this.outputLength = inputLength - filterSize + 1; // Valid convolution
        
        // Initialize weights: [numFilters][inputChannels][filterSize]
        Random rand = new Random();
        weights = new float[numFilters][inputChannels][filterSize];
        bias = new float[numFilters];
        
        // Xavier initialization
        float scale = (float)Math.sqrt(2.0 / (inputChannels * filterSize));
        for(int f = 0; f < numFilters; f++) {
            for(int c = 0; c < inputChannels; c++) {
                for(int i = 0; i < filterSize; i++) {
                    weights[f][c][i] = (float)(rand.nextGaussian() * scale);
                }
            }
            bias[f] = 0;
        }
        
        System.err.println(name + ": " + (inputChannels * inputLength) + " -> " + 
                          (numFilters * outputLength) + 
                          " (channels=" + inputChannels + "→" + numFilters + 
                          ", length=" + inputLength + "→" + outputLength + ")");
    }
    
    @Override
    public float[] forward(float[] input) {
        lastInput = input;
        float[] output = new float[numFilters * outputLength];
        
        // For each output filter
        for(int f = 0; f < numFilters; f++) {
            // For each output position
            for(int pos = 0; pos < outputLength; pos++) {
                float sum = bias[f];
                
                // Sum across all input channels
                for(int c = 0; c < inputChannels; c++) {
                    int channelOffset = c * inputLength;
                    // Convolve with the filter
                    for(int i = 0; i < filterSize; i++) {
                        sum += input[channelOffset + pos + i] * weights[f][c][i];
                    }
                }
                
                // ReLU activation
                output[f * outputLength + pos] = Math.max(0, sum);
            }
        }
        
        lastOutput = output;
        return output;
    }
    
    @Override
    public float[] backward(float[] gradientIn) {
        // TODO: Implement in Phase 3
        return new float[inputChannels * inputLength];
    }
    
    @Override
    public void updateWeights(float learningRate) {
        // TODO: Implement in Phase 3
    }
    
    @Override
    public int getOutputSize() {
        return numFilters * outputLength;
    }

    public int getParameterCount() {
        return numFilters * inputChannels * filterSize + numFilters;
    }    
    
    public int getOutputChannels() { return numFilters; }
    public int getOutputLength() { return outputLength; }
    
    private final int inputChannels;
    private final int inputLength;
    private final int numFilters;
    private final int filterSize;
    private final int outputLength;
    private float[][][] weights;  // [numFilters][inputChannels][filterSize]
    private float[] bias;
}