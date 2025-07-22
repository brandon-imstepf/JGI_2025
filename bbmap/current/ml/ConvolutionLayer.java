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
        // Initialize gradient for input
        float[] gradientOut = new float[inputChannels * inputLength];
        
        // Initialize weight gradients if needed
        if(weightGradients == null) {
            weightGradients = new float[numFilters][inputChannels][filterSize];
            biasGradients = new float[numFilters];
        }
        
        // For each filter
        for(int f = 0; f < numFilters; f++) {
            // For each output position
            for(int i = 0; i < outputLength; i++) {
                // Get gradient from this output position
                float grad = gradientIn[f * outputLength + i];
                
                // Skip if gradient is zero (ReLU killed it)
                if(grad == 0) continue;
                
                // For each input channel
                for(int c = 0; c < inputChannels; c++) {
                    // For each position in the filter
                    for(int k = 0; k < filterSize; k++) {
                        int inputIdx = c * inputLength + i + k;
                        
                        // Gradient with respect to input
                        gradientOut[inputIdx] += grad * weights[f][c][k];
                        
                        // Gradient with respect to weights
                        weightGradients[f][c][k] += grad * lastInput[inputIdx];
                    }
                }
                
                // Gradient with respect to bias
                biasGradients[f] += grad;
            }
        }
        
        // CLIP GRADIENTS AFTER ALL ACCUMULATION IS DONE
        for(int f = 0; f < numFilters; f++) {
            for(int c = 0; c < inputChannels; c++) {
                for(int k = 0; k < filterSize; k++) {
                    // Clip gradients to prevent explosion
                    weightGradients[f][c][k] = Math.max(-5.0f, Math.min(5.0f, weightGradients[f][c][k]));
                }
            }
            // Also clip bias gradients
            biasGradients[f] = Math.max(-5.0f, Math.min(5.0f, biasGradients[f]));
        }
        
        return gradientOut;
    }

    @Override
    public void updateWeights(float learningRate) {
        if(weightGradients == null) return;
        
        // Update weights
        for(int f = 0; f < numFilters; f++) {
            for(int c = 0; c < inputChannels; c++) {
                for(int k = 0; k < filterSize; k++) {
                    weights[f][c][k] -= learningRate * weightGradients[f][c][k];
                    weightGradients[f][c][k] = 0;  // Reset
                }
            }
            bias[f] -= learningRate * biasGradients[f];
            biasGradients[f] = 0;  // Reset
        }
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
    
    int inputChannels;
    int inputLength;
    int numFilters;
    int filterSize;
    int outputLength;
    float[][][] weights;  // [numFilters][inputChannels][filterSize]
    float[] bias;
    float[][][] weightGradients;
    float[] biasGradients;
    int stride = 1; // default stride
}