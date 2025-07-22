package ml;

/**
 * @author Brandon Imstepf
 * Last edited: @date 7-21-2025
 * Max pooling layer - reduces dimensionality by taking maximum values.
 * Handles multiple input channels/filters correctly.
 */
public class MaxPoolingLayer extends Layer {
    
    public MaxPoolingLayer(String name, int inputChannels, int inputLength, int poolSize) {
        this.name = name;
        this.inputChannels = inputChannels;
        this.inputLength = inputLength;
        this.poolSize = poolSize;
        this.outputLength = inputLength / poolSize;
        this.outputSize = inputChannels * outputLength;
        
        System.err.println(name + ": " + (inputChannels * inputLength) + " -> " + 
                          outputSize + " (channels=" + inputChannels + 
                          ", pool=" + poolSize + ")");
    }
        
    @Override
    public float[] forward(float[] input) {
        lastInput = input;
        float[] output = new float[inputChannels * outputLength];  
        maxIndices = new int[inputChannels * outputLength];      
        
        // For each channel
        for(int c = 0; c < inputChannels; c++) {
            int channelOffset = c * inputLength;
            int outputOffset = c * outputLength;
            
            // For each output position
            for(int i = 0; i < outputLength; i++) {
                int startIdx = i * poolSize;
                float maxVal = input[channelOffset + startIdx];
                int maxIdx = startIdx;
                
                // Find max in pool window
                for(int j = 1; j < poolSize && startIdx + j < inputLength; j++) {
                    float val = input[channelOffset + startIdx + j];
                    if(val > maxVal) {
                        maxVal = val;
                        maxIdx = startIdx + j;
                    }
                }
                
                output[outputOffset + i] = maxVal;
                maxIndices[outputOffset + i] = maxIdx;
            }
        }
        
        lastOutput = output;
        return output;
    }
    
    @Override
    public float[] backward(float[] gradientIn) {
        float[] gradientOut = new float[inputChannels * inputLength];
        
        // Route gradients back to max positions
        for(int c = 0; c < inputChannels; c++) {
            for(int i = 0; i < outputLength; i++) {
                int maxIdx = maxIndices[c * outputLength + i];
                int inputIdx = c * inputLength + maxIdx;
                int outputIdx = c * outputLength + i;
                gradientOut[inputIdx] = gradientIn[outputIdx];
            }
        }
        
        return gradientOut;
    }

    @Override
    public void updateWeights(float learningRate) {
        // No weights to update in pooling layer
    }
    
    @Override
    public int getOutputSize() {
        return outputSize;
    }
    
    // Need to know number of channels for dense layer
    public int getOutputChannels() { return inputChannels; }
    public int getOutputLength() { return outputLength; }
    
    int inputChannels;
    int inputLength;
    int poolSize;
    int outputLength;
    int outputSize;
    int[] maxIndices;  // Store indices of max values for backpropagation
}