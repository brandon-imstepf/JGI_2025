package ml;

/**
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
        float[] output = new float[outputSize];
        
        // For each channel/filter
        for(int c = 0; c < inputChannels; c++) {
            int channelOffset = c * inputLength;
            int outputOffset = c * outputLength;
            
            // For each output position in this channel
            for(int i = 0; i < outputLength; i++) {
                float maxVal = Float.NEGATIVE_INFINITY;
                // Find max in the pooling window
                for(int j = 0; j < poolSize; j++) {
                    int idx = channelOffset + (i * poolSize + j);
                    if(idx < channelOffset + inputLength) {
                        maxVal = Math.max(maxVal, input[idx]);
                    }
                }
                output[outputOffset + i] = maxVal;
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
        // No weights to update in pooling layer
    }
    
    @Override
    public int getOutputSize() {
        return outputSize;
    }
    
    // Need to know number of channels for dense layer
    public int getOutputChannels() { return inputChannels; }
    public int getOutputLength() { return outputLength; }
    
    private final int inputChannels;
    private final int inputLength;
    private final int poolSize;
    private final int outputLength;
    private final int outputSize;
}