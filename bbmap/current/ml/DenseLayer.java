package ml;

import java.util.Random;

/**
 * Fully connected (dense) layer.
 * Every input connects to every output.
 */
public class DenseLayer extends Layer {
    
    public DenseLayer(String name, int inputSize, int outputSize) {
        this.name = name;
        this.inputSize = inputSize;
        this.outputSize = outputSize;
        
        // Initialize weights and biases
        Random rand = new Random();
        weights = new float[outputSize][inputSize];
        bias = new float[outputSize];
        
        // Xavier initialization
        float scale = (float)Math.sqrt(2.0 / inputSize);
        for(int i = 0; i < outputSize; i++) {
            for(int j = 0; j < inputSize; j++) {
                weights[i][j] = (float)(rand.nextGaussian() * scale);
            }
            bias[i] = 0;
        }
        
        System.err.println(name + ": " + inputSize + " -> " + outputSize);
    }
    
    @Override
    public float[] forward(float[] input) {
        lastInput = input;
        float[] output = new float[outputSize];
        
        // Matrix multiplication: output = weights * input + bias
        for(int i = 0; i < outputSize; i++) {
            float sum = bias[i];
            for(int j = 0; j < inputSize; j++) {
                sum += weights[i][j] * input[j];
            }
            // ReLU activation for hidden layers
            output[i] = Math.max(0, sum);
        }
        
        lastOutput = output;
        return output;
    }
    
    @Override
    public float[] backward(float[] gradientIn) {
        // TODO: Implement in Phase 3
        return new float[inputSize];
    }
    
    @Override
    public void updateWeights(float learningRate) {
        // TODO: Implement in Phase 3
    }
    
    @Override
    public int getOutputSize() {
        return outputSize;
    }
    
    private final int inputSize;
    private final int outputSize;
    private float[][] weights;
    private float[] bias;
}