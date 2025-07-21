package ml;

/**
 * Output layer for binary classification.
 * Uses sigmoid activation instead of ReLU.
 */
public class OutputLayer extends DenseLayer {
    
    public OutputLayer(String name, int inputSize) {
        super(name, inputSize, 1);  // Binary classification = 1 output
    }
    
    @Override
    public float[] forward(float[] input) {
        lastInput = input;
        float[] output = new float[1];
        
        // Linear combination
        float sum = bias[0];
        for(int j = 0; j < inputSize; j++) {
            sum += weights[0][j] * input[j];
        }
        
        // Sigmoid activation: 1 / (1 + e^-x)
        output[0] = 1.0f / (1.0f + (float)Math.exp(-sum));
        
        lastOutput = output;
        return output;
    }
}