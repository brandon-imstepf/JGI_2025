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
        lastZ = new float[]{sum};  // Store pre-activation for backward pass
        return output;
    }
    
    @Override
    public float[] backward(float[] gradientIn) {
        // gradientIn[0] contains the loss gradient
        // We need to multiply by sigmoid derivative: sigmoid'(z) = sigmoid(z) * (1 - sigmoid(z))
        float sigmoidDerivative = lastOutput[0] * (1.0f - lastOutput[0]);
        float outputGradient = gradientIn[0] * sigmoidDerivative;
        
        // Gradient with respect to inputs
        float[] gradientOut = new float[inputSize];
        
        for(int j = 0; j < inputSize; j++) {
            gradientOut[j] = outputGradient * weights[0][j];
        }
        
        // Store gradients for weight updates
        if(weightGradients == null) {
            weightGradients = new float[outputSize][inputSize];
            biasGradients = new float[outputSize];
        }
        
        // Accumulate gradients
        for(int j = 0; j < inputSize; j++) {
            weightGradients[0][j] += outputGradient * lastInput[j];
        }
        biasGradients[0] += outputGradient;
        
        // Clip gradients here too
        for(int j = 0; j < inputSize; j++) {
            weightGradients[0][j] = Math.max(-5.0f, Math.min(5.0f, weightGradients[0][j]));
        }
        biasGradients[0] = Math.max(-5.0f, Math.min(5.0f, biasGradients[0]));
        
        return gradientOut;
    }
    
    private float[] lastZ;  // Pre-activation values
}