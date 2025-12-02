package ml;

import java.util.Random;

/**
 * @author Brandon Imstepf
 * Last edited: @date 7-21-2025
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
        // Gradient with respect to inputs
        float[] gradientOut = new float[inputSize];
        
        // Apply ReLU derivative
        float[] reluGradient = new float[outputSize];
        for(int i = 0; i < outputSize; i++) {
            reluGradient[i] = lastOutput[i] > 0 ? gradientIn[i] : 0;
        }
        
        // Backpropagate through weights
        for(int j = 0; j < inputSize; j++) {
            float sum = 0;
            for(int i = 0; i < outputSize; i++) {
                sum += weights[i][j] * reluGradient[i];
            }
            gradientOut[j] = sum;
        }
        
        // Store gradients for weight updates
        if(weightGradients == null) {
            weightGradients = new float[outputSize][inputSize];
            biasGradients = new float[outputSize];
        }
        
        // Accumulate gradients
        for(int i = 0; i < outputSize; i++) {
            for(int j = 0; j < inputSize; j++) {
                weightGradients[i][j] += reluGradient[i] * lastInput[j];
            }
            biasGradients[i] += reluGradient[i];
        }
        
        // CLIP GRADIENTS AFTER ALL ACCUMULATION IS DONE
        for(int i = 0; i < outputSize; i++) {
            for(int j = 0; j < inputSize; j++) {
                // Clip weight gradients
                weightGradients[i][j] = Math.max(-5.0f, Math.min(5.0f, weightGradients[i][j]));
            }
            // Clip bias gradients
            biasGradients[i] = Math.max(-5.0f, Math.min(5.0f, biasGradients[i]));
        }
        
        return gradientOut;  // Note: it's gradientOut, not outputGradient
    }

    @Override
    public void updateWeights(float learningRate) {
        if(weightGradients == null) return;
        
        // Update weights and biases
        for(int i = 0; i < outputSize; i++) {
            for(int j = 0; j < inputSize; j++) {
                weights[i][j] -= learningRate * weightGradients[i][j];
                weightGradients[i][j] = 0;  // Reset gradient
            }
            bias[i] -= learningRate * biasGradients[i];
            biasGradients[i] = 0;  // Reset gradient
        }
    }


    
    @Override
    public int getOutputSize() {
        return outputSize;
    }

    public int getParameterCount() {
        return outputSize * inputSize + outputSize;
    }
    
    int inputSize;
    int outputSize;
    float[][] weights;
    float[] bias;
    float[][] weightGradients;
    float[] biasGradients;
    boolean useReLU = true; //
}