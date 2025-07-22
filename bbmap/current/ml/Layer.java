package ml;

/**
 * @author Brandon Imstepf
 * @date 7-11-2025
 * Base class for neural network layers.
 * Following BBTools pattern of simple, efficient classes.
 */
public abstract class Layer {
    
    /** Forward pass - processes input and produces output */
    public abstract float[] forward(float[] input);
    
    /** Get the output size of this layer */
    public abstract int getOutputSize();
    
    /** Layer name for debugging */
    protected String name;
    
    /** Store input for backward pass */
    protected float[] lastInput;
    
    /** Store output for backward pass */
    protected float[] lastOutput;

    /**
     * Backward pass - compute gradients.
     * @param gradientIn Gradient from the next layer
     * @return Gradient to pass to previous layer
     */
    public abstract float[] backward(float[] gradientIn);

    /**
     * Update weights using accumulated gradients.
     * @param learningRate Learning rate for gradient descent
     */
    public void updateWeights(float learningRate) {
        // Default: do nothing (for layers without parameters)
    }
}