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
    
    /** Backward pass - calculates gradients for backpropagation */
    public abstract float[] backward(float[] gradientIn);
    
    /** Update weights using calculated gradients */
    public abstract void updateWeights(float learningRate);
    
    /** Get the output size of this layer */
    public abstract int getOutputSize();
    
    /** Layer name for debugging */
    protected String name;
    
    /** Store input for backward pass */
    protected float[] lastInput;
    
    /** Store output for backward pass */
    protected float[] lastOutput;
}