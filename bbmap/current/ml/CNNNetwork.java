package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Convolutional Neural Network implementation.
 * 
 * @author Brandon Imstepf
 * @date 7-11-2025
 */
public class CNNNetwork {
    
    /*--------------------------------------------------------------*/
    /*----------------        Initialization        ----------------*/
    /*--------------------------------------------------------------*/
    
    /**
     * Constructor.
     * @param inputs Number of input features
     * @param outputs Number of output classes
     */
    public CNNNetwork(int inputs, int outputs) {
        this.numInputs = inputs;
        this.numOutputs = outputs;
        this.outstream = System.err;
        
        outstream.println("CNNNetwork initialized with " + inputs + " inputs and " + outputs + " outputs");
    }
    
    /**
     * Initialize the network architecture.
     */
    public void initialize() {
        outstream.println("Building network architecture...");
        
        layers = new ArrayList<>();
        int currentSize = numInputs;
        int currentChannels = 1;  // Start with 1 channel
        int currentLength = numInputs;
        
        // Add convolutional layers
        for(int i = 0; i < filterCounts.length; i++) {
            // Convolution layer
            int filterSize = i < filterSizes.length ? filterSizes[i] : 3;
            ConvolutionLayer conv = new ConvolutionLayer(
                "Conv" + (i+1), currentChannels, currentLength, 
                filterCounts[i], filterSize
            );
            layers.add(conv);
            currentChannels = conv.getOutputChannels();
            currentLength = conv.getOutputLength();
            currentSize = conv.getOutputSize();
            
            // Pooling layer (if specified)
            if(i < poolSizes.length && poolSizes[i] > 1) {
                MaxPoolingLayer pool = new MaxPoolingLayer(
                    "Pool" + (i+1), currentChannels, currentLength, poolSizes[i]
                );
                layers.add(pool);
                currentLength = pool.getOutputLength();
                currentSize = pool.getOutputSize();
            }
        }
        
        // Add dense layers
        for(int i = 0; i < denseLayers.length; i++) {
            DenseLayer dense = new DenseLayer(
                "Dense" + (i+1), currentSize, denseLayers[i]
            );
            layers.add(dense);
            currentSize = dense.getOutputSize();
        }
        
        // Final output layer (sigmoid for binary classification)
        OutputLayer outputLayer = new OutputLayer("Output", currentSize);
        layers.add(outputLayer);
        
        outstream.println("\nNetwork architecture built with " + layers.size() + " layers");
        outstream.println("Total parameters: " + countParameters());
    }

    private int countParameters() {
        int total = 0;
        for(Layer layer : layers) {
            if(layer instanceof ConvolutionLayer) {
                ConvolutionLayer conv = (ConvolutionLayer)layer;
                // weights: filters * channels * filterSize + biases
                total += conv.getParameterCount();
            } else if(layer instanceof DenseLayer) {
                DenseLayer dense = (DenseLayer)layer;
                total += dense.getParameterCount();
            }
            // MaxPooling has no parameters
        }
        return total;
    }
    
    /**
     * Train the network.
     */
    public void train(SampleSet trainData, SampleSet valData) {
        outstream.println("You're calling CNNNetwork.train()");
        outstream.println("Training samples: " + trainData.samples.length);
        outstream.println("Validation samples: " + (valData != null ? valData.samples.length : 0));
    
        // Update train method to test forward pass

        outstream.println("Testing forward pass...");
        
        // Test on first sample
        Sample firstSample = trainData.samples[0];
        float[] output = forward(firstSample.in);  // <-- Changed to use .in
        
        outstream.println("First sample forward pass:");
        outstream.println("  Input size: " + firstSample.in.length);
        outstream.println("  Output: " + output[0]);
        outstream.println("  Target: " + firstSample.goal[0]);
        outstream.println("  Initial prediction: " + (output[0] > 0.5f ? "positive" : "negative"));
    
        // TODO: Implement actual training loop

    }

    /**
     * Forward pass through the entire network.
     * @param input The input features
     * @return The network output (probability for binary classification)
     */
    public float[] forward(float[] input) {
        float[] current = input;
        for(Layer layer : layers) {
            current = layer.forward(current);
        }
        return current;
    }

    /**
     * Set the network architecture.
     */
    public void setArchitecture(int[] filterCounts, int[] filterSizes, 
                               int[] poolSizes, int[] denseLayers) {
        this.filterCounts = filterCounts;
        this.filterSizes = filterSizes;
        this.poolSizes = poolSizes;
        this.denseLayers = denseLayers;
        
        outstream.println("Architecture set:");
        outstream.println("  Conv layers: " + filterCounts.length);
        outstream.println("  Filter counts: " + Arrays.toString(filterCounts));
        outstream.println("  Filter sizes: " + Arrays.toString(filterSizes));
        outstream.println("  Pool sizes: " + Arrays.toString(poolSizes));
        outstream.println("  Dense layers: " + Arrays.toString(denseLayers));
    }

    /**
     * Set training parameters.
     */
    public void setTrainingParams(int epochs, int batchSize, 
                                 float learningRate, float dropout) {
        this.epochs = epochs;
        this.batchSize = batchSize;
        this.learningRate = learningRate;
        this.dropout = dropout;
        
        outstream.println("Training parameters:");
        outstream.println("  Epochs: " + epochs);
        outstream.println("  Batch size: " + batchSize);
        outstream.println("  Learning rate: " + learningRate);
        outstream.println("  Dropout: " + dropout);
    }
    
    /*--------------------------------------------------------------*/
    /*----------------            Fields            ----------------*/
    /*--------------------------------------------------------------*/
    
    private final int numInputs; // Number of input features
    private final int numOutputs; // Number of input features and output classes
    private PrintStream outstream; // Output stream for logging
    private int[] filterCounts; // Number of filters for each convolution layer
    private int[] filterSizes; // Filter sizes for convolution layers
    private int[] poolSizes; // Pool sizes for max pooling layers
    private int[] denseLayers; // Dense layers after convolutions
    private int epochs; // Number of training epochs
    private int batchSize; // Batch size for training
    private float learningRate; // Learning rate for weight updates
    private float dropout; // Dropout rate for regularization
    private ArrayList<Layer> layers; // List of layers in the network
}