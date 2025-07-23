package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.PrintWriter;
import java.io.IOException;


/**
 * Convolutional Neural Network implementation.
 * 
 * @author Brandon Imstepf
 * Last Edited: @date 7-21-2025
 */
public class CNNNetwork {
    
    /*--------------------------------------------------------------*/
    /*----------------        Initialization        ----------------*/
    /*--------------------------------------------------------------*/

    private float baseLearningRate;
    private float currentLearningRate;
    
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
    
    /*--------------------------------------------------------------*/
    /*----------------        Training Methods       ----------------*/
    /*--------------------------------------------------------------*/

    public void train(SampleSet trainData, SampleSet valData) {
        this.baseLearningRate = this.learningRate;
        
        outstream.println("Starting training for " + epochs + " epochs...");
        outstream.println("Batch size: " + batchSize);
        outstream.println("Base learning rate: " + baseLearningRate);

        outstream.println();
        // Get sample counts
        int trainSize = trainData.samples.length;
        int valSize = valData.samples.length;
        
        // Calculate positive/negative counts for metrics
        int trainPositive = 0, trainNegative = 0;
        for(Sample s : trainData.samples) {
            if(s.goal[0] > 0.5f) trainPositive++;
            else trainNegative++;
        }
        
        // Progress tracking
        long startTime = System.currentTimeMillis();
        int totalBatches = 0;
        int samplesProcessed = 0;
        
        // Training loop
        for(int epoch = 0; epoch < epochs; epoch++) {
            currentLearningRate = baseLearningRate * (float)Math.pow(0.95, epoch);
            outstream.println("\n=== Epoch " + (epoch + 1) + "/" + epochs + " ===");
            outstream.println("Current learning rate: " + String.format("%.6f", currentLearningRate) +  " (base * 0.95^" + epoch + ")");

            float epochLoss = 0;
            int epochCorrect = 0;
            int epochTP = 0, epochFP = 0, epochTN = 0, epochFN = 0;
            
            // Shuffle training data
            shuffleSamples(trainData.samples);
            
            // Process mini-batches
            for(int b = 0; b < trainSize; b += batchSize) {
                int batchEnd = Math.min(b + batchSize, trainSize);
                float batchLoss = 0;
                int batchCorrect = 0;
                int batchTP = 0, batchFP = 0, batchTN = 0, batchFN = 0;
                
                // Process each sample in batch
                for(int i = b; i < batchEnd; i++) {
                    Sample sample = trainData.samples[i];
                    
                    // Forward pass
                    float[] output = forward(sample.in);
                    
            // Check output for NaN
                if(Float.isNaN(output[0]) || Float.isInfinite(output[0])) {
                    outstream.println("\n!!! NaN/Inf OUTPUT !!!");
                    outstream.println("Epoch: " + (epoch+1) + ", Batch: " + totalBatches);
                    outstream.println("Sample index: " + i);
                    outstream.println("Output value: " + output[0]);
                    outstream.println("Target value: " + sample.goal[0]);
                    
                    // Check some input values
                    outstream.print("Input sample (first 10): [");
                    for(int j = 0; j < Math.min(10, sample.in.length); j++) {
                        outstream.print(sample.in[j]);
                        if(j < 9) outstream.print(", ");
                    }
                    outstream.println("...]");
                    
                    // Skip this sample
                    continue;
                }
                
                // Calculate loss
                float loss = BinaryCrossEntropyLoss.calculateLoss(output[0], sample.goal[0]);
                
                // Check loss for NaN
                if(Float.isNaN(loss) || Float.isInfinite(loss)) {
                    outstream.println("\n!!! NaN/Inf LOSS !!!");
                    outstream.println("Epoch: " + (epoch+1) + ", Batch: " + totalBatches);
                    outstream.println("Sample index: " + i);
                    outstream.println("Output: " + output[0] + ", Target: " + sample.goal[0]);
                    outstream.println("Loss calculation breakdown:");
                    outstream.println("  -[" + sample.goal[0] + " * log(" + output[0] + ") + " + 
                                    (1-sample.goal[0]) + " * log(" + (1-output[0]) + ")]");
                    
                    // Early stopping if this happens in first epoch
                    if(epoch == 0 && totalBatches < 10) {
                        outstream.println("Early NaN detected - stopping training!");
                        return;
                    }
                    continue;
                }
                    batchLoss += loss;
                    epochLoss += loss;
                    
                    // Track accuracy and confusion matrix
                    boolean predicted = output[0] > 0.5f;
                    boolean actual = sample.goal[0] > 0.5f;
                    
                    if(predicted == actual) {
                        batchCorrect++;
                        epochCorrect++;
                    }
                    
                    // Update confusion matrix
                    if(actual && predicted) {
                        batchTP++; epochTP++;
                    } else if(!actual && !predicted) {
                        batchTN++; epochTN++;
                    } else if(!actual && predicted) {
                        batchFP++; epochFP++;
                    } else {
                        batchFN++; epochFN++;
                    }
                    
                    // Calculate loss gradient
                    float lossGrad = BinaryCrossEntropyLoss.calculateGradient(
                        output[0], sample.goal[0]);

                    // Check loss gradient for NaN
                    if(Float.isNaN(lossGrad) || Float.isInfinite(lossGrad)) {
                        outstream.println("\n!!! NaN/Inf GRADIENT !!!");
                        outstream.println("Loss gradient: " + lossGrad);
                        outstream.println("From output: " + output[0] + ", target: " + sample.goal[0]);
                        continue;
                    }
                    
                    // Backward pass
                    backward(new float[]{lossGrad});
                }
                
                // Update weights after batch
                updateWeights();
                
                totalBatches++;
                samplesProcessed += (batchEnd - b);
                
                // Print progress every N batches or at specific milestones
                boolean shouldPrint = false;
                int printInterval = 1000; // Every 1000 batches
                
                // Check if we should print (similar to BBTools pattern)
                if(totalBatches <= 10 || 
                totalBatches == 25 || totalBatches == 50 || totalBatches == 100 ||
                totalBatches == 250 || totalBatches == 500 ||
                (totalBatches % printInterval == 0)) {
                    shouldPrint = true;
                }
                
                if(shouldPrint) {
                    // Calculate metrics
                    float err = 1.0f - (float)epochCorrect / samplesProcessed;
                    float fpr = epochFP > 0 ? (float)epochFP / (epochFP + epochTN) : 0;
                    float fnr = epochFN > 0 ? (float)epochFN / (epochFN + epochTP) : 0;
                    float avgLoss = epochLoss / samplesProcessed;

                    // Sanity Check
                    if(err > 0.9) {
                        outstream.println("WARNING: Very high error rate: " + err);
                        outstream.println("  Correct: " + epochCorrect + " / " + samplesProcessed);
                        outstream.println("  This might indicate a problem!");
                    }
                                
                    // Time elapsed
                    long elapsed = System.currentTimeMillis() - startTime;
                    String timeStr = formatTime(elapsed);
                    
                    // Format batch number
                    String batchStr = totalBatches >= 1000 ? 
                        String.format("%dk", totalBatches/1000) : 
                        String.valueOf(totalBatches);
                    
                    outstream.println(String.format(
                        "Batch %-8s ERR= %.6f  LOSS= %.6f  FPR= %.6f  FNR= %.6f  %s",
                        batchStr + ":", err, avgLoss, fpr, fnr, timeStr));
                }
            }
            
            // End of epoch - validation
            float valLoss = 0;
            int valCorrect = 0;
            int valTP = 0, valFP = 0, valTN = 0, valFN = 0;
            
            for(int i = 0; i < valSize; i++) {
                Sample sample = valData.samples[i];
                float[] output = forward(sample.in);
                valLoss += BinaryCrossEntropyLoss.calculateLoss(
                    output[0], sample.goal[0]);
                
                boolean predicted = output[0] > 0.5f;
                boolean actual = sample.goal[0] > 0.5f;
                
                if(predicted == actual) valCorrect++;
                
                if(actual && predicted) valTP++;
                else if(!actual && !predicted) valTN++;
                else if(!actual && predicted) valFP++;
                else valFN++;
            }
            
            // Calculate epoch metrics
            float trainErr = 1.0f - (float)epochCorrect / trainSize;
            float trainFPR = epochFP > 0 ? (float)epochFP / (epochFP + epochTN) : 0;
            float trainFNR = epochFN > 0 ? (float)epochFN / (epochFN + epochTP) : 0;
            
            float valErr = 1.0f - (float)valCorrect / valSize;
            float valFPR = valFP > 0 ? (float)valFP / (valFP + valTN) : 0;
            float valFNR = valFN > 0 ? (float)valFN / (valFN + valTP) : 0;
            
            outstream.println(String.format(
                "\nEpoch %d/%d Complete:",
                epoch + 1, epochs));
            outstream.println(String.format(
                "  Train - ERR: %.4f  FPR: %.4f  FNR: %.4f  Loss: %.4f",
                trainErr, trainFPR, trainFNR, epochLoss/trainSize));
            outstream.println(String.format(
                "  Valid - ERR: %.4f  FPR: %.4f  FNR: %.4f  Loss: %.4f",
                valErr, valFPR, valFNR, valLoss/valSize));
            outstream.println();
        }
        
        long totalTime = System.currentTimeMillis() - startTime;
        outstream.println("Training complete! Total time: " + formatTime(totalTime));
    }

    /**
     * Format time in seconds or minutes.
     */
    private String formatTime(long millis) {
        if(millis < 60000) {
            return String.format("%.1fs", millis / 1000.0);
        } else {
            return String.format("%.1fm", millis / 60000.0);
        }
    }

    /**
     * Simple array shuffle using Fisher-Yates algorithm.
     */
    private void shuffleSamples(Sample[] samples) {
        java.util.Random rand = new java.util.Random();
        for(int i = samples.length - 1; i > 0; i--) {
            int j = rand.nextInt(i + 1);
            Sample temp = samples[i];
            samples[i] = samples[j];
            samples[j] = temp;
        }
    }
    /**
     * Backward pass through the entire network.
     */
    private void backward(float[] lossGradient) {
        float[] gradient = lossGradient;
        
        // Backward through layers in reverse order
        for(int i = layers.size() - 1; i >= 0; i--) {
            Layer layer = layers.get(i);
            gradient = layer.backward(gradient);
            
            // ADD: Check for NaN after each layer
            if(checkForNaN("After " + layer.getClass().getSimpleName() + " backward", gradient)) {
                outstream.println("Layer index: " + i);
                outstream.println("Stopping backward pass due to NaN");
                // You might want to skip weight update for this batch
                return;
            }
        }
    }

    /**
     * Update all weights in the network.
     */
    private void updateWeights() {
        for(Layer layer : layers) {
            layer.updateWeights(currentLearningRate);  // Use current learning rate
        }
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

    public static CNNNetwork loadNetwork(String filename) {
        // Implementation depends on your save format
        // This is a placeholder - you'd parse the .bbnet file
        // and reconstruct the network architecture and weights
        try {
            // Load and parse the file
            // Reconstruct layers and weights
            // Return new CNNNetwork instance
            return null; // Placeholder
        } catch(Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    
    // Check for NaN values in a float array
    private boolean checkForNaN(String location, float[] values) {
        for(int i = 0; i < values.length; i++) {
            if(Float.isNaN(values[i]) || Float.isInfinite(values[i])) {
                outstream.println("\n!!! NaN/Inf DETECTED !!!");
                outstream.println("Location: " + location);
                outstream.println("Index: " + i + " of " + values.length);
                outstream.println("Value: " + values[i]);
                
                // Print some surrounding values
                int start = Math.max(0, i - 2);
                int end = Math.min(values.length, i + 3);
                outstream.print("Context values: [");
                for(int j = start; j < end; j++) {
                    if(j == i) outstream.print("*" + values[j] + "*");
                    else outstream.print(values[j]);
                    if(j < end - 1) outstream.print(", ");
                }
                outstream.println("]");
                return true;
            }
        }
        return false;
    }

    /* --------------------------------------------------------------*/
    /*----------------        BBNet Saving           ----------------*/
    /*--------------------------------------------------------------*/
    /**
     * Save the network architecture and weights to a file.
     * @param filename The output file name
     */
    public void saveNetwork(String filename) {
        try (PrintWriter writer = new PrintWriter(filename)) {
            // Header
            writer.println("##bbnet");
            writer.println("#version 1");
            writer.println("#cnn");  // Mark as CNN instead of dense
            writer.println("#layers " + layers.size());
            writer.println("#epochs " + epochs);
            writer.println("#samples " + totalSamplesProcessed);
            
            // Architecture description - use numInputs instead of inputSize
            writer.print("#dims " + numInputs);
            for(Layer layer : layers) {
                writer.print(" " + layer.getOutputSize());
            }
            writer.println();
            
            // Training stats
            writer.println("##err " + String.format("%.6f", lastTrainError));
            writer.println("##fpr " + String.format("%.6f", lastTrainFPR));
            writer.println("##fnr " + String.format("%.6f", lastTrainFNR));
            writer.println();
            
            // Save each layer
            int layerNum = 1;
            for(Layer layer : layers) {
                writer.println("##layer " + layerNum++);
                saveLayerWeights(writer, layer);
                writer.println();
            }
            
        } catch(IOException e) {
            e.printStackTrace();
        }
    }

    private void saveLayerWeights(PrintWriter writer, Layer layer) {
        if(layer instanceof ConvolutionLayer) {
            ConvolutionLayer conv = (ConvolutionLayer) layer;
            // Direct field access since they're package-private
            writer.println("CONV " + conv.inputChannels + " " + conv.numFilters + 
                        " " + conv.filterSize + " " + 1); // stride always 1
            
            // Write filter weights
            for(int f = 0; f < conv.numFilters; f++) {
                writer.print("F" + f + " ");
                // Bias first
                writer.print(conv.bias[f] + " ");
                // Then all filter weights
                for(int c = 0; c < conv.inputChannels; c++) {
                    for(int k = 0; k < conv.filterSize; k++) {
                        writer.print(conv.weights[f][c][k] + " ");
                    }
                }
                writer.println();
            }
            
        } else if(layer instanceof DenseLayer) {
            DenseLayer dense = (DenseLayer) layer;
            // Always ReLU based on your implementation
            String activation = "RELU";
            writer.println("D" + (dense.inputSize + 1) + " " + activation);
            
            // Write weights for each output neuron
            for(int o = 0; o < dense.outputSize; o++) {
                // Bias first
                writer.print(dense.bias[o] + " ");
                // Then all input weights
                for(int i = 0; i < dense.inputSize; i++) {
                    writer.print(dense.weights[o][i] + " ");
                }
                writer.println();
            }
            
        } else if(layer instanceof MaxPoolingLayer) {
            MaxPoolingLayer pool = (MaxPoolingLayer) layer;
            // poolSize is used as both pool size and stride
            writer.println("POOL " + pool.poolSize + " " + pool.poolSize);
        }
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
    private float lastTrainError = 0;
    private float lastTrainFPR = 0; 
    private float lastTrainFNR = 0;
    private int totalSamplesProcessed = 0;
}