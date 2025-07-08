package ml;

import ml.SampleSet;
import ml.DataLoader;
import fileIO.FileFormat;
import shared.Shared;

public class CNNTrainer {
    // --- Only data loading for debugging ---

    /**
     * Loads data using the DataLoader class and prints debug output.
     * @param dataPath Path to the data file
     * @return SampleSet object
     */
    public SampleSet loadData(String dataPath) {
        System.err.println("[DEBUG] Loading data from: " + dataPath);
        SampleSet[] sets = DataLoader.load(dataPath, Shared.MAX_ARRAY_LEN, false, 0f, 0, true, 0f);
        SampleSet set = (sets != null && sets.length > 0) ? sets[0] : null;
        System.err.println("[DEBUG] Loaded SampleSet: " + set);
        // Print more info if possible
        return set;
    }

    // --- CNN network architecture fields ---
    // These define the structure and parameters of the CNN.
    private int numEpochs = 50;         // Number of training epochs
    private int batchSize = 32;         // Batch size for training
    private double learningRate = 0.001;// Learning rate for optimizer
    private int numFilters1 = 32;       // Number of filters in first conv layer
    private int filterSize1 = 7;        // Size of filters in first conv layer
    private int poolSize1 = 2;          // Pooling size after first conv layer
    private int numFilters2 = 64;       // Number of filters in second conv layer
    private int filterSize2 = 5;        // Size of filters in second conv layer
    private int poolSize2 = 2;          // Pooling size after second conv layer
    private int numHidden = 128;        // Number of units in dense hidden layer
    private int numOutputs = 1;         // Number of outputs (binary classification)

    // --- Network weights (parameters) ---
    // These arrays hold the weights for each layer of the CNN.
    private double[][][] convFilters1;   // [numFilters1][1][filterSize1]
    private double[] convBiases1;
    private double[][][] convFilters2;   // [numFilters2][numFilters1][filterSize2]
    private double[] convBiases2;
    private double[][] denseWeights;     // [flattenedSize][numHidden]
    private double[] denseBiases;        // [numHidden]
    private double[][] outWeights;       // [numHidden][numOutputs]
    private double[] outBiases;          // [numOutputs]

    private java.util.Random rand = new java.util.Random(42); // For reproducibility

    /**
     * Constructor: Initializes the CNN weights.
     * This uses He initialization for ReLU layers and Xavier for sigmoid output.
     */
    public CNNTrainer() {
        initializeWeights();
    }

    /**
     * Initializes all network weights with random values.
     * Uses best practices for deep learning initialization.
     */
    private void initializeWeights() {
        // First conv layer: processes raw input
        convFilters1 = new double[numFilters1][1][filterSize1];
        convBiases1 = new double[numFilters1];
        double stdDev1 = Math.sqrt(2.0 / filterSize1); // He initialization
        for (int i = 0; i < numFilters1; i++) {
            for (int j = 0; j < filterSize1; j++) {
                convFilters1[i][0][j] = rand.nextGaussian() * stdDev1;
            }
        }
        // Second conv layer: processes feature maps from first layer
        convFilters2 = new double[numFilters2][numFilters1][filterSize2];
        convBiases2 = new double[numFilters2];
        double stdDev2 = Math.sqrt(2.0 / (numFilters1 * filterSize2));
        for (int i = 0; i < numFilters2; i++) {
            for (int j = 0; j < numFilters1; j++) {
                for (int k = 0; k < filterSize2; k++) {
                    convFilters2[i][j][k] = rand.nextGaussian() * stdDev2;
                }
            }
        }
        // Calculate size after convolutions and pooling (assuming input length 356)
        int sizeAfterConv1 = (356 - filterSize1 + 1) / poolSize1;
        int sizeAfterConv2 = (sizeAfterConv1 - filterSize2 + 1) / poolSize2;
        int flattenedSize = sizeAfterConv2 * numFilters2;
        // Dense layer
        denseWeights = new double[flattenedSize][numHidden];
        denseBiases = new double[numHidden];
        double stdDevDense = Math.sqrt(2.0 / flattenedSize);
        for (int i = 0; i < flattenedSize; i++) {
            for (int j = 0; j < numHidden; j++) {
                denseWeights[i][j] = rand.nextGaussian() * stdDevDense;
            }
        }
        // Output layer (Xavier initialization for sigmoid)
        outWeights = new double[numHidden][numOutputs];
        outBiases = new double[numOutputs];
        double stdDevOut = Math.sqrt(1.0 / numHidden);
        for (int i = 0; i < numHidden; i++) {
            outWeights[i][0] = rand.nextGaussian() * stdDevOut;
        }
    }

    // --- CNN architecture fields (example, adjust as needed) ---
    // Number of input features (length of input vector)
    int inputLength = 100; // TODO: set based on your data
    // Convolution layer parameters
    int numFilters = 8;
    int filterSize = 5;
    float[][] convWeights; // [numFilters][filterSize]
    float[] convBiases;    // [numFilters]
    // Fully connected layer parameters
    int fcOutputSize = 1; // For regression; change for classification
    float[][] fcWeights;   // [numFilters * convOutputLength][fcOutputSize]
    float[] fcBiases;      // [fcOutputSize]

    // --- Forward pass for 1D CNN ---
    /**
     * Performs a forward pass through the CNN for a single input sample.
     * @param input 1D input array (length = inputLength)
     * @return output array (length = fcOutputSize)
     */
    public float[] forward(float[] input) {
        // 1. Convolution layer
        int convOutputLength = input.length - filterSize + 1;
        float[][] convOutput = new float[numFilters][convOutputLength];
        for (int f = 0; f < numFilters; f++) {
            for (int i = 0; i < convOutputLength; i++) {
                float sum = 0f;
                for (int k = 0; k < filterSize; k++) {
                    sum += input[i + k] * convWeights[f][k];
                }
                sum += convBiases[f];
                // ReLU activation
                convOutput[f][i] = Math.max(0, sum);
            }
        }
        // 2. Flatten conv output for fully connected layer
        float[] flat = new float[numFilters * convOutputLength];
        for (int f = 0; f < numFilters; f++) {
            for (int i = 0; i < convOutputLength; i++) {
                flat[f * convOutputLength + i] = convOutput[f][i];
            }
        }
        // 3. Fully connected layer
        float[] output = new float[fcOutputSize];
        for (int o = 0; o < fcOutputSize; o++) {
            float sum = 0f;
            for (int j = 0; j < flat.length; j++) {
                sum += flat[j] * fcWeights[j][o];
            }
            sum += fcBiases[o];
            output[o] = sum; // No activation for regression; use softmax/sigmoid for classification
        }
        return output;
    }

    /**
     * Main method for CLI use. Accepts a data file path as argument.
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java ml.CNNTrainer <data_file>");
            return;
        }
        CNNTrainer trainer = new CNNTrainer();
        SampleSet data = trainer.loadData(args[0]);
        if (data == null || data.samples == null) {
            System.err.println("[ERROR] No data loaded. Exiting.");
            return;
        }
        int nSamples = data.samples.length;
        System.err.println("[DEBUG] Number of samples: " + nSamples);
        System.err.println("[DEBUG] Data loading complete. Ready for training.");
        // ...add training loop here in next steps...
    }
}