# Neural Network Architecture Analysis

This document details the architecture of the Feed-Forward Neural Network (FFNN) and Convolutional Neural Network (CNN) found in the `ml/` directory.

## Feed-Forward Neural Network (FFNN)

The FFNN is implemented in the `ml/CellNet.java` file. It is a flexible, layer-based network that can be configured as either dense or sparse.

### Architecture

*   **Layers**: The network is composed of an input layer, one or more hidden layers, and an output layer. The dimensions of these layers are specified by the `dims` array in the constructor ([`ml/CellNet.java:38`](ml/CellNet.java:38)).
*   **Connectivity**: The network can be either fully-connected (dense) or sparsely-connected. This is controlled by the `DENSE` static boolean flag ([`ml/CellNet.java:1118`](ml/CellNet.java:1118)).
    *   In **dense mode**, every neuron in a layer is connected to every neuron in the previous layer ([`ml/CellNet.java:231`](ml/CellNet.java:231)).
    *   In **sparse mode**, connections are created based on a `density` parameter, which controls the probability of a connection forming ([`ml/CellNet.java:275`](ml/CellNet.java:275)).
*   **Weight Initialization**: Weights are initialized randomly. The `randomWeight` method uses a combination of uniform and quadratic distributions to generate initial weight values ([`ml/CellNet.java:385`](ml/CellNet.java:385)). Biases are also initialized randomly ([`ml/CellNet.java:372`](ml/CellNet.java:372)).

### Activation Functions

*   The default activation function for hidden layers is configurable. The `Cell.defaultActivationType` variable determines the default function ([`ml/CellNet.java:211`](ml/CellNet.java:211)).
*   The final layer has its own activation type, `Cell.finalLayerType` ([`ml/CellNet.java:211`](ml/CellNet.java:211)).
*   The network supports a variety of activation functions, which are defined in `ml/Functions.java`. These include:
    *   **Sigmoid**: ([`ml/Functions.java:9`](ml/Functions.java:9))
    *   **Extended Sigmoid**: ([`ml/Functions.java:21`](ml/Functions.java:21))
    *   **TanH**: ([`ml/Functions.java:46`](ml/Functions.java:46))
    *   **Swish**: ([`ml/Functions.java:76`](ml/Functions.java:76))
    *   **RSLog** (Root-Signed Logarithm): ([`ml/Functions.java:109`](ml/Functions.java:109))
    *   **Mirrored Sigmoid**: ([`ml/Functions.java:139`](ml/Functions.java:139))
    *   **Gaussian (Bell curve)**: ([`ml/Functions.java:199`](ml/Functions.java:199))
*   Activation functions for hidden layers can be randomized during network initialization if configured ([`ml/CellNet.java:145`](ml/CellNet.java:145)).

### Training Strategy

*   **Backpropagation**: The network is trained using backpropagation. Separate methods exist for dense and sparse networks: `backPropDense` ([`ml/CellNet.java:532`](ml/CellNet.java:532)) and `backPropSparse` ([`ml/CellNet.java:505`](ml/CellNet.java:505)).
*   **Loss Function**: The specific loss function is not explicitly defined in a separate class like the CNN, but the error calculation in the `calcError` method ([`ml/CellNet.java:696`](ml/CellNet.java:696)) suggests a form of squared error.
*   **Weight Updates**: Weight updates are applied in the `applyChanges` method ([`ml/CellNet.java:668`](ml/CellNet.java:668)). The updates are scaled by the number of samples in the batch and a learning rate (`alpha`).

## Convolutional Neural Network (CNN)

The CNN is implemented across several files, with the main network structure in `ml/CNNNetwork.java` and the training process orchestrated by `ml/CNNTrainer.java`.

### Architecture

The CNN architecture is configurable through command-line arguments parsed in `ml/CNNTrainer.java`. The network consists of a sequence of convolutional, pooling, and dense layers.

*   **Convolutional Layers** (`ml/ConvolutionLayer.java`):
    *   **1D Convolutions**: The network uses 1D convolutions, suitable for sequence data ([`ml/ConvolutionLayer.java:8`](ml/ConvolutionLayer.java:8)).
    *   **Kernel Size**: The size of the convolutional filters (kernels) is configurable via the `filtersizes` parameter ([`ml/CNNTrainer.java:118`](ml/CNNTrainer.java:118)).
    *   **Number of Filters**: The number of filters in each convolutional layer is set by the `filters` parameter ([`ml/CNNTrainer.java:110`](ml/CNNTrainer.java:110)).
    *   **Activation Function**: Convolutional layers use the **ReLU** (Rectified Linear Unit) activation function (`Math.max(0, sum)`) ([`ml/ConvolutionLayer.java:65`](ml/ConvolutionLayer.java:65)).
    *   **Weight Initialization**: Weights are initialized using Xavier initialization to maintain variance across layers ([`ml/ConvolutionLayer.java:27`](ml/ConvolutionLayer.java:27)).

*   **Max Pooling Layers** (`ml/MaxPoolingLayer.java`):
    *   **Pooling Size**: The size of the max pooling window is configurable via the `poolsizes` parameter ([`ml/CNNTrainer.java:126`](ml/CNNTrainer.java:126)). The stride is equal to the pool size.

*   **Dense Layers** (`ml/DenseLayer.java`):
    *   **Fully Connected**: After the convolutional and pooling layers, the flattened output is passed through one or more dense layers. The size of these layers is configured by the `dense` parameter ([`ml/CNNTrainer.java:134`](ml/CNNTrainer.java:134)).
    *   **Activation Function**: Dense layers also use the **ReLU** activation function ([`ml/DenseLayer.java:47`](ml/DenseLayer.java:47)).
    *   **Weight Initialization**: Dense layer weights are also initialized using Xavier initialization ([`ml/DenseLayer.java:24`](ml/DenseLayer.java:24)).

*   **Output Layer** (`ml/OutputLayer.java`):
    *   The final layer is an `OutputLayer` that applies a **Sigmoid** activation function, which is appropriate for binary classification tasks, to produce a probability score between 0 and 1.

### Training Strategy

The training process is managed by `ml/CNNTrainer.java` and `ml/CNNNetwork.java`.

*   **Loss Function**: The network uses **Binary Cross-Entropy** as its loss function, which is standard for binary classification problems. The implementation is in `ml/BinaryCrossEntropyLoss.java` ([`ml/BinaryCrossEntropyLoss.java:5`](ml/BinaryCrossEntropyLoss.java:5)).
*   **Optimization**:
    *   **Mini-batch Gradient Descent**: The network is trained using mini-batches of data. The `batchSize` is configurable ([`ml/CNNTrainer.java:146`](ml/CNNTrainer.java:146)).
    *   **Learning Rate**: The learning rate is configurable (`learningRate`) and is decayed over time, multiplied by `0.95` each epoch ([`ml/CNNNetwork.java:137`](ml/CNNNetwork.java:137)).
    *   **Gradient Clipping**: Gradients in the dense and convolutional layers are clipped to a range of `[-5.0, 5.0]` to prevent exploding gradients during backpropagation ([`ml/DenseLayer.java:92`](ml/DenseLayer.java:92), [`ml/ConvolutionLayer.java:123`](ml/ConvolutionLayer.java:123)).
*   **Regularization**:
    *   **Dropout**: Dropout is used to prevent overfitting, with the dropout rate being a configurable parameter ([`ml/CNNTrainer.java:150`](ml/CNNTrainer.java:150)). *Note: While the `dropout` parameter is parsed, its application in the layer code is not immediately apparent from the provided files.*
*   **Evaluation**: The training process includes evaluation on a validation set after each epoch, tracking metrics like error rate, FPR, FNR, and loss ([`ml/CNNNetwork.java:309`](ml/CNNNetwork.java:309)). The `CNNTrainer` also includes a separate evaluation mode to test a trained network ([`ml/CNNTrainer.java:392`](ml/CNNTrainer.java:392)).