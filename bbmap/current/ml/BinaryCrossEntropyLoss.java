package ml;

/**
 * Binary cross-entropy loss for binary classification.
 * Loss = -[y*log(p) + (1-y)*log(1-p)]
 */
public class BinaryCrossEntropyLoss {
    
    /**
     * Calculate the loss for a single prediction.
     * @param predicted The predicted probability (0-1)
     * @param target The true label (0 or 1)
     * @return The loss value
     */
    public static float calculateLoss(float prediction, float target) {
        // Clip predictions to avoid log(0)
        float epsilon = 1e-7f;
        float clippedPred = Math.max(epsilon, Math.min(1.0f - epsilon, prediction));
        
        // Binary cross-entropy: -[y*log(p) + (1-y)*log(1-p)]
        return -(target * (float)Math.log(clippedPred) + 
                (1 - target) * (float)Math.log(1 - clippedPred));
    }
    
    /**
     * Calculate the gradient of the loss with respect to the prediction.
     * @param predicted The predicted probability (0-1)
     * @param target The true label (0 or 1)
     * @return The gradient
     */
    public static float calculateGradient(float prediction, float target) {
        // Clip to avoid division by zero
        float epsilon = 1e-7f;
        float clippedPred = Math.max(epsilon, Math.min(1.0f - epsilon, prediction));
        
        // Gradient: (p - y) / (p * (1 - p))
        return (clippedPred - target) / (clippedPred * (1 - clippedPred));
    }
}