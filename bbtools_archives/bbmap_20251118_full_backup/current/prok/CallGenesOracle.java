package prok;

import java.util.Locale;
import prok.CallGenesHelper.GeneQuad;

/**
 * Oracle utility for CallGenes that can operate in multiple advisory modes:
 * <ul>
 *     <li>NEURAL – use the neural network output (current behavior).</li>
 *     <li>NEUTRAL – emit the cutoff value so the multiplier becomes 1.0.</li>
     *    <li>TRUTH – emit perfect binary scores derived from the reference GFF.</li>
 * </ul>
 */
public class CallGenesOracle {

    public enum Mode {
        NEURAL,
        NEUTRAL,
        TRUTH;

        public static Mode fromString(String value) {
            if (value == null) {
                return NEURAL;
            }
            return switch (value.trim().toLowerCase(Locale.ROOT)) {
                case "neutral", "none", "off" -> NEUTRAL;
                case "truth", "oracle", "perfect" -> TRUTH;
                case "nn", "neural", "default" -> NEURAL;
                default -> throw new IllegalArgumentException("Unknown oracle mode: " + value);
            };
        }
    }

    public interface TruthLookup {
        boolean isTrueGene(GeneQuad quad);
    }

    private Mode mode = Mode.NEURAL;
    private float neutralScore = 0.5f;
    private TruthLookup truthLookup = null;

    public Mode getMode() {
        return mode;
    }

    public void setMode(Mode mode) {
        this.mode = (mode == null) ? Mode.NEURAL : mode;
    }

    public void setMode(String value) {
        this.mode = Mode.fromString(value);
    }

    public void setNeutralScore(float nnCutoff) {
        this.neutralScore = nnCutoff;
    }

    public void setTruthLookup(TruthLookup lookup) {
        this.truthLookup = lookup;
    }

    public boolean needsNeuralScore() {
        return mode == Mode.NEURAL;
    }

    public float scoreNeutral() {
        return neutralScore;
    }

    public float scoreNeural(float nnScore) {
        return (mode == Mode.NEUTRAL) ? neutralScore : nnScore;
    }

    public float scoreTruth(CallGenesHelper.GeneQuad quad) {
        if (mode != Mode.TRUTH) {
            return neutralScore;
        }
        if (truthLookup == null || quad == null) {
            return neutralScore;
        }
        return truthLookup.isTrueGene(quad) ? 1.0f : 0.0f;
    }
}
