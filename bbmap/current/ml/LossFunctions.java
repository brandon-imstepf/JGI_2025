package ml;

import shared.Tools;

/** Brandon Imstepf
 * Last Modified: November 24, 2025
 * Centralizes loss calculations for CellNet training so new objectives can be wired
 * in without rewriting the trainer. Loss gradients are computed with respect to the
 * network output (post-activation) so the existing backprop plumbing can remain intact.
 */
public final class LossFunctions {

	public enum LossType {
		MSE,
		WEIGHTED_BCE,
		TVERSKY
	}

	private static LossType lossType=LossType.MSE;

	private static float posWeight=1f;
	private static float negWeight=1f;
	private static float tverskyAlpha=0.3f;
	private static float tverskyBeta=0.7f;
	private static float epsilon=1e-6f;
	private static float tverskySmooth=1e-6f;
	private static boolean legacyWeighting=true;

	public static void setLossType(LossType type) {
		lossType=type;
		legacyWeighting=(type==LossType.MSE);
	}

	public static void setLossType(String name) {
		if(name==null) {return;}
		switch(name.toLowerCase()) {
		case "mse":
		case "ls":
		case "l2":
			setLossType(LossType.MSE);
			break;
		case "wbce":
		case "weighted_bce":
		case "bce":
		case "crossentropy":
			setLossType(LossType.WEIGHTED_BCE);
			break;
		case "tversky":
			setLossType(LossType.TVERSKY);
			break;
		default:
			throw new RuntimeException("Unknown loss type: "+name);
		}
	}

	public static LossType lossType() {return lossType;}
	public static boolean useLegacyWeighting() {return legacyWeighting;}
	public static void setLegacyWeighting(boolean useLegacy) {legacyWeighting=useLegacy;}

	public static void setPositiveWeight(float w) {posWeight=Tools.max(1e-6f, w);}
	public static void setNegativeWeight(float w) {negWeight=Tools.max(1e-6f, w);}
	public static void setEpsilon(float e) {epsilon=Tools.max(1e-12f, e);}
	public static void setTverskyAlpha(float a) {tverskyAlpha=Tools.max(0f, a);}
	public static void setTverskyBeta(float b) {tverskyBeta=Tools.max(0f, b);}
	public static void setTverskySmooth(float s) {tverskySmooth=Tools.max(1e-12f, s);}

	public static float loss(float target, float pred) {
		switch(lossType) {
		case WEIGHTED_BCE:
			return weightedBCELoss(target, pred);
		case TVERSKY:
			return tverskyLoss(target, pred);
		case MSE:
		default:
			return mseLoss(target, pred);
		}
	}

	public static float gradient(float pred, float target) {
		switch(lossType) {
		case WEIGHTED_BCE:
			return weightedBCEGradient(target, pred);
		case TVERSKY:
			return tverskyGradient(target, pred);
		case MSE:
		default:
			return mseGradient(target, pred);
		}
	}

	/* ---------------- Loss implementations ---------------- */

	private static float mseLoss(float target, float pred) {
		final float diff=pred-target;
		return 0.5f*diff*diff;
	}

	private static float mseGradient(float target, float pred) {
		return pred-target;
	}

	private static float weightedBCELoss(float target, float pred) {
		final float p=clampProbability(pred);
		final float lossPos=-posWeight*target*(float)Math.log(p);
		final float lossNeg=-negWeight*(1f-target)*(float)Math.log(1f-p);
		return lossPos+lossNeg;
	}

	private static float weightedBCEGradient(float target, float pred) {
		final float p=clampProbability(pred);
		final float gradPos=-posWeight*target/p;
		final float gradNeg=negWeight*(1f-target)/(1f-p);
		return gradPos+gradNeg;
	}

	private static float tverskyLoss(float target, float pred) {
		final float p=clampProbability(pred);
		final float tp=target*p;
		final float fp=(1f-target)*p;
		final float fn=target*(1f-p);
		final float numerator=tp+tverskySmooth;
		final float denominator=tp+tverskyAlpha*fp+tverskyBeta*fn+tverskySmooth;
		return 1f-(numerator/denominator);
	}

	private static float tverskyGradient(float target, float pred) {
		final float p=clampProbability(pred);
		final float tp=target*p;
		final float fp=(1f-target)*p;
		final float fn=target*(1f-p);
		final float numerator=tp+tverskySmooth;
		final float denominator=tp+tverskyAlpha*fp+tverskyBeta*fn+tverskySmooth;

		float dNumerator=target;
		float dDenominator;
		if(target>=0.5f) {
			dDenominator=1f-tverskyBeta;
		}else {
			dDenominator=tverskyAlpha;
		}
		final float dFraction=(dNumerator*denominator - numerator*dDenominator)/(denominator*denominator);
		return -dFraction;
	}

	private static float clampProbability(float pred) {
		return Tools.mid(pred, epsilon, 1f-epsilon);
	}

	private LossFunctions() {}
}
