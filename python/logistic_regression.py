import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import sys

def run_baseline_test(filepath):
    """
    Loads a TSV training set, trains a logistic regression model,
    and reports its performance as a simple baseline.
    """
    print(f"--- Loading data from: {filepath} ---")
    try:
        data = pd.read_csv(filepath, sep='\t')
    except FileNotFoundError:
        print(f"FATAL ERROR: File not found at '{filepath}'. Aborting.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        print("Please ensure it is a tab-separated file with a header.")
        sys.exit(1)

    # For this baseline, let's treat all non-positives as a single negative class.
    # We'll map anything that isn't '1' to '0'.
    data.iloc[:, -1] = data.iloc[:, -1].apply(lambda x: 1 if x == 1 else 0)

    X = data.iloc[:, :-1]
    y = data.iloc[:, -1]
    
    if len(y.unique()) < 2:
        print("ERROR: The dataset contains only one class. Cannot train a classifier.")
        sys.exit(1)

    print(f"Loaded {len(data)} samples.")
    print(f"Class distribution:\n{y.value_counts()}\n")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    print("--- Training Logistic Regression Model ---")
    model = LogisticRegression(C=0.1, max_iter=1000, solver='liblinear')
    model.fit(X_train, y_train)
    print("Model training complete.\n")

    print("--- Evaluating Model Performance on Test Set ---")
    y_pred = model.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)
    print(f"Overall Accuracy: {accuracy:.4f}")
    
    # --- NEW: Explicit Calculation and Printout ---
    cm = confusion_matrix(y_test, y_pred)
    tn, fp, fn, tp = cm.ravel()

    print("\nConfusion Matrix Breakdown:")
    print(f"  - True Positives (Gene -> Gene):       {tp}")
    print(f"  - True Negatives (Not Gene -> Not Gene): {tn}")
    print(f"  - False Positives (Not Gene -> Gene):    {fp}  <-- Model is hallucinating genes")
    print(f"  - False Negatives (Gene -> Not Gene):    {fn}  <-- Model is missing real genes")
    # --- END OF NEW BLOCK ---

    print("\nFull Classification Report:")
    print(classification_report(y_test, y_pred, target_names=['Not a Gene (0)', 'Gene (1)']))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python logistic_regression.py <path_to_your_training_set.tsv>")
        sys.exit(1)
    
    training_file = sys.argv[1]
    run_baseline_test(training_file)
