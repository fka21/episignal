# Import required libraries
import pandas as pd
from pycaret.classification import *

# Load the dataset
data = pd.read_csv(
    "/Users/ferenc.kagan/Documents/Projects/episignal/input/6_samples_tpm_human.csv",
    delimiter="\t",
)

# Define the vector of normal sample names
normal_samples = ["CG-in_S31", "Nav-in_S36", "QB-in_S33"]

# Create the sample information DataFrame
sample_info = pd.DataFrame(
    {
        "sampleId": data.columns[1:],  # Skip the first column (gene IDs)
        "status": [
            "normal" if sample in normal_samples else "inversed"
            for sample in data.columns[1:]
        ],
    }
)


# Add the sample status to the dataset for clustering
data_transformed = data.set_index("Geneid").T  # Transpose to have samples as rows
data_transformed["status"] = sample_info["status"].values

# Initialize PyCaret setup for classification with feature selection and PCA
setup_instance = setup(
    data=data_transformed,
    target=data_transformed[
        "status"
    ],  # Binary target with 'normal' and 'inversed' classes
    session_id=42,  # Random seed for reproducibility
    feature_selection=True,  # Reduces number of features to avoid overfitting
    normalize=True,  # Scales data to improve model training
    transformation=True,  # Normalizes skewed distributions (Box-Cox/Quantile)
    ignore_features=["status"],  # Exclude the target variable from features
    numeric_features=data_transformed.columns.drop("status").tolist(),
    data_split_shuffle=True,  # Shuffle to ensure random splitting
    data_split_stratify=True,  # Stratify split by the 'status' classes
    fold_strategy="kfold",
    fold=2,
    remove_multicollinearity=True,  # Remove redundant correlated features
    remove_outliers=False,  # Avoid removing data points due to small sample size
    html=False,  # Non-interactive output mode
)

best = compare_models()

plot_model(best)
plot_model(best, plot="confusion_matrix")
plot_model(best, plot="feature")
