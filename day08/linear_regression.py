"""
Linear Regression with Gradient Descent
Complete function library extracted from the original notebook implementation
Includes: preprocessing, gradient descent, pseudoinverse, feature selection
"""

import numpy as np
import itertools


# ============================================================================
# PREPROCESSING FUNCTIONS
# ============================================================================

def preprocess(X, y):
    """
    Perform min-max scaling for both the data and the targets.
    
    Input:
    - X: Inputs (n features, m instances).
    - y: True labels (1 target, m instances).

    Output:
    - X: The scaled inputs.
    - y: The scaled labels.
    """
    X_min = X.min(axis=0)
    Y_min = y.min()
    X_max = X.max(axis=0)
    Y_max = y.max()
    # Avoid division by zero
    X_max[X_max == X_min] = X_min[X_max == X_min] + 1
    if Y_max == Y_min:
        Y_max = Y_min + 1    
    X = (X - X_min) / (X_max - X_min)
    y = (y - Y_min) / (Y_max - Y_min)
    
    return X, y


def reverse_scaling(y_scaled, y_original):
    """
    Reverse the min-max scaling to get original values.
    
    Input:
    - y_scaled: Scaled target values.
    - y_original: Original target values (before scaling).
    
    Output:
    - y: Values in original scale.
    """
    y_min = y_original.min()
    y_max = y_original.max()
    y = y_scaled * (y_max - y_min) + y_min
    return y


# ============================================================================
# CORE LINEAR REGRESSION FUNCTIONS
# ============================================================================

def compute_cost(X, y, theta):
    """
    Computes the MSE cost for linear regression.

    Input:
    - X: Inputs (n features over m instances).
    - y: True labels (1 value over m instances).
    - theta: The parameters (weights) of the model being learned.

    Output:
    - J: the cost associated with the current set of parameters (single number).
    """
    m = X.shape[0]
    predictions = X @ theta
    squared_errors = (predictions - y) ** 2
    J = (1 / (2 * m)) * np.sum(squared_errors)
    return J


def gradient_descent(X, y, theta, alpha, num_iters):
    """
    Learn the parameters of the model using gradient descent with early stopping.
    
    Input:
    - X: Inputs (n features over m instances).
    - y: True labels (1 value over m instances).
    - theta: The parameters (weights) of the model being learned.
    - alpha: The learning rate of the model.
    - num_iters: The number of iterations performed.

    Output:
    - theta: The learned parameters of the model.
    - J_history: the loss value in each iteration.
    """
    J_history = []
    m = X.shape[0]
    XT = X.T

    for i in range(num_iters):
        predictions = X.dot(theta)
        errors = predictions - y
        gradient = (1 / m) * (XT.dot(errors))
        theta = theta - alpha * gradient
        J = compute_cost(X, y, theta)
        J_history.append(J)
        
        # Early stopping: check if improvement is less than 1e-8
        if i > 0 and abs(J_history[i-1] - J_history[i]) < 1e-8:
            break

        if np.isnan(J) or np.isinf(J):
            print("Gradient descent diverged. Try reducing alpha.")
            break
    
    return theta, J_history


def pinv(X, y):
    """
    Calculate the optimal parameters using pseudoinverse (normal equation).

    Input:
    - X: Inputs (n features over m instances).
    - y: True labels (1 value over m instances).

    Output:
    - theta: The optimal parameters of your model.
    """
    pinv_theta = np.linalg.inv(X.T @ X) @ X.T @ y
    return pinv_theta


# ============================================================================
# HYPERPARAMETER TUNING
# ============================================================================

def find_best_alpha(X, y, iterations):
    """
    Iterate over predefined alpha values and find the best learning rate.
    
    Input:
    - X: Inputs (n features over m instances).
    - y: True labels (1 value over m instances).
    - iterations: The number of iterations performed.

    Output:
    - alpha_dict: Dictionary with alpha as key and final loss as value.
    """
    alphas = [0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 2, 3]
    alpha_dict = {}
    
    for alpha in alphas:
        np.random.seed(42)
        theta = np.random.random(size=X.shape[1])
        final_theta, J_history = gradient_descent(X, y, theta, alpha, iterations)
        final_cost = J_history[-1]
        alpha_dict[alpha] = final_cost
    
    return alpha_dict


# ============================================================================
# FEATURE SELECTION
# ============================================================================

def generate_triplets(X):
    """
    Generate all possible sets of three features from dataset X.
    
    Input:
    - X: Inputs (n features over m instances, including bias at column 0).

    Output:
    - triplets: List of feature triplets as tuples (excluding bias).
    """
    # Generate triplets from actual features only (exclude bias column at index 0)
    n_features = X.shape[1] - 1
    triplets = list(itertools.combinations(range(n_features), 3))
    return triplets


def find_best_triplet(X_train, y_train, X_val, y_val, triplets, alpha, num_iter):
    """
    Find the best triplet of features that minimizes validation cost.

    Input:
    - X_train: Training dataset (preprocessed with bias column at index 0).
    - y_train: Training labels (preprocessed).
    - X_val: Validation dataset (preprocessed with bias column at index 0).
    - y_val: Validation labels (preprocessed).
    - triplets: List of three-feature combinations.
    - alpha: Learning rate.
    - num_iter: Number of gradient descent iterations.

    Output:
    - best_triplet: The best triplet of features.
    """
    best_triplet = None
    best_val_cost = float('inf')
    
    # Initialize theta (same for all triplets for consistency)
    np.random.seed(42)
    theta = np.random.random(size=4)  # 3 features + bias term = 4
    
    for triplet in triplets:
        # Extract the triplet features + bias term
        X_train_triplet = X_train[:, [0] + [i + 1 for i in triplet]]
        X_val_triplet = X_val[:, [0] + [i + 1 for i in triplet]]
        
        # Train the model
        final_theta, J_history = gradient_descent(X_train_triplet, y_train, theta, alpha, num_iter)
        
        # Compute validation cost
        val_cost = compute_cost(X_val_triplet, y_val, final_theta)
        
        # Update best triplet if needed
        if best_triplet is None or val_cost < best_val_cost:
            best_triplet = triplet
            best_val_cost = val_cost
    
    return best_triplet


def forward_feature_selection(X_train, y_train, X_val, y_val, alpha, num_iter, max_features=3):
    """
    Perform greedy forward feature selection.
    Iteratively add features that minimize validation cost.
    
    Input:
    - X_train: Training dataset (with bias column at index 0).
    - y_train: Training labels.
    - X_val: Validation dataset (with bias column at index 0).
    - y_val: Validation labels.
    - alpha: Learning rate.
    - num_iter: Number of gradient descent iterations.
    - max_features: Maximum number of features to select (default 3).
    
    Output:
    - selected_features: List of selected feature indices (excluding bias).
    - validation_costs: List of validation costs at each step.
    """
    # Number of features (excluding bias column at index 0)
    n_features = X_train.shape[1] - 1
    
    # Track selected features (indices without bias)
    selected_features = []
    
    # Track remaining features to consider
    remaining_features = list(range(n_features))
    
    # Track validation costs at each step
    validation_costs = []
    
    # Iteratively select features until we reach max_features
    for step in range(max_features):
        best_feature = None
        best_val_cost = float('inf')
        
        # Try adding each remaining feature
        for feature_idx in remaining_features:
            # Create feature set: bias + selected features + candidate feature
            current_feature_indices = [0] + [f + 1 for f in selected_features] + [feature_idx + 1]
            
            # Extract columns
            X_train_subset = X_train[:, current_feature_indices]
            X_val_subset = X_val[:, current_feature_indices]
            
            # Initialize theta
            np.random.seed(42)
            theta = np.random.random(size=len(current_feature_indices))
            
            # Train model
            final_theta, J_history = gradient_descent(X_train_subset, y_train, theta, alpha, num_iter)
            
            # Compute validation cost
            val_cost = compute_cost(X_val_subset, y_val, final_theta)
            
            # Update best feature
            if val_cost < best_val_cost:
                best_val_cost = val_cost
                best_feature = feature_idx
        
        # Add best feature to selected set
        selected_features.append(best_feature)
        remaining_features.remove(best_feature)
        validation_costs.append(best_val_cost)
    
    return selected_features, validation_costs
