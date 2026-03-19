# =============================================================================
# Pipeline: Computer Vision and Machine Learning for Chicken Fecal Analysis
# Author: Susanta Behura
#
# Description:
# This pipeline extracts image-based features from chicken fecal samples and
# applies machine learning to classify disease status (Disease vs Healthy).
#
# Feature categories:
# - Color features (RGB + HSV statistics)
# - Shape features (contour-based morphology)
# - Texture features (Haralick + Local Binary Patterns)
#
# Outputs:
# - Feature dataset (CSV)
# - Classification performance report
# - Feature importance rankings (candidate biomarkers)
#
# Requirements:
# - OpenCV (cv2)
# - NumPy, Pandas
# - scikit-image
# - scikit-learn
# =============================================================================


import cv2
import glob
import numpy as np
import pandas as pd

from skimage.feature import graycomatrix, graycoprops, local_binary_pattern

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report


# =============================================================================
# PARAMETERS
# =============================================================================

LBP_RADIUS = 2
LBP_POINTS = 8 * LBP_RADIUS


# =============================================================================
# FEATURE EXTRACTION FUNCTIONS
# =============================================================================

def extract_color_features(img):
    """
    Extract color statistics from RGB and HSV color spaces.
    """
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)

    features = {}

    # RGB statistics
    for i, name in enumerate(["B", "G", "R"]):
        channel = img[:, :, i]
        features[f"{name}_mean"] = np.mean(channel)
        features[f"{name}_std"] = np.std(channel)

    # HSV statistics
    for i, name in enumerate(["H", "S", "V"]):
        channel = hsv[:, :, i]
        features[f"{name}_mean"] = np.mean(channel)
        features[f"{name}_std"] = np.std(channel)

    return features


def segment_sample(gray):
    """
    Segment fecal region using adaptive thresholding and morphology.
    """
    blur = cv2.GaussianBlur(gray, (5, 5), 0)
    equalized = cv2.equalizeHist(blur)

    thresh = cv2.adaptiveThreshold(
        equalized,
        255,
        cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
        cv2.THRESH_BINARY_INV,
        51,
        2
    )

    kernel = np.ones((5, 5), np.uint8)
    mask = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)

    return mask


def extract_shape_features(mask):
    """
    Extract morphological features from segmented object.
    """
    contours, _ = cv2.findContours(
        mask,
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_SIMPLE
    )

    if len(contours) == 0:
        return {}

    c = max(contours, key=cv2.contourArea)

    area = cv2.contourArea(c)
    perimeter = cv2.arcLength(c, True)

    x, y, w, h = cv2.boundingRect(c)
    hull = cv2.convexHull(c)
    convex_area = cv2.contourArea(hull)

    features = {
        "area": area,
        "perimeter": perimeter,
        "aspect_ratio": w / h if h > 0 else 0,
        "solidity": area / convex_area if convex_area > 0 else 0,
        "circularity": (4 * np.pi * area) / (perimeter ** 2) if perimeter > 0 else 0
    }

    return features


def extract_haralick(gray):
    """
    Extract Haralick texture features using GLCM.
    """
    glcm = graycomatrix(
        gray,
        distances=[1],
        angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
        levels=256,
        symmetric=True,
        normed=True
    )

    features = {}
    props = ["contrast", "dissimilarity", "homogeneity", "energy", "correlation", "ASM"]

    for p in props:
        features[f"haralick_{p}"] = graycoprops(glcm, p).mean()

    return features


def extract_lbp(gray):
    """
    Extract Local Binary Pattern (LBP) histogram features.
    """
    lbp = local_binary_pattern(gray, LBP_POINTS, LBP_RADIUS, method="uniform")

    hist, _ = np.histogram(
        lbp.ravel(),
        bins=np.arange(0, LBP_POINTS + 3),
        range=(0, LBP_POINTS + 2)
    )

    hist = hist.astype("float")
    hist /= hist.sum() + 1e-6  # Normalize histogram

    features = {f"lbp_{i}": val for i, val in enumerate(hist)}

    return features


# =============================================================================
# LOAD AND PROCESS IMAGE DATA
# =============================================================================

data = []

datasets = [
    ("Disease", 1),
    ("Healthy", 0)
]

for folder, label in datasets:
    files = glob.glob(f"{folder}/*.jpg")

    for path in files:
        img = cv2.imread(path)

        if img is None:
            continue

        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        mask = segment_sample(gray)

        features = {
            "image": path,
            "label": label
        }

        # Combine all feature types
        features.update(extract_color_features(img))
        features.update(extract_shape_features(mask))
        features.update(extract_haralick(gray))
        features.update(extract_lbp(gray))

        data.append(features)


# =============================================================================
# CREATE FEATURE DATASET
# =============================================================================

df = pd.DataFrame(data)

df.to_csv("chicken_fecal_features_dataset.csv", index=False)

print("Dataset shape:", df.shape)


# =============================================================================
# MACHINE LEARNING CLASSIFICATION
# =============================================================================

X = df.drop(columns=["image", "label"])
y = df["label"]

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, random_state=42
)

# Train Random Forest classifier
model = RandomForestClassifier(
    n_estimators=300,
    random_state=42
)

model.fit(X_train, y_train)

# Evaluate model
pred = model.predict(X_test)

print("\nClassification Report:")
print(classification_report(y_test, pred))


# =============================================================================
# BIOMARKER DISCOVERY (FEATURE IMPORTANCE)
# =============================================================================

importance = pd.Series(
    model.feature_importances_,
    index=X.columns
).sort_values(ascending=False)

importance.to_csv("feature_importance_biomarkers.csv")

print("\nTop Potential Biomarkers:")
print(importance.head(20))