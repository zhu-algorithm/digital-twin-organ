from pathlib import Path
import pandas as pd
from tdln_pd1.model import train_evaluate


def test_train_evaluate_writes_reproducible_artifacts(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4],
        "frac_cd8_t": [0.2, 0.3, 0.2, 0.3, 0.7, 0.8, 0.7, 0.8],
    })
    labels = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)], "response": [0]*4 + [1]*4})
    metrics = train_evaluate(features, labels, tmp_path, folds=2)
    assert metrics["n_patients"] == 8
    assert (tmp_path / "model.joblib").exists()
    assert (tmp_path / "cv_predictions.csv").exists()
