#!/usr/bin/env python

# Test the learn.auto_ml module
#
# h2o needs a running Java server, so it is not installed for the tests.
# Instead a small in-memory fake of the parts of ``h2o`` and ``h2o.automl``
# used by teemi is injected into ``sys.modules`` before the module is imported.

import importlib
import re
import sys
import types

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import teemi.learn


class FakeColumn:
    def __init__(self, name, is_factor=False):
        self.name = name
        self.is_factor = is_factor

    def asfactor(self):
        return FakeColumn(self.name, is_factor=True)


class FakeH2OFrame:
    """Mimics h2o.H2OFrame: built from a pandas DataFrame, string column names,
    ``len()`` and column re-assignment (used for ``asfactor``)."""

    created = []

    def __init__(self, python_obj):
        self.data = python_obj.copy()
        self.factor_columns = []
        FakeH2OFrame.created.append(self)

    @property
    def columns(self):
        # H2O always exposes column names as strings
        return [str(c) for c in self.data.columns]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, column):
        return FakeColumn(column)

    def __setitem__(self, column, value):
        if value.is_factor:
            self.factor_columns.append(column)


class FakeTwoDimTable:
    def __init__(self, df):
        self.df = df

    def as_data_frame(self):
        return self.df


class FakeModel:
    """Best model whose metrics are a deterministic function of the number of
    rows it was trained on, so the results can be traced back."""

    def __init__(self, n_rows):
        self.n_rows = n_rows
        self.model_id = f"GBM_model_{n_rows}_rows"

    def mae(self):
        return 10.0 / self.n_rows

    def cross_validation_metrics_summary(self):
        # Same layout as H2O: first row is MAE, first columns '', mean, sd
        return FakeTwoDimTable(
            pd.DataFrame(
                {
                    "": ["mae", "mean_residual_deviance", "mse"],
                    "mean": [20.0 / self.n_rows, 99.0, 99.0],
                    "sd": [1.0 / self.n_rows, 99.0, 99.0],
                    "cv_1_valid": [0.0, 0.0, 0.0],
                    "cv_2_valid": [0.0, 0.0, 0.0],
                }
            )
        )


class StopAfterTraining(Exception):
    """Raised by the fake to stop the run once the training phase is done."""


class FakeH2OAutoML:
    instances = []
    stop_after_training = False

    def __init__(self, **kwargs):
        self.init_kwargs = kwargs
        self.train_kwargs = None
        FakeH2OAutoML.instances.append(self)

    def train(self, x=None, y=None, training_frame=None):
        self.train_kwargs = {"x": x, "y": y, "training_frame": training_frame}

    def get_best_model(self):
        if FakeH2OAutoML.stop_after_training:
            raise StopAfterTraining
        return FakeModel(len(self.train_kwargs["training_frame"]))


@pytest.fixture
def auto_ml(monkeypatch):
    """Import teemi.learn.auto_ml against the fake h2o and clean up afterwards."""
    FakeH2OFrame.created = []
    FakeH2OAutoML.instances = []
    FakeH2OAutoML.stop_after_training = False

    fake_h2o = types.ModuleType("h2o")
    fake_automl = types.ModuleType("h2o.automl")
    fake_h2o.H2OFrame = FakeH2OFrame
    fake_automl.H2OAutoML = FakeH2OAutoML
    fake_h2o.automl = fake_automl
    monkeypatch.setitem(sys.modules, "h2o", fake_h2o)
    monkeypatch.setitem(sys.modules, "h2o.automl", fake_automl)
    sys.modules.pop("teemi.learn.auto_ml", None)

    module = importlib.import_module("teemi.learn.auto_ml")
    yield module

    sys.modules.pop("teemi.learn.auto_ml", None)
    if getattr(teemi.learn, "auto_ml", None) is module:
        delattr(teemi.learn, "auto_ml")


@pytest.fixture
def ml_df():
    # Target column first, categorical features after it
    return pd.DataFrame(
        {
            "yield": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            "promoter": ["p1", "p2", "p3", "p1", "p2", "p3", "p1"],
            "terminator": ["t1", "t1", "t2", "t2", "t1", "t2", "t1"],
        }
    )


def test_auto_ml_uses_the_h2o_api(auto_ml):
    assert auto_ml.h2o is sys.modules["h2o"]
    assert auto_ml.H2OAutoML is FakeH2OAutoML


def test_autoML_trains_on_growing_partitions(auto_ml, ml_df, capsys):
    FakeH2OAutoML.stop_after_training = True
    feature_cols = ["promoter", "terminator"]

    with pytest.raises(StopAfterTraining):
        auto_ml.autoML_on_partitioned_data(
            feature_cols, "yield", ml_df, partitions=3
        )

    # 7 rows in 3 partitions -> nested prefixes of 3, 6 and all 7 rows
    frames = FakeH2OFrame.created
    assert [len(frame) for frame in frames] == [3, 6, 7]
    for frame in frames:
        assert_frame_equal(frame.data, ml_df.iloc[: len(frame)])
        # the features, and only the features, are made categorical
        assert sorted(frame.factor_columns) == sorted(feature_cols)

    # One AutoML run per partition, each trained on its own frame
    assert len(FakeH2OAutoML.instances) == 3
    for automl, frame in zip(FakeH2OAutoML.instances, frames):
        assert automl.init_kwargs == {
            "max_runtime_secs": 5,
            "max_models": None,
            "nfolds": 10,
            "seed": 1,
            "sort_metric": "MAE",
            "keep_cross_validation_predictions": True,
        }
        assert automl.train_kwargs["x"] == feature_cols
        assert automl.train_kwargs["y"] == "yield"
        assert automl.train_kwargs["training_frame"] is frame

    out = capsys.readouterr().out
    assert "len of dataframes that are being trained on : 3" in out
    assert "len of dataframes that are being trained on : 6" in out
    assert "len of dataframes that are being trained on : 7" in out


def test_autoML_passes_training_time_and_nfold(auto_ml, ml_df):
    FakeH2OAutoML.stop_after_training = True

    with pytest.raises(StopAfterTraining):
        auto_ml.autoML_on_partitioned_data(
            ["promoter"],
            "yield",
            ml_df,
            partitions=2,
            training_time=60,
            nfold=3,
        )

    # 7 rows in 2 partitions -> 4 rows and all 7 rows
    assert [len(frame) for frame in FakeH2OFrame.created] == [4, 7]
    assert len(FakeH2OAutoML.instances) == 2
    for automl in FakeH2OAutoML.instances:
        assert automl.init_kwargs["max_runtime_secs"] == 60
        assert automl.init_kwargs["nfolds"] == 3
        assert automl.train_kwargs["x"] == ["promoter"]


def test_autoML_writes_results_csv(auto_ml, ml_df, tmp_path):
    auto_ml.autoML_on_partitioned_data(
        ["promoter", "terminator"],
        "yield",
        ml_df,
        path=str(tmp_path) + "/",
        partitions=3,
    )

    written = list(tmp_path.iterdir())
    assert len(written) == 1
    # e.g. 2026_10_07_15-42_ml_models_running_over_partioned_data.csv
    assert re.fullmatch(
        r"\d{4}_\d{2}_\d{2}_\d{2}-\d{2}_ml_models_running_over_partioned_data\.csv",
        written[0].name,
    )

    results = pd.read_csv(written[0], index_col=0)
    # one row per partition, indexed by the partition size
    assert list(results.index) == [3, 6, 7]
    assert list(results.columns) == ["0", "CV_mean_MAE", "CV_SD_MAE", "Model_name"]
    assert results["0"].tolist() == pytest.approx([10 / 3, 10 / 6, 10 / 7])
    assert results["CV_mean_MAE"].tolist() == pytest.approx([20 / 3, 20 / 6, 20 / 7])
    assert results["CV_SD_MAE"].tolist() == pytest.approx([1 / 3, 1 / 6, 1 / 7])
    assert results["Model_name"].tolist() == [
        "GBM_model_3_rows",
        "GBM_model_6_rows",
        "GBM_model_7_rows",
    ]


def test_autoML_factors_features_with_target_last(auto_ml, capsys):
    # The notebooks' layout: name, part-number features, then the target.
    df = pd.DataFrame(
        {
            "Line_name": [f"yp49_{i}" for i in range(6)],
            "0": [1, 2, 1, 2, 1, 2],
            "1": [3, 3, 4, 4, 5, 5],
            "Amt_norm": [0.1, 0.4, 0.2, 0.8, 0.5, 0.9],
        }
    )
    FakeH2OAutoML.stop_after_training = True

    with pytest.raises(StopAfterTraining):
        auto_ml.autoML_on_partitioned_data(["0", "1"], "Amt_norm", df, partitions=2)

    for frame in FakeH2OFrame.created:
        assert sorted(frame.factor_columns) == ["0", "1"]
