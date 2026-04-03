import pytest
import pandas as pd
import numpy as np

from cellstatemapping.stats import (
    compute_patient_weights,
    compute_cell_proportions,
    weighted_t_test,
    weighted_cluster_comparison
)


@pytest.fixture
def sample_metadata():
    # 20 cells, 4 patients (2 R, 2 NR), 2 clusters
    df = pd.DataFrame({
        "patient": ["P1"]*5 + ["P2"]*5 + ["P3"]*5 + ["P4"]*5,
        "group": ["R"]*10 + ["NR"]*10,
        "cluster": ["C1", "C2", "C1", "C1", "C2"] * 4
    })
    return df


def test_compute_patient_weights():
    # counts = [1000, 2000, 800, 500]
    counts = np.array([1000, 2000, 800, 500])
    weights = compute_patient_weights(counts)
    
    # expected: count * n / sum
    total = sum(counts)
    assert np.isclose(weights[0], 1000 * 4 / total)
    assert sum(weights) == 4.0


def test_compute_cell_proportions(sample_metadata):
    props = compute_cell_proportions(sample_metadata)
    assert "C1" in props.columns
    assert "C2" in props.columns
    assert "group" in props.columns
    
    # P1 has 3 C1 and 2 C2 total 5 cells -> C1: 0.6, C2: 0.4
    assert np.isclose(props.loc["P1", "C1"], 0.6)
    assert np.isclose(props.loc["P1", "C2"], 0.4)
    assert props.loc["P1", "group"] == "R"


def test_weighted_t_test():
    group_a = np.array([0.6, 0.7])
    group_b = np.array([0.2, 0.3])
    weights_a = np.array([1.0, 1.0])
    weights_b = np.array([1.0, 1.0])
    
    t_stat, p_val = weighted_t_test(group_a, group_b, weights_a, weights_b)
    # A > B, so p_val should be small
    assert t_stat > 0
    assert p_val < 0.05


def test_weighted_cluster_comparison(sample_metadata):
    meta = sample_metadata.copy()
    meta.loc[0:4, "cluster"] = ["C1", "C1", "C1", "C1", "C2"] # P1 R: 4 C1, 1 C2
    meta.loc[5:9, "cluster"] = ["C1", "C1", "C1", "C2", "C2"] # P2 R: 3 C1, 2 C2
    meta.loc[10:14, "cluster"] = ["C1", "C2", "C2", "C2", "C2"] # P3 NR: 1 C1, 4 C2
    meta.loc[15:19, "cluster"] = ["C1", "C1", "C2", "C2", "C2"] # P4 NR: 2 C1, 3 C2
    
    res = weighted_cluster_comparison(
        meta,
        n_bootstrap=0 # disable bootstrap for faster test
    )
    
    # Should test C1 and C2
    assert len(res) == 2
    assert "t_statistic" in res.columns
    assert "p_value" in res.columns
    
    # C1 should be significant (R > NR)
    res_c1 = res[res["cluster"] == "C1"].iloc[0]
    assert res_c1["t_statistic"] > 0
    assert res_c1["p_value"] < 0.1
