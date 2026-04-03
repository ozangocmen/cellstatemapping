import pytest
import numpy as np
import pandas as pd

from cellstatemapping.distance import (
    bhattacharyya_distance,
    compute_distance_matrix
)


def test_bhattacharyya_distance():
    p = np.array([0.3, 0.3, 0.4])
    q = np.array([0.2, 0.5, 0.3])
    dist = bhattacharyya_distance(p, q)
    assert dist > 0
    assert dist < 1.0  # Actually between 0 and inf, but this pair is close

    # Same dist -> 0
    dist_same = bhattacharyya_distance(p, p)
    assert np.isclose(dist_same, 0.0, atol=1e-5)

    # Completely disjoint -> high distance (limited by epsilon)
    p_dis = np.array([1.0, 0.0])
    q_dis = np.array([0.0, 1.0])
    dist_dis = bhattacharyya_distance(p_dis, q_dis, epsilon=1e-10)
    assert dist_dis > 5.0  # -ln(1e-10) = 23.0


def test_compute_distance_matrix():
    factor_matrix = pd.DataFrame({
        "Cluster1": [0.5, 0.5],
        "Cluster2": [0.1, 0.9],
        "Cluster3": [0.5, 0.5]
    }, index=["F1", "F2"])
    
    dist_mat = compute_distance_matrix(factor_matrix)
    
    # 3x3 matrix
    assert dist_mat.shape == (3, 3)
    
    # Diagonal should be 0
    assert np.isclose(dist_mat.loc["Cluster1", "Cluster1"], 0.0)
    
    # Cluster1 and Cluster3 are identical -> 0
    assert np.isclose(dist_mat.loc["Cluster1", "Cluster3"], 0.0)
    
    # Cluster1 and Cluster2 are different -> distance > 0
    assert dist_mat.loc["Cluster1", "Cluster2"] > 0
