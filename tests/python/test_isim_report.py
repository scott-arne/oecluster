"""Tests for isim() and the approximate iSIM cluster report."""

import math

import numpy as np
import pytest

oefp = pytest.importorskip("oefp")

import oecluster


def _batch(bits, rows):
    return oefp.OEFPBatch.from_fingerprints(
        [oefp.OEFP.from_on_bits(bits, list(on)) for on in rows])


def _random_rows(n, bits, seed):
    rng = np.random.default_rng(seed)
    mask = rng.random((n, bits)) < 0.125
    return [np.flatnonzero(row).tolist() for row in mask]


def _result(labels):
    clusters = {}
    for i, label in enumerate(labels):
        if label >= 0:
            clusters.setdefault(label, []).append(i)
    members = [clusters[k] for k in range(len(clusters))]
    return oecluster.ClusteringResult(labels, members)


def _mixed_labels():
    labels = [i % 7 - 1 for i in range(40)]
    labels[39] = 6
    return labels


def _ratio(inter, union):
    return 1.0 if union == 0 else inter / union


class _Oracle:
    """Brute-force pair sums over an integer bit matrix."""

    def __init__(self, rows, bits):
        self.x = np.zeros((len(rows), bits), dtype=np.int64)
        for i, on in enumerate(rows):
            self.x[i, on] = 1
        self.inter = self.x @ self.x.T
        p = self.x.sum(axis=1)
        self.union = p[:, None] + p[None, :] - self.inter

    def dist(self, i, j):
        return 1.0 - _ratio(int(self.inter[i, j]), int(self.union[i, j]))

    def pair_sums(self, left, right, *, same):
        inter = union = 0
        for a, i in enumerate(left):
            for j in (left[a + 1:] if same else right):
                inter += int(self.inter[i, j])
                union += int(self.union[i, j])
        return inter, union

    def medoid(self, members):
        best, best_num, best_den = None, 0, 0
        for i in members:
            num = sum(int(self.inter[i, j]) for j in members if j != i)
            den = sum(int(self.union[i, j]) for j in members if j != i)
            if den == 0:
                num, den = 1, 1
            if best is None or num * best_den > best_num * den or (
                    num * best_den == best_num * den and i < best):
                best, best_num, best_den = i, num, den
        return best


def _close(actual, expected):
    if math.isnan(expected):
        return math.isnan(actual)
    return actual == pytest.approx(expected, rel=1e-14, abs=0.0)


def test_isim_matches_the_oracle():
    rows = _random_rows(30, 200, 3)
    oracle = _Oracle(rows, 200)
    members = list(range(30))
    inter, union = oracle.pair_sums(members, None, same=True)
    assert _close(oecluster.isim(_batch(200, rows)), _ratio(inter, union))


def test_isim_edge_values():
    assert math.isnan(oecluster.isim(_batch(8, [[1]])))
    assert oecluster.isim(_batch(8, [[], [], []])) == 1.0


def test_isim_validation_order():
    with pytest.raises(TypeError, match="OEFPBatch"):
        oecluster.isim([[1]], metric="dice")
    with pytest.raises(ValueError, match="tanimoto"):
        oecluster.isim(_batch(8, [[1], [2]]), metric="dice")


def test_report_matches_the_oracle_with_noise():
    labels = _mixed_labels()
    rows = _random_rows(40, 70, 2024)
    oracle = _Oracle(rows, 70)
    report = oecluster.isim_report(
        _result(labels), _batch(70, rows), compute_centroid_indices=True,
        compute_per_cluster_records=True, coverage_thresholds=[0.6, 0.75, 0.9])
    k_count = max(labels) + 1
    clusters = [[i for i, lab in enumerate(labels) if lab == k] for k in range(k_count)]
    clustered = [i for i, lab in enumerate(labels) if lab >= 0]

    intra = [oracle.pair_sums(c, None, same=True) for c in clusters]
    expected_intra = 1.0 - _ratio(sum(a for a, _ in intra), sum(u for _, u in intra))
    assert _close(report.isim_intra_distance, expected_intra)

    cross_i = cross_u = 0
    for k in range(k_count):
        for m in range(k + 1, k_count):
            a, u = oracle.pair_sums(clusters[k], clusters[m], same=False)
            cross_i += a
            cross_u += u
    assert _close(report.isim_inter_distance, 1.0 - _ratio(cross_i, cross_u))

    medoids = [oracle.medoid(c) for c in clusters]
    assert [r.medoid for r in report.records] == medoids
    for k, record in enumerate(report.records):
        outside = [j for j in clustered if labels[j] != k]
        a, u = oracle.pair_sums(clusters[k], outside, same=False)
        assert _close(record.isim_separation, 1.0 - _ratio(a, u))
        dists = [oracle.dist(medoids[k], j) for j in clusters[k]]
        assert _close(record.radius, max(dists))
        expected_mean = (0.0 if len(clusters[k]) < 2
                         else sum(d for j, d in zip(clusters[k], dists) if j != medoids[k])
                         / (len(clusters[k]) - 1))
        assert _close(record.mean_medoid_distance, expected_mean)

        best, nearest = -math.inf, -1
        for m in range(k_count):
            if m == k:
                continue
            a, u = oracle.pair_sums(clusters[k], clusters[m], same=False)
            if _ratio(a, u) > best:
                best, nearest = _ratio(a, u), m
        assert record.nearest_cluster == nearest
        assert _close(record.nearest_cluster_similarity, best)
        expected_record_intra = (math.nan if len(clusters[k]) < 2
                                 else 1.0 - _ratio(*intra[k]))
        assert _close(record.isim_intra_distance, expected_record_intra)

    radii, means, scatter, square = [], [], [], []
    for k in range(k_count):
        dists = [oracle.dist(medoids[k], j) for j in clusters[k]]
        radii.append(max(dists))
        means.append(0.0 if len(clusters[k]) < 2
                     else sum(d for j, d in zip(clusters[k], dists) if j != medoids[k])
                     / (len(clusters[k]) - 1))
        scatter.append(sum(dists) / len(clusters[k]))
        square.append(sum(d * d for d in dists))
    assert _close(report.median_radius, float(np.median(radii)))
    assert _close(report.median_medoid_member_distance, float(np.median(means)))

    global_medoid = oracle.medoid(clustered)
    between = sum(len(clusters[k]) * oracle.dist(medoids[k], global_medoid) ** 2
                  for k in range(k_count))
    within = sum(square)
    expected_ch = ((between / (k_count - 1)) / (within / (len(clustered) - k_count))
                   if within != 0.0 else math.nan)
    assert _close(report.calinski_harabasz_medoid, expected_ch)

    sil_total = 0.0
    for k in range(k_count):
        sil_sum = 0.0
        for i in clusters[k]:
            s = 0.0
            if len(clusters[k]) >= 2:
                a = 1.0 - _ratio(*oracle.pair_sums([i], [j for j in clusters[k] if j != i],
                                                   same=False))
                b = min(1.0 - _ratio(*oracle.pair_sums([i], clusters[m], same=False))
                        for m in range(k_count) if m != k)
                top = max(a, b)
                s = 0.0 if top == 0.0 else (b - a) / top
            sil_sum += s
        sil_total += sil_sum
        assert _close(report.records[k].isim_silhouette, sil_sum / len(clusters[k]))
    assert _close(report.isim_silhouette, sil_total / len(clustered))

    db_total, min_sep, max_spread = 0.0, math.inf, 0.0
    for a in range(k_count):
        worst = 0.0
        for b in range(k_count):
            if a == b:
                continue
            sep = oracle.dist(medoids[a], medoids[b])
            min_sep = min(min_sep, sep)
            worst = max(worst, math.inf if sep == 0.0 else (scatter[a] + scatter[b]) / sep)
        db_total += worst
        max_spread = max(max_spread, 2.0 * scatter[a])
    assert _close(report.davies_bouldin_medoid, db_total / k_count)
    assert _close(report.dunn_medoid_separation_medoid_spread,
                  math.nan if max_spread == 0.0 else min_sep / max_spread)

    num_noise = labels.count(-1)
    for t, threshold in enumerate(report.coverage_thresholds):
        near = [min(oracle.dist(i, m) for m in medoids) for i in range(len(labels))]
        covered = [i for i in range(len(labels)) if near[i] <= threshold]
        assert report.coverage_at[t] == len(covered) / len(labels)
        assert report.noise_coverage_at[t] == (
            sum(1 for i in covered if labels[i] < 0) / num_noise)


def test_report_is_read_only_and_typed():
    report = oecluster.isim_report(
        _result([0, 0, 1, 1]), _batch(16, [[1], [1, 2], [5], [5, 6]]),
        compute_per_cluster_records=True)
    assert isinstance(report, oecluster.ISimReport)
    assert isinstance(report.records[0], oecluster.ISimClusterRecord)
    assert report.requested == oecluster.ISimReportRequested(
        centroid_indices=False, per_cluster_records=True)
    assert report.method == ""
    assert report.coverage_thresholds == (0.25, 0.35, 0.45)
    assert report.coverage_at == ()
    assert math.isnan(report.isim_silhouette)
    with pytest.raises(AttributeError):
        report.num_clusters = 3
    assert repr(report).startswith("ISimReport(method=''")
    assert oecluster.ISimClusterRecord().nearest_cluster == -1
    assert oecluster.ISimReportRequested() == oecluster.ISimReportRequested(
        centroid_indices=False, per_cluster_records=False)
    assert math.isnan(oecluster.ISimClusterRecord().isim_separation)


def test_no_clusters_records_the_request():
    report = oecluster.isim_report(
        _result([-1, -1]), _batch(8, [[1], [2]]),
        compute_centroid_indices=True, compute_per_cluster_records=True)
    assert report.num_clusters == 0
    assert report.records == ()
    assert report.requested.centroid_indices is True
    assert math.isnan(report.isim_intra_distance)


def test_preset_seeds_thresholds():
    report = oecluster.isim_report(
        _result([0, 0]), _batch(8, [[1], [2]]), preset="tight")
    assert report.coverage_thresholds == (0.20, 0.30, 0.40)


def test_numpy_bool_flags_are_accepted():
    report = oecluster.isim_report(
        _result([0, 0, 1]), _batch(8, [[1], [2], [3]]),
        compute_centroid_indices=np.bool_(True))
    assert report.requested.centroid_indices is True


@pytest.mark.parametrize("kwargs, exc, message", [
    ({"metric": "dice"}, ValueError, "tanimoto"),
    ({"preset": "loose"}, ValueError, "preset"),
    ({"coverage_thresholds": [float("nan")]}, ValueError,
     "coverage thresholds must not be NaN"),
    ({"coverage_thresholds": [-0.1]}, ValueError,
     "coverage thresholds must be non-negative"),
    ({"treat_noise_as_singletons": "no"}, TypeError,
     "treat_noise_as_singletons must be True or False, not str"),
    ({"compute_centroid_indices": 1}, TypeError,
     "compute_centroid_indices must be True or False, not int"),
    ({"compute_per_cluster_records": None}, TypeError,
     "compute_per_cluster_records must be True or False, not NoneType"),
    ({"num_threads": -1}, ValueError, "num_threads must be non-negative"),
])
def test_report_validation(kwargs, exc, message):
    with pytest.raises(exc, match=message):
        oecluster.isim_report(_result([0, 0]), _batch(8, [[1], [2]]), **kwargs)


def test_report_validation_order():
    batch = _batch(8, [[1], [2]])
    with pytest.raises(TypeError, match="ClusteringResult"):
        oecluster.isim_report([0, 0], [[1]], metric="dice")
    with pytest.raises(TypeError, match="OEFPBatch"):
        oecluster.isim_report(_result([0, 0]), [[1]], metric="dice")
    with pytest.raises(ValueError, match="tanimoto"):
        oecluster.isim_report(_result([0, 0]), batch, metric="dice", preset="loose")
    # A bad flag outranks a negative num_threads and a cardinality mismatch.
    with pytest.raises(TypeError, match="compute_centroid_indices"):
        oecluster.isim_report(_result([0, 0, 0]), batch,
                              compute_centroid_indices="yes", num_threads=-1)
    with pytest.raises(ValueError, match="num_threads"):
        oecluster.isim_report(_result([0, 0, 0]), batch, num_threads=-1)
    with pytest.raises(ValueError, match="3 samples and the batch 2"):
        oecluster.isim_report(_result([0, 0, 0]), batch)
