"""Reference transcription of the published MODI and RMODI definitions.

Transcribed from ``src/rogi/modi.py`` in coleygroup/rogi (MIT licence), which
implements Golbraikh et al. (2014) for MODI and Ruiz and Gomez-Nieto (2018)
for RMODI. rogi is not a dependency of this project, and nothing imports this
module: it exists so the parity literals pinned in
``tests/python/test_sar_coherence.py`` can be re-derived by hand if the C++
implementation ever has to be argued from first principles again.

Run it directly to print the values the parity tests pin::

    /Users/johnss51/Applications/micromamba/envs/main/bin/python \\
        tests/assets/sar_coherence_reference.py
"""

import numpy as np


def _square(condensed, n):
    """Expand a condensed upper-triangle vector into an n x n matrix."""
    square = np.zeros((n, n), dtype=float)
    rows, cols = np.triu_indices(n, k=1)
    square[rows, cols] = condensed
    square[cols, rows] = condensed
    return square


def _nearest_neighbours(square):
    """Index of each row's nearest other sample, ties to the lowest index.

    The diagonal is masked with infinity rather than skipped, so a zero
    self-distance cannot win the argmin. numpy's argmin returns the first
    minimum, which is the tie rule the C++ sweep applies as well.
    """
    masked = square.copy()
    np.fill_diagonal(masked, np.inf)
    return np.argmin(masked, axis=1)


def modi(condensed, classes):
    """MODI: the mean over classes of the same-class nearest-neighbour rate.

    :param condensed: Condensed upper-triangle distances.
    :param classes: One class label per sample.
    :returns: ``(modi, {class: fraction})``, the classes in first-appearance
        order.
    """
    labels = [str(value) for value in classes]
    neighbours = _nearest_neighbours(_square(condensed, len(labels)))
    fractions = {}
    for name in dict.fromkeys(labels):
        members = [i for i, label in enumerate(labels) if label == name]
        same = sum(1 for i in members if labels[neighbours[i]] == name)
        fractions[name] = same / len(members)
    return sum(fractions.values()) / len(fractions), fractions


def rmodi(condensed, activity, delta=0.625):
    """RMODI: the continuous-activity analogue, over an activity band.

    A sample counts when its nearest neighbour inside the band -- within
    ``delta`` population standard deviations of its own activity -- is strictly
    closer than its nearest neighbour outside it.

    :param condensed: Condensed upper-triangle distances.
    :param activity: One measurement per sample.
    :param delta: Band half-width in population standard deviations.
    :returns: The fraction of samples that count, as a float.
    """
    values = np.asarray(activity, dtype=float)
    square = _square(condensed, len(values))
    # numpy's std is the population form, which is the one the band is
    # defined in.
    band = delta * float(values.std())
    inside = np.abs(values[:, None] - values[None, :]) <= band
    np.fill_diagonal(inside, False)
    outside = ~inside
    np.fill_diagonal(outside, False)
    same_min = np.where(inside, square, np.inf).min(axis=1)
    diff_min = np.where(outside, square, np.inf).min(axis=1)
    return float(np.mean(same_min < diff_min))


if __name__ == "__main__":
    # Six points on a line at 0.0, 0.60, 0.40, 0.45, 0.80 and 0.85, so every
    # distance below is |xi - xj| and the whole fixture is a metric.
    CONDENSED = [0.60, 0.40, 0.45, 0.80, 0.85, 0.20, 0.15, 0.20, 0.25,
                 0.05, 0.40, 0.45, 0.35, 0.40, 0.05]
    ACTIVITY = [0.0, 0.2, 1.0, 1.2, 2.0, 2.2]
    CLASSES = ["A", "B", "A", "A", "B", "B"]
    print(modi(CONDENSED, CLASSES))
    # The parity tests pin the band's standard deviation as well, so print it
    # rather than leave it buried inside rmodi. numpy's summation order puts
    # this one ulp above the sqrt(4.06 / 6) that test pins, far inside its
    # abs=1e-12 tolerance.
    print(float(np.std(ACTIVITY)))
    print(rmodi(CONDENSED, ACTIVITY))
