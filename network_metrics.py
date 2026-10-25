"""Compute system segregation at 20% density and 21 network integrations.

Install dependencies with: pip install numpy nibabel bctpy

Example:
    python network_metrics_20pct.py --timecourses data/timeseries --roi-atlas gm10_84.nii --network-atlas 7_rsn_resampled.nii --output network_metrics.csv

Each text file must be named sub-<ID>-ses<NN>.txt. Its columns must follow
the ascending nonzero labels in the ROI atlas; rows are BOLD time points.
The two NIfTI atlases must share a voxel grid. Use --exclude SUBJECT:SESSION
to omit a scan after quality control (for example, sub-001:ses03).
The source analysis excluded sub-060278:ses03, sub-070336:ses03, and
sub-071636:ses03; pass these through --exclude to reproduce that cohort.
ROIs without overlap with networks 1-7 receive label 8, as in the source code.

Segregation follows segregation.py: negative Pearson correlations are set to
zero, the matrix is proportionally thresholded at 20% density, and only
retained positive edges enter the within- and between-network means.
Integration follows all_pair.py: positive correlations receive a Fisher-z
transform before Louvain community detection (gamma=1), without thresholding.
Louvain and label shuffling are seeded for repeatable results.
"""

import argparse
import csv
import re
from itertools import combinations
from pathlib import Path

import bct
import nibabel as nib
import numpy as np


NETWORKS = {
    1: "VIS", 2: "SMN", 3: "DAN", 4: "VAN",
    5: "LIM", 6: "FPN", 7: "DMN",
}
PAIRS = list(combinations(NETWORKS, 2))
PAIR_COLUMNS = [f"{NETWORKS[a]}_{NETWORKS[b]}_Integration" for a, b in PAIRS]
FILE_PATTERN = re.compile(r"(sub-\d+)-(ses\d+)\.txt")


def roi_network_labels(roi_path, network_path):
    """Assign each ROI to the seven-network label with greatest voxel overlap."""
    roi_img, network_img = nib.load(str(roi_path)), nib.load(str(network_path))
    if roi_img.shape != network_img.shape or not np.allclose(
        roi_img.affine, network_img.affine
    ):
        raise ValueError("ROI and network atlases must have the same voxel grid")

    roi = roi_img.get_fdata()
    network = network_img.get_fdata().astype(int)
    roi_ids = np.unique(roi[roi > 0])
    if not len(roi_ids):
        raise ValueError("The ROI atlas contains no positive labels")
    if not np.allclose(roi_ids, np.rint(roi_ids)):
        raise ValueError("ROI atlas labels must be integers")
    roi_ids = roi_ids.astype(int)

    labels = []
    for roi_id in roi_ids:
        voxels = network[roi == roi_id]
        counts = np.bincount(voxels[(voxels >= 1) & (voxels <= 7)], minlength=8)
        labels.append(int(np.argmax(counts[1:]) + 1) if counts[1:].any() else 8)

    labels = np.asarray(labels)
    if not set(NETWORKS).issubset(set(labels)):
        raise ValueError("At least one of the seven networks has no ROI")
    if np.any(labels == 8):
        print(f"Warning: {np.sum(labels == 8)} ROIs have no seven-network overlap")
    return labels


def load_timecourses(path, n_rois):
    """Load a time-by-ROI matrix, accepting its transpose when unambiguous."""
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] != n_rois and data.shape[0] == n_rois:
        data = data.T
    if data.shape[1] != n_rois or data.shape[0] < 3:
        raise ValueError(f"Expected at least 3 time points and {n_rois} ROI columns")
    if not np.isfinite(data).all() or np.any(np.std(data, axis=0) == 0):
        raise ValueError("Time courses contain nonfinite values or constant ROIs")
    return data


def fisher_connectivity(timecourses):
    """Pearson correlations transformed with Fisher's atanh."""
    correlations = np.corrcoef(timecourses, rowvar=False)
    weights = np.arctanh(np.clip(correlations, -0.9999, 0.9999))
    np.fill_diagonal(weights, 0)
    return weights


def system_segregation(timecourses, labels):
    """Legacy segregation from retained positive edges at 20% density."""
    correlations = np.corrcoef(timecourses, rowvar=False)
    positive = np.maximum(correlations, 0)
    np.fill_diagonal(positive, 0)
    positive = np.nan_to_num(positive)
    thresholded = bct.threshold_proportional(positive, p=0.2)

    i, j = np.triu_indices(len(labels), k=1)
    weights = thresholded[i, j]
    retained = weights > 0
    within = retained & (labels[i] == labels[j])
    between = retained & (labels[i] != labels[j])
    mean_within = weights[within].mean() if within.any() else 0.0
    mean_between = weights[between].mean() if between.any() else 0.0
    return (mean_within - mean_between) / mean_within if mean_within > 0 else 0.0


def pair_integrations(weights, labels, runs=100, permutations=1000, seed=42):
    """Return null-normalized allegiance for all 21 between-network pairs."""
    graph = np.maximum(weights, 0)
    if not np.any(graph):
        raise ValueError("No positive connections are available for Louvain")

    rng = np.random.default_rng(seed)
    allegiance = np.zeros_like(graph)
    for _ in range(runs):
        communities, _ = bct.community_louvain(
            graph, gamma=1.0, seed=int(rng.integers(0, 2**31 - 1))
        )
        allegiance += communities[:, None] == communities[None, :]
    allegiance /= runs

    a, b = np.asarray(PAIRS).T
    counts = np.bincount(labels, minlength=9)
    denominator = counts[a] * counts[b]
    identity = np.eye(9)

    def pair_means(assignment):
        membership = identity[assignment]
        block_sums = membership.T @ allegiance @ membership
        return block_sums[a, b] / denominator

    observed = pair_means(labels)
    null = np.asarray([pair_means(rng.permutation(labels)) for _ in range(permutations)])
    spread = null.std(axis=0)
    return np.divide(
        observed - null.mean(axis=0), spread,
        out=np.full(len(PAIRS), np.nan), where=spread > 0,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timecourses", type=Path, required=True)
    parser.add_argument("--roi-atlas", type=Path, required=True)
    parser.add_argument("--network-atlas", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclude", action="append", default=[], metavar="SUBJECT:SESSION")
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.runs < 1 or args.permutations < 2:
        parser.error("--runs must be positive and --permutations must be at least 2")

    excluded = set()
    for item in args.exclude:
        if ":" not in item:
            parser.error("--exclude must use SUBJECT:SESSION")
        excluded.add(tuple(item.split(":", 1)))

    labels = roi_network_labels(args.roi_atlas, args.network_atlas)
    rows = []
    for path in sorted(args.timecourses.glob("*.txt")):
        match = FILE_PATTERN.fullmatch(path.name)
        if not match or match.group(2) not in {"ses01", "ses02", "ses03"}:
            continue
        subject, session = match.groups()
        if (subject, session) in excluded:
            continue

        timecourses = load_timecourses(path, len(labels))
        weights = fisher_connectivity(timecourses)
        integration = pair_integrations(
            weights, labels, args.runs, args.permutations, args.seed
        )
        row = {
            "subject": subject,
            "session": session,
            "System_Segregation": system_segregation(timecourses, labels),
        }
        row.update(zip(PAIR_COLUMNS, integration))
        rows.append(row)
        print(f"Processed {path.name}")

    if not rows:
        raise ValueError("No matching time-course files were found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["subject", "session", "System_Segregation"] + PAIR_COLUMNS
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} scans to {args.output}")


if __name__ == "__main__":
    main()
