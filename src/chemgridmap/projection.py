"""Two-dimensional projection helpers."""

from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler


def project_2d(
    representation: np.ndarray,
    method: str = "pca",
    random_state: int = 42,
    umap_metric: str = "euclidean",
) -> np.ndarray:
    """Project a representation matrix with PCA, t-SNE or UMAP."""
    matrix = np.asarray(representation, dtype=np.float32)
    if matrix.ndim != 2:
        raise ValueError("Representation must be a two-dimensional matrix.")
    if len(matrix) < 2:
        raise ValueError("At least two molecules are required for projection.")
    if not np.isfinite(matrix).all():
        raise ValueError("Representation contains missing or non-finite values.")

    scaled = StandardScaler().fit_transform(matrix).astype(np.float32, copy=False)
    if scaled.shape[1] == 1:
        scaled = np.column_stack([scaled, np.zeros(len(scaled))])

    method = method.lower()
    if method == "pca":
        return PCA(n_components=2, random_state=random_state).fit_transform(scaled)

    if method == "tsne":
        pre_components = min(50, scaled.shape[1], max(1, len(scaled) - 1))
        reduced = PCA(
            n_components=pre_components, random_state=random_state
        ).fit_transform(scaled)
        perplexity = min(30.0, max(1.0, (len(scaled) - 1) / 3.0))
        return TSNE(
            n_components=2,
            perplexity=perplexity,
            init="pca",
            learning_rate="auto",
            random_state=random_state,
        ).fit_transform(reduced)

    if method == "umap":
        try:
            import umap
        except ImportError as exc:
            raise ImportError(
                "UMAP support is optional. Install it with "
                "`pip install -e '.[umap]'` or `pip install umap-learn`."
            ) from exc
        metric = str(umap_metric).strip().lower()
        if metric not in {"euclidean", "jaccard", "cosine"}:
            raise ValueError(
                "umap_metric must be 'euclidean', 'jaccard' or 'cosine'."
            )
        umap_input = matrix if metric == "jaccard" else scaled
        return umap.UMAP(
            n_components=2,
            n_neighbors=min(15, max(2, len(scaled) - 1)),
            min_dist=0.1,
            metric=metric,
            random_state=random_state,
        ).fit_transform(umap_input)

    raise ValueError(
        "Unsupported projection {!r}. Choose pca, tsne or umap.".format(method)
    )
