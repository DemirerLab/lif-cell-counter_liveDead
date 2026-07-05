#!/usr/bin/env python3
"""
Count live and dead Arabidopsis thaliana cells from Leica .lif confocal files.

Expected channels:
- C=0: fluorescent cytosol of live cells
- C=1: fluorescent nuclei/puncta of dead cells
- C=2: optional overlay/reference channel, not analyzed by default
"""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from readlif.reader import LifFile
from scipy import ndimage as ndi
from skimage import feature, filters, measure, morphology, segmentation, util


def to_float01(img: np.ndarray, low_pct: float = 1, high_pct: float = 99.5) -> np.ndarray:
    img = util.img_as_float32(img)
    p_low, p_high = np.percentile(img, (low_pct, high_pct))
    if p_high <= p_low:
        return np.clip(img, 0, 1)
    return np.clip((img - p_low) / (p_high - p_low), 0, 1)


def background_subtract(img01: np.ndarray, sigma: float) -> np.ndarray:
    background = filters.gaussian(img01, sigma=sigma, preserve_range=True)
    high_pass = np.clip(img01 - background, 0, None)
    max_val = float(high_pass.max())
    return high_pass / max_val if max_val > 0 else high_pass


def circularity(area: float, perimeter: float) -> float:
    if perimeter <= 0:
        return 0.0
    return float(4.0 * np.pi * area / (perimeter * perimeter))


def get_frame_2d(lif_image, channel: int) -> np.ndarray:
    return np.asarray(lif_image.get_frame(c=channel))


def segment_live_cells(
    img2d: np.ndarray,
    *,
    bg_sigma: float = 25,
    smooth_sigma: float = 1.0,
    sauvola_k: float = 0.22,
    sauvola_window: Optional[int] = None,
    min_area: int = 300,
    max_area: int = 30000,
    hmax_h: float = 2.0,
    fallback_peak_footprint: int = 21,
) -> Tuple[np.ndarray, Dict[str, float]]:
    img01 = to_float01(img2d)
    high_pass = background_subtract(img01, sigma=bg_sigma)
    smoothed = filters.gaussian(high_pass, sigma=smooth_sigma, preserve_range=True)

    if sauvola_window is None:
        sauvola_window = int(max(61, (min(smoothed.shape) // 18) // 2 * 2 + 1))
    if sauvola_window % 2 == 0:
        sauvola_window += 1

    threshold = filters.threshold_sauvola(smoothed, window_size=sauvola_window, k=sauvola_k)
    foreground = smoothed > threshold

    foreground = morphology.remove_small_objects(foreground, min_size=max(1, min_area // 2))
    foreground = morphology.binary_opening(foreground, morphology.disk(2))
    foreground = morphology.binary_closing(foreground, morphology.disk(3))
    foreground = ndi.binary_fill_holes(foreground)

    distance = ndi.distance_transform_edt(foreground)
    maxima = morphology.h_maxima(distance, h=hmax_h)
    markers = measure.label(maxima)

    if markers.max() < 2:
        peaks = feature.peak_local_max(
            distance,
            footprint=np.ones((fallback_peak_footprint, fallback_peak_footprint)),
            labels=foreground,
        )
        marker_img = np.zeros_like(distance, dtype=np.int32)
        for i, (row, col) in enumerate(peaks, start=1):
            marker_img[row, col] = i
        markers = measure.label(marker_img > 0)
        if markers.max() == 0:
            markers = measure.label(foreground)

    labels = segmentation.watershed(-distance, markers, mask=foreground, compactness=0.0)

    keep = np.zeros(labels.max() + 1, dtype=bool)
    for region in measure.regionprops(labels):
        if min_area <= region.area <= max_area:
            keep[region.label] = True

    filtered = labels.copy()
    filtered[~keep[filtered]] = 0
    filtered = measure.label(filtered > 0)

    return filtered, {
        "live_count": int(filtered.max()),
        "live_bg_sigma": bg_sigma,
        "live_smooth_sigma": smooth_sigma,
        "live_sauvola_k": sauvola_k,
        "live_sauvola_window": sauvola_window,
        "live_min_area": min_area,
        "live_max_area": max_area,
        "live_hmax_h": hmax_h,
    }


def segment_dead_nuclei(
    img2d: np.ndarray,
    *,
    bg_sigma: float = 10,
    smooth_sigma: float = 0.8,
    tophat_radius: int = 10,
    min_area: int = 30,
    max_area: int = 4000,
    min_circularity: float = 0.55,
    watershed_peak_footprint: int = 7,
) -> Tuple[np.ndarray, Dict[str, float], np.ndarray]:
    img01 = to_float01(img2d)
    high_pass = background_subtract(img01, sigma=bg_sigma)
    smoothed = filters.gaussian(high_pass, sigma=smooth_sigma, preserve_range=True)

    top_hat = morphology.white_tophat(smoothed, footprint=morphology.disk(tophat_radius))
    top_hat = to_float01(top_hat)

    try:
        threshold = filters.threshold_yen(top_hat)
        threshold_method = "yen"
    except Exception:
        threshold = filters.threshold_otsu(top_hat)
        threshold_method = "otsu"

    foreground = top_hat > threshold
    foreground = morphology.remove_small_objects(foreground, min_size=max(1, min_area - 5))
    foreground = morphology.binary_opening(foreground, morphology.disk(1))
    foreground = ndi.binary_fill_holes(foreground)

    distance = ndi.distance_transform_edt(foreground)
    peaks = feature.peak_local_max(
        distance,
        footprint=np.ones((watershed_peak_footprint, watershed_peak_footprint)),
        labels=foreground,
    )

    marker_img = np.zeros_like(distance, dtype=np.int32)
    for i, (row, col) in enumerate(peaks, start=1):
        marker_img[row, col] = i
    markers = measure.label(marker_img > 0)
    if markers.max() == 0:
        markers = measure.label(foreground)

    labels = segmentation.watershed(-distance, markers, mask=foreground, compactness=0.0)

    keep = np.zeros(labels.max() + 1, dtype=bool)
    for region in measure.regionprops(labels):
        circ = circularity(region.area, region.perimeter)
        if min_area <= region.area <= max_area and circ >= min_circularity:
            keep[region.label] = True

    filtered = labels.copy()
    filtered[~keep[filtered]] = 0
    filtered = measure.label(filtered > 0)

    return filtered, {
        "dead_count": int(filtered.max()),
        "dead_bg_sigma": bg_sigma,
        "dead_smooth_sigma": smooth_sigma,
        "dead_tophat_radius": tophat_radius,
        "dead_threshold": float(threshold),
        "dead_threshold_method": threshold_method,
        "dead_min_area": min_area,
        "dead_max_area": max_area,
        "dead_min_circularity": min_circularity,
    }, top_hat


def annotate_labels(ax, labels: np.ndarray, *, fontsize: int = 6, max_labels: int = 500) -> None:
    for region in measure.regionprops(labels)[:max_labels]:
        y, x = region.centroid
        ax.text(x, y, str(region.label), color="yellow", fontsize=fontsize, ha="center", va="center")


def save_qc_preview(
    output_path: Path,
    *,
    live_raw: np.ndarray,
    dead_raw: np.ndarray,
    live_labels: np.ndarray,
    dead_labels: np.ndarray,
    dead_tophat: np.ndarray,
    title: str,
    annotation_fontsize: int = 6,
) -> None:
    live_norm = to_float01(live_raw)
    dead_norm = to_float01(dead_raw)
    live_boundaries = segmentation.find_boundaries(live_labels, mode="outer")
    dead_boundaries = segmentation.find_boundaries(dead_labels, mode="outer")

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), dpi=130)
    fig.suptitle(title, fontsize=11)

    axes[0, 0].imshow(live_norm, cmap="gray")
    axes[0, 0].set_title("Live C0 normalized")
    axes[0, 0].axis("off")

    axes[0, 1].imshow(live_norm, cmap="gray")
    axes[0, 1].imshow(live_boundaries, alpha=0.9)
    axes[0, 1].set_title(f"Live boundaries n={live_labels.max()}")
    axes[0, 1].axis("off")
    annotate_labels(axes[0, 1], live_labels, fontsize=annotation_fontsize)

    axes[0, 2].imshow(live_labels, cmap="nipy_spectral")
    axes[0, 2].set_title("Live labels")
    axes[0, 2].axis("off")

    axes[1, 0].imshow(dead_norm, cmap="gray")
    axes[1, 0].set_title("Dead C1 normalized")
    axes[1, 0].axis("off")

    axes[1, 1].imshow(dead_tophat, cmap="gray")
    axes[1, 1].set_title("Dead top-hat nuclei enhanced")
    axes[1, 1].axis("off")

    axes[1, 2].imshow(dead_norm, cmap="gray")
    axes[1, 2].imshow(dead_boundaries, alpha=0.9)
    axes[1, 2].set_title(f"Dead boundaries n={dead_labels.max()}")
    axes[1, 2].axis("off")
    annotate_labels(axes[1, 2], dead_labels, fontsize=annotation_fontsize)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def iter_lif_files(input_path: Path) -> Iterable[Path]:
    if input_path.is_file() and input_path.suffix.lower() == ".lif":
        yield input_path
    elif input_path.is_dir():
        yield from sorted(input_path.glob("*.lif"))
    else:
        raise FileNotFoundError(f"Input must be a .lif file or directory containing .lif files: {input_path}")


def process_lif_file(
    lif_path: Path,
    *,
    preview_dir: Path,
    channel_live: int = 0,
    channel_dead: int = 1,
    annotation_fontsize: int = 6,
    live_params: Optional[Dict[str, float]] = None,
    dead_params: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    live_params = live_params or {}
    dead_params = dead_params or {}

    lif = LifFile(str(lif_path))
    images = list(lif.get_iter_image())
    rows = []

    for idx, lif_image in enumerate(images, start=1):
        series_name = lif_image.name
        live_raw = get_frame_2d(lif_image, channel=channel_live)
        dead_raw = get_frame_2d(lif_image, channel=channel_dead)

        live_labels, live_metadata = segment_live_cells(live_raw, **live_params)
        dead_labels, dead_metadata, dead_tophat = segment_dead_nuclei(dead_raw, **dead_params)

        preview_path = preview_dir / f"{lif_path.stem}_Series{idx:02d}_{series_name}_qc.png"
        save_qc_preview(
            preview_path,
            live_raw=live_raw,
            dead_raw=dead_raw,
            live_labels=live_labels,
            dead_labels=dead_labels,
            dead_tophat=dead_tophat,
            title=f"{lif_path.name} | {series_name}",
            annotation_fontsize=annotation_fontsize,
        )

        row = {
            "file": lif_path.name,
            "series_index": idx,
            "series_name": series_name,
            "live_C0_count": int(live_labels.max()),
            "dead_C1_count": int(dead_labels.max()),
        }
        row.update(live_metadata)
        row.update(dead_metadata)
        rows.append(row)

    return pd.DataFrame(rows)


def zip_previews(preview_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for png in sorted(preview_dir.glob("*.png")):
            zf.write(png, arcname=png.name)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Count live and dead Arabidopsis cells from Leica .lif files.")
    parser.add_argument("input", type=Path, help="Path to a .lif file or a directory of .lif files.")
    parser.add_argument("--out", type=Path, default=Path("counts.xlsx"), help="Output Excel file.")
    parser.add_argument("--preview-dir", type=Path, default=Path("previews"), help="Directory for QC preview PNGs.")
    parser.add_argument("--zip-previews", type=Path, default=None, help="Optional output .zip file for preview PNGs.")

    parser.add_argument("--channel-live", type=int, default=0)
    parser.add_argument("--channel-dead", type=int, default=1)
    parser.add_argument("--annotation-fontsize", type=int, default=6)

    parser.add_argument("--live-bg-sigma", type=float, default=25)
    parser.add_argument("--live-smooth-sigma", type=float, default=1.0)
    parser.add_argument("--live-sauvola-k", type=float, default=0.22)
    parser.add_argument("--live-sauvola-window", type=int, default=None)
    parser.add_argument("--live-min-area", type=int, default=300)
    parser.add_argument("--live-max-area", type=int, default=30000)
    parser.add_argument("--live-hmax-h", type=float, default=2.0)

    parser.add_argument("--dead-bg-sigma", type=float, default=10)
    parser.add_argument("--dead-smooth-sigma", type=float, default=0.8)
    parser.add_argument("--dead-tophat-radius", type=int, default=10)
    parser.add_argument("--dead-min-area", type=int, default=30)
    parser.add_argument("--dead-max-area", type=int, default=4000)
    parser.add_argument("--dead-min-circularity", type=float, default=0.55)

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    args.preview_dir.mkdir(parents=True, exist_ok=True)

    live_params = {
        "bg_sigma": args.live_bg_sigma,
        "smooth_sigma": args.live_smooth_sigma,
        "sauvola_k": args.live_sauvola_k,
        "sauvola_window": args.live_sauvola_window,
        "min_area": args.live_min_area,
        "max_area": args.live_max_area,
        "hmax_h": args.live_hmax_h,
    }
    dead_params = {
        "bg_sigma": args.dead_bg_sigma,
        "smooth_sigma": args.dead_smooth_sigma,
        "tophat_radius": args.dead_tophat_radius,
        "min_area": args.dead_min_area,
        "max_area": args.dead_max_area,
        "min_circularity": args.dead_min_circularity,
    }

    all_results = []
    for lif_path in iter_lif_files(args.input):
        all_results.append(
            process_lif_file(
                lif_path,
                preview_dir=args.preview_dir,
                channel_live=args.channel_live,
                channel_dead=args.channel_dead,
                annotation_fontsize=args.annotation_fontsize,
                live_params=live_params,
                dead_params=dead_params,
            )
        )

    result = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_excel(args.out, index=False)

    if args.zip_previews is not None:
        args.zip_previews.parent.mkdir(parents=True, exist_ok=True)
        zip_previews(args.preview_dir, args.zip_previews)

    print(f"Wrote counts: {args.out}")
    print(f"Wrote QC previews to: {args.preview_dir}")
    if args.zip_previews is not None:
        print(f"Wrote preview archive: {args.zip_previews}")


if __name__ == "__main__":
    main()
