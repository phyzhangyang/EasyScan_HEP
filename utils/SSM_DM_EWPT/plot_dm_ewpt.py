#!/usr/bin/env python3
"""Plot the scalar-singlet DM and EWPT observables from an EasyScan result."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


RELIC_CONTOURS = (
    (0.012, "#0b3d91", "--", r"$\Omega h^2=0.012$"),
    (0.12, "#d73027", "-", r"$\Omega h^2=0.12$"),
)


def read_scan_result(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    required = {"mS", "lambdaHS", "Omega_h2", "Tc", "true_h_Tc"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise SystemExit(f"Missing required columns in {path}: {', '.join(missing)}")
    data = data.copy()
    data["v_over_t"] = np.where(data["Tc"] > 0, np.abs(data["true_h_Tc"]) / data["Tc"], np.nan)
    return data


def grid_values(data: pd.DataFrame, column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = (
        data.pivot_table(index="lambdaHS", columns="mS", values=column, aggfunc="mean")
        .sort_index()
        .sort_index(axis=1)
    )
    return table.columns.to_numpy(float), table.index.to_numpy(float), table.to_numpy(float)


def grid_edges(values: np.ndarray, log: bool = False) -> np.ndarray:
    if values.size == 1:
        step = values[0] * 0.1 if log else 0.5
        return np.array([values[0] - step, values[0] + step])

    work = np.log(values) if log else values
    edges = np.empty(values.size + 1)
    edges[1:-1] = 0.5 * (work[:-1] + work[1:])
    edges[0] = work[0] - 0.5 * (work[1] - work[0])
    edges[-1] = work[-1] + 0.5 * (work[-1] - work[-2])
    return np.exp(edges) if log else edges


def draw_combined_figure(data: pd.DataFrame, output: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "easyscan-matplotlib"))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output.parent.mkdir(parents=True, exist_ok=True)

    m_grid, lam_grid, omega_grid = grid_values(data, "Omega_h2")
    _, _, tc_grid = grid_values(data, "Tc")
    _, _, vt_grid = grid_values(data, "v_over_t")

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=160)

    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad(color="lightgray")
    vt_for_color = np.ma.masked_where((tc_grid <= 0) | ~np.isfinite(vt_grid), vt_grid)
    mesh = ax.pcolormesh(
        grid_edges(m_grid),
        grid_edges(lam_grid),
        vt_for_color,
        shading="auto",
        cmap=cmap,
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(r"$v_C/T_C$")

    omega_min = float(np.nanmin(omega_grid))
    omega_max = float(np.nanmax(omega_grid))
    for level, color, linestyle, label in RELIC_CONTOURS:
        if not (omega_min <= level <= omega_max):
            continue
        ax.contour(
            m_grid,
            lam_grid,
            omega_grid,
            levels=[level],
            colors="white",
            linestyles=linestyle,
            linewidths=4.0,
            zorder=4,
        )
        contour = ax.contour(
            m_grid,
            lam_grid,
            omega_grid,
            levels=[level],
            colors=color,
            linestyles=linestyle,
            linewidths=2.0,
            zorder=5,
        )
        ax.clabel(
            contour,
            fmt={level: label},
            inline=True,
            fontsize=8,
        )

    ax.scatter(
        data["mS"],
        data["lambdaHS"],
        marker="o",
        s=18,
        facecolors="none",
        edgecolors="black",
        linewidth=0.6,
        alpha=0.55,
        zorder=3,
    )

    ax.set_xlabel(r"$m_S$ [GeV]")
    ax.set_ylabel(r"$\lambda_{HS}$")
    ax.set_xlim(float(np.nanmin(data["mS"])), float(np.nanmax(data["mS"])))
    ax.set_ylim(float(np.nanmin(data["lambdaHS"])), float(np.nanmax(data["lambdaHS"])))
    fig.tight_layout()
    fig.savefig(output)
    if output.suffix.lower() != ".pdf":
        fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "scan_result",
        nargs="?",
        default="SSM_DM_EWPT_grid/ScanResult.txt",
        help="EasyScan ScanResult.txt produced by templates/scan_SSM_DM_EWPT.ini.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="SSM_DM_EWPT_grid/Figures/SSM_DM_EWPT_combined.png",
        help="Output figure path. A PDF copy is also written for non-PDF outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = read_scan_result(Path(args.scan_result))
    draw_combined_figure(data, Path(args.output))


if __name__ == "__main__":
    main()
