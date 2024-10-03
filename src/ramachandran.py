from importlib import resources
from itertools import cycle
from typing import List

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from ramachandraw.parser import get_phi_psi


def ramachandran_canvas(
        ax: plt.Axes,
        line_width: float = 1,
        line_color: str = "black",
        bg_alpha: float = 0.9,
        bg_cmap: str | None | mcolors.Colormap = mcolors.LinearSegmentedColormap.from_list(
            "ramachandran_cmap",
            [
                (0.17, 0.17, 0.17),
                (0.46, 0.46, 0.46),
                (0.85, 0.85, 0.85),
                (1.00, 1.00, 1.00),
            ]),
):
    density_file = resources.files(package="ramachandraw") / "kde.dat"
    with density_file.open("r") as kde_data:
        z = np.fromfile(kde_data)
    z = np.reshape(z, (100, 100))

    # Normalize
    if bg_cmap is not None:
        data = np.log10(np.rot90(z))
        img = ax.imshow(
            data, cmap=plt.get_cmap(bg_cmap), extent=(-180, 180, -180, 180), alpha=bg_alpha
        )
        cbar = plt.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.text(1.15, 1.02, 'favored', transform=cbar.ax.transAxes, ha='right', va='bottom')
        cbar.set_ticks([])
        cbar.ax.tick_params(labelsize=0)

    # Add contour lines
    data = np.rot90(np.fliplr(z))
    ax.contour(
        data,
        colors=line_color,
        linewidths=line_width,
        levels=[10 ** i for i in range(-7, 0)],
        antialiased=True,
        extent=[-180, 180, -180, 180],
        alpha=0.5,
    )

    ticks = list(range(-180, 181, 45))
    ax.set(
        aspect="equal",
        xlabel="\u03C6",
        ylabel="\u03C8",
        xlim=(-180, 180),
        ylim=(-180, 180),
        xticks=ticks,
        yticks=ticks,
    )
    plt.axhline(y=0, color="k", lw=0.6)
    plt.axvline(x=0, color="k", lw=0.6)
    plt.grid(visible=None, which="major", axis="both", color="black", alpha=0.3)
    plt.xticks(rotation=90)

    for spine in ax.spines.values():
        spine.set_edgecolor("k")
        spine.set_linewidth(0.6)


def ramachandran_pdbs(
        pdb: str | List[str],
        ax: plt.Axes,
        dot_alpha: float = 0.7,
        dot_color: str | tuple = "black",
        dot_size: float = 2,
        dot_marker: str = "x",
        legend: bool = False,
):
    angles = get_phi_psi(pdb_filepath=pdb)

    def draw(data: dict[str, list], color: str) -> None:
        """
        Plot the aminoacid residues given their Φ-Ψ torsion angles.
        :param dict[str, list] data: angles data (e.g. {"A:ARG156": [-29.8, -32.4]})
        :param str color: marker color, defaults to "k" (black)
        """
        x = [torsion[0] for torsion in data.values()]
        y = [torsion[1] for torsion in data.values()]
        ax.scatter(x=x, y=y, marker=dot_marker, s=dot_size, color=color, alpha=dot_alpha)

    # Single PDB
    if isinstance(angles, dict):
        draw(data=angles, color=dot_color)

    # Multiple PDBs
    elif isinstance(angles, list):
        ax.set_title(f"Batch ({len(pdb)} files)")
        colors = cycle(mcolors.TABLEAU_COLORS.values())
        custom_legend = []
        for pdb, data in zip(pdb, angles):
            color = next(colors)
            draw(data=data, color=dot_color)
            point = Line2D(
                [0],
                [0],
                label=pdb,
                marker="o",
                markerfacecolor=color,
                markeredgewidth=0,
                markersize=5,
                linestyle="",
            )
            custom_legend.append(point)
        handles, _ = plt.gca().get_legend_handles_labels()
        handles.extend(custom_legend)
        if legend:
            ax.legend(handles=handles, loc=1)
