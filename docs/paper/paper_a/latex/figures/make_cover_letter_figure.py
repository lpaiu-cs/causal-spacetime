"""Build the cover-letter overview as a deterministic vector figure.

AI-generated concept sketches were used only to compare compositions.  The
submitted PDF/PNG contains no generated raster pixels: every mark below is
declared in code so its scientific meaning can be inspected and reproduced.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

OUT = Path(__file__).resolve().parent
BLUE = "#005A9C"
ORANGE = "#D97706"
INK = "#263238"
MID = "#87929D"
LIGHT = "#D5DEE6"
PALE = "#F7F9FA"


def arrow(ax, start, end, color=INK, lw=0.8, scale=7):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=scale,
        linewidth=lw, color=color, shrinkA=2, shrinkB=2,
    ))


def dot(ax, xy, color=BLUE, size=10):
    ax.scatter(*xy, s=size, color=color, zorder=4, linewidths=0)


def causal_order(ax):
    chains = (
        (((0.45, 2.02), (0.74, 2.66), (0.60, 3.27)), ((0, 1), (1, 2))),
        (((1.15, 2.05), (1.52, 2.60), (1.31, 3.19), (1.83, 3.30)),
         ((0, 1), (1, 2), (1, 3))),
        (((2.04, 2.02), (2.31, 2.63), (2.17, 3.26)), ((0, 1), (1, 2))),
    )
    for points, edges in chains:
        for point in points:
            dot(ax, point)
        for lo, hi in edges:
            arrow(ax, points[lo], points[hi], lw=0.65, scale=6)


def density_icon(ax, center, radius):
    ax.add_patch(Circle(center, radius, fill=False, edgecolor=BLUE, lw=0.9))
    offsets = (
        (-0.19, 0.13), (0.03, 0.19), (0.20, 0.09), (-0.14, -0.05),
        (0.08, -0.02), (0.22, -0.15), (-0.06, -0.20),
    )
    for i, offset in enumerate(offsets):
        color = ORANGE if i in (1, 4) else BLUE
        dot(ax, (center[0] + offset[0], center[1] + offset[1]),
            color=color, size=7)


def clock_icon(ax, center, radius):
    ax.add_patch(Circle(center, radius, fill=False, edgecolor=BLUE, lw=0.9))
    for angle in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        p0 = (center[0] + 0.84 * radius * np.cos(angle),
              center[1] + 0.84 * radius * np.sin(angle))
        p1 = (center[0] + radius * np.cos(angle),
              center[1] + radius * np.sin(angle))
        ax.plot((p0[0], p1[0]), (p0[1], p1[1]), color=MID, lw=0.55)
    ax.plot((center[0], center[0]),
            (center[1], center[1] + 0.22), color=BLUE, lw=1.0)
    ax.plot((center[0], center[0] + 0.18),
            (center[1], center[1] + 0.10), color=BLUE, lw=1.0)
    dot(ax, center, size=8)


def orientation_icon(ax, center, radius):
    ax.add_patch(Circle(center, radius, fill=False, edgecolor=BLUE, lw=0.9))
    for angle, color in ((np.pi / 2, BLUE), (7 * np.pi / 6, ORANGE),
                         (11 * np.pi / 6, MID)):
        end = (center[0] + 0.30 * np.cos(angle),
               center[1] + 0.30 * np.sin(angle))
        arrow(ax, center, end, color=color, lw=1.0, scale=7)
    dot(ax, center, size=8)


def measure_icon(ax, center):
    cx, cy = center
    top = np.array(((cx, cy + 0.34), (cx + 0.32, cy + 0.17),
                    (cx, cy), (cx - 0.32, cy + 0.17)))
    left = np.array(((cx - 0.32, cy + 0.17), (cx, cy),
                     (cx, cy - 0.34), (cx - 0.32, cy - 0.17)))
    right = np.array(((cx + 0.32, cy + 0.17), (cx, cy),
                      (cx, cy - 0.34), (cx + 0.32, cy - 0.17)))
    ax.add_patch(Polygon(top, closed=True, facecolor=PALE,
                         edgecolor=BLUE, lw=0.8))
    ax.add_patch(Polygon(left, closed=True, facecolor="white",
                         edgecolor=BLUE, lw=0.8))
    ax.add_patch(Polygon(right, closed=True, facecolor="#EDF3F8",
                         edgecolor=BLUE, lw=0.8))
    ax.add_patch(Rectangle((cx - 0.13, cy - 0.13), 0.26, 0.26,
                           facecolor="#F6C98D", edgecolor=ORANGE,
                           lw=0.65, alpha=0.8))


def geometry_mesh(ax):
    u_values = np.linspace(-1, 1, 9)
    v_values = np.linspace(-1, 1, 7)
    for u in u_values:
        v = np.linspace(-1, 1, 100)
        x = 9.92 + 1.48 * u * (0.82 + 0.18 * v**2)
        y = 2.66 + 0.68 * v * (0.92 + 0.08 * u**2)
        ax.plot(x, y, color=LIGHT if abs(u) < 0.99 else BLUE,
                lw=0.5 if abs(u) < 0.99 else 0.9)
    for v in v_values:
        u = np.linspace(-1, 1, 100)
        x = 9.92 + 1.48 * u * (0.82 + 0.18 * v**2)
        y = 2.66 + 0.68 * v * (0.92 + 0.08 * u**2)
        ax.plot(x, y, color=LIGHT if abs(v) < 0.99 else BLUE,
                lw=0.5 if abs(v) < 0.99 else 0.9)


def validation_circle(ax, center):
    ax.add_patch(Circle(center, 0.42, facecolor="white",
                        edgecolor=MID, lw=0.65, zorder=2))


def plane_wave(ax, center):
    validation_circle(ax, center)
    cx, cy = center
    x = np.linspace(cx - 0.30, cx + 0.30, 80)
    for shift in (-0.14, -0.04, 0.06, 0.16):
        y = cy + shift + 0.025 * np.sin((x - cx) * 18)
        ax.plot(x, y, color="#73A8D4", lw=0.45)
    for dx in (-0.18, 0.16):
        ax.plot((cx + dx, cx + dx - 0.09, cx + dx),
                (cy - 0.25, cy, cy + 0.25), color=INK, lw=0.45)
        for dy in (-0.11, 0.0, 0.11):
            dot(ax, (cx + dx - 0.02, cy + dy), color=ORANGE, size=3)


def schwarzschild(ax, center):
    validation_circle(ax, center)
    cx, cy = center
    for radius in (0.18, 0.28, 0.37):
        theta = np.linspace(-1.1, 1.1, 80)
        ax.plot(cx + radius * np.sin(theta),
                cy + 0.75 * radius * np.cos(theta), color=LIGHT, lw=0.45)
    t = np.linspace(-1, 1, 80)
    ax.plot(cx - 0.18 + 0.05 * t**2, cy + 0.28 * t,
            color=BLUE, lw=0.9)
    ax.plot(cx + 0.18 - 0.05 * t**2, cy + 0.28 * t,
            color=ORANGE, lw=0.9)
    for x0, color in ((cx - 0.13, BLUE), (cx + 0.13, ORANGE)):
        dot(ax, (x0, cy - 0.28), color=color, size=5)
        dot(ax, (x0, cy + 0.28), color=color, size=5)


def count_volume(ax, center):
    validation_circle(ax, center)
    cx, cy = center
    centers = ((cx - 0.16, cy + 0.14), (cx + 0.16, cy + 0.14),
               (cx - 0.16, cy - 0.15), (cx + 0.16, cy - 0.15))
    for i, (x0, y0) in enumerate(centers):
        width = 0.20 + 0.018 * i
        height = 0.18 + 0.012 * i
        ax.add_patch(Rectangle((x0 - width / 2, y0 - height / 2),
                               width, height, facecolor=PALE,
                               edgecolor=BLUE, lw=0.45))
        offsets = ((-0.05, 0.04), (0.04, 0.03), (-0.02, -0.04),
                   (0.06, -0.03))
        for dx, dy in offsets[: 2 + i // 2]:
            dot(ax, (x0 + dx, y0 + dy), size=2.5)


def build():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(6, 2), dpi=180)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.text(1.40, 3.73, "causal order", ha="center", va="center",
            color=INK, fontsize=9, fontweight="semibold")
    ax.text(5.32, 3.73, "declared structure", ha="center", va="center",
            color=INK, fontsize=9, fontweight="semibold")
    ax.text(9.92, 3.73, "operational geometry", ha="center", va="center",
            color=INK, fontsize=9, fontweight="semibold")

    causal_order(ax)
    arrow(ax, (2.55, 2.66), (3.08, 2.66), color=BLUE, lw=1.1, scale=9)

    icon_centers = ((3.45, 2.66), (4.65, 2.66),
                    (5.85, 2.66), (7.05, 2.66))
    density_icon(ax, icon_centers[0], 0.38)
    clock_icon(ax, icon_centers[1], 0.38)
    orientation_icon(ax, icon_centers[2], 0.38)
    measure_icon(ax, icon_centers[3])
    for center, label in zip(icon_centers,
                             ("density", "clock", "orientation",
                              "local\nmeasure"),
                             strict=True):
        ax.text(center[0], 2.10, label, ha="center", va="center",
                fontsize=8.2, color=MID, linespacing=0.85)

    arrow(ax, (7.52, 2.66), (7.98, 2.66), color=BLUE, lw=1.1, scale=9)
    geometry_mesh(ax)

    ax.text(6.0, 1.42, "preregistered tests at finite density",
            ha="center", va="center", fontsize=8.2, color=INK,
            fontweight="semibold")
    centers = ((4.28, 0.78), (6.00, 0.78), (7.72, 0.78))
    ax.plot((centers[0][0], centers[-1][0]),
            (centers[0][1], centers[-1][1]), color=MID, lw=0.55, zorder=1)
    plane_wave(ax, centers[0])
    schwarzschild(ax, centers[1])
    count_volume(ax, centers[2])
    for center, label in zip(centers,
                             ("plane wave", "Schwarzschild", "count / volume"),
                             strict=True):
        ax.text(center[0], 0.22, label, ha="center", va="center",
                fontsize=8.2, color=MID)

    fig.subplots_adjust(left=0.005, right=0.995, bottom=0.01, top=0.99)
    fig.savefig(
        OUT / "cover_letter_reconstruction.pdf",
        facecolor="white",
        metadata={"Creator": "causal-spacetime-lab",
                  "CreationDate": None, "ModDate": None},
    )
    fig.savefig(OUT / "cover_letter_reconstruction.png", dpi=300,
                facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    build()
