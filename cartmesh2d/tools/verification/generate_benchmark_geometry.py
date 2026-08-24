#!/usr/bin/env python3
"""Generate deterministic analytic 2-D benchmark boundaries."""

import argparse
import math
from pathlib import Path


def circle(segments: int, radius: float, centre: tuple[float, float], phase_deg: float):
    phase = math.radians(phase_deg)
    return [(centre[0] + radius * math.cos(phase + 2.0 * math.pi * i / segments),
             centre[1] + radius * math.sin(phase + 2.0 * math.pi * i / segments))
            for i in range(segments)]


def naca0012(points_per_surface: int):
    def thickness(x):
        # Closed trailing-edge NACA 00xx polynomial, chord=1 and t/c=0.12.
        return 5.0 * 0.12 * (0.2969 * math.sqrt(x) - 0.1260 * x
                             - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)

    xs = [(1.0 - math.cos(math.pi * i / points_per_surface)) * 0.5
          for i in range(points_per_surface + 1)]
    upper = [(x, thickness(x)) for x in reversed(xs)]
    lower = [(x, -thickness(x)) for x in xs[1:-1]]
    return upper + lower


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("circle", "naca0012"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--segments", type=int, default=128)
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--centre", type=float, nargs=2, default=(0.0, 0.0))
    parser.add_argument("--phase-deg", type=float, default=0.0)
    parser.add_argument("--points-per-surface", type=int, default=128)
    args = parser.parse_args()
    if args.kind == "circle":
        if args.segments < 16:
            raise ValueError("circle benchmark requires at least 16 segments")
        points = circle(args.segments, args.radius, tuple(args.centre), args.phase_deg)
    else:
        if args.points_per_surface < 32:
            raise ValueError("NACA benchmark requires at least 32 points per surface")
        points = naca0012(args.points_per_surface)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("# deterministic analytic benchmark boundary\n" +
                           "".join(f"{x:.17g} {y:.17g}\n" for x, y in points),
                           encoding="utf-8")


if __name__ == "__main__":
    main()
