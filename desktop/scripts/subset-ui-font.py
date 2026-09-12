#!/usr/bin/env python3
"""Build the bundled CartMesh UI WOFF from the official Noto Sans CJK SC OTF."""

from __future__ import annotations

import argparse
import hashlib
import pathlib

from fontTools import subset
from fontTools.ttLib import TTFont


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
DESKTOP_DIR = SCRIPT_DIR.parent
PROJECT_DIR = DESKTOP_DIR.parent
DEFAULT_SOURCE_DIR = PROJECT_DIR / "outputs" / "font-source"
DEFAULT_OUTPUT_DIR = DESKTOP_DIR / "src" / "renderer" / "assets" / "fonts"
SOURCE_URL = (
    "https://github.com/notofonts/noto-cjk/blob/main/Sans/OTF/"
    "SimplifiedChinese/NotoSansCJKsc-Regular.otf"
)
LICENSE_URL = "https://github.com/notofonts/noto-cjk/blob/main/Sans/LICENSE"


def gb2312_characters() -> set[str]:
    characters: set[str] = set()
    for lead in range(0xA1, 0xF8):
        for trail in range(0xA1, 0xFF):
            try:
                characters.update(bytes((lead, trail)).decode("gb2312"))
            except UnicodeDecodeError:
                pass
    return characters


def renderer_characters() -> set[str]:
    characters: set[str] = set()
    for source in sorted((DESKTOP_DIR / "src" / "renderer").rglob("*")):
        if source.suffix in {".css", ".html", ".js"}:
            characters.update(
                character
                for character in source.read_text(encoding="utf-8")
                if ord(character) >= 0x20
            )
    return characters


def rename_font(font: TTFont) -> None:
    names = font["name"]
    values = {
        1: "CartMesh UI",
        2: "Regular",
        3: "2.004;CartMesh2D;CartMeshUI-Regular",
        4: "CartMesh UI",
        6: "CartMeshUI-Regular",
        16: "CartMesh UI",
        17: "Regular",
    }
    for name_id, value in values.items():
        names.removeNames(nameID=name_id)
        names.setName(value, name_id, 3, 1, 0x409)
        names.setName(value, name_id, 1, 0, 0)


def build(source_dir: pathlib.Path, output_dir: pathlib.Path) -> None:
    source_font = source_dir / "NotoSansCJKsc-Regular.otf"
    source_license = source_dir / "OFL.txt"
    if not source_font.is_file() or not source_license.is_file():
        raise SystemExit(
            "Missing official source font or OFL license under "
            f"{source_dir}. See {SOURCE_URL} and {LICENSE_URL}."
        )

    characters = gb2312_characters() | renderer_characters()
    characters.update(chr(codepoint) for codepoint in range(0x20, 0x7F))

    options = subset.Options()
    options.flavor = "woff"
    options.layout_features = ["*"]
    options.name_IDs = [0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 17]
    options.name_languages = ["*"]
    options.recommended_glyphs = True
    font = subset.load_font(str(source_font), options)
    subsetter = subset.Subsetter(options=options)
    subsetter.populate(unicodes=sorted(ord(character) for character in characters))
    subsetter.subset(font)
    rename_font(font)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_font = output_dir / "cartmesh-ui-regular.woff"
    subset.save_font(font, str(output_font), options)

    source_hash = hashlib.sha256(source_font.read_bytes()).hexdigest()
    notice = (
        "CartMesh UI is a renamed, character-subsetted derivative of "
        "Noto Sans CJK SC Regular version 2.004.\n"
        f"Source: {SOURCE_URL}\n"
        f"Source SHA-256: {source_hash}\n"
        "The derivative contains the GB2312 repertoire, ASCII, and characters "
        "used by the CartMesh2D renderer.\n\n"
    )
    (output_dir / "LICENSE.txt").write_text(
        notice + source_license.read_text(encoding="utf-8"), encoding="utf-8"
    )

    cmap = set().union(*(table.cmap.keys() for table in font["cmap"].tables))
    missing = sorted(ord(character) for character in characters if ord(character) not in cmap)
    if missing:
        preview = ", ".join(f"U+{codepoint:04X}" for codepoint in missing[:12])
        raise SystemExit(f"Source font lacks {len(missing)} requested characters: {preview}")
    print(
        f"wrote {output_font} ({output_font.stat().st_size} bytes, "
        f"{len(cmap)} encoded characters)"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=pathlib.Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output-dir", type=pathlib.Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    build(args.source_dir.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
