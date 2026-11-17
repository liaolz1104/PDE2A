#!/usr/bin/env python3
"""Compose the manuscript-development Fig. 2a-j from frozen standalone panels.

This script only lays out existing PNG panels. It does not recalculate or alter
any statistic. The generated source_data.tsv records the provenance of every
component placed in the composite figure.
"""

from __future__ import annotations

import csv
import base64
import html
import os
import sys
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont


FIGURE_ROOT = Path("results/07_figures/manuscript_redraw_v3_single_panel")
OUTPUT_DIR = Path("results/07_figures/manuscript_composites/Figure_02_Transcriptome_Overview")
OUTPUT_NAME = "Figure_02_Transcriptome_Overview.png"
OUTPUT_SVG_NAME = "Figure_02_Transcriptome_Overview.svg"


def font_path(bold: bool = False) -> str | None:
    candidates = [
        Path(r"C:/Windows/Fonts/arialbd.ttf" if bold else r"C:/Windows/Fonts/arial.ttf"),
        Path(r"C:/Windows/Fonts/Arial_Bold.ttf" if bold else r"C:/Windows/Fonts/Arial.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    path = font_path(bold=bold)
    if path is not None:
        return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def crop_to_content(image: Image.Image, padding: int = 24) -> Image.Image:
    """Trim only the uniform white border while preserving all plot content."""
    rgb = image.convert("RGB")
    background = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, background)
    bbox = diff.getbbox()
    if bbox is None:
        return rgb
    left = max(0, bbox[0] - padding)
    top = max(0, bbox[1] - padding)
    right = min(rgb.width, bbox[2] + padding)
    bottom = min(rgb.height, bbox[3] + padding)
    return rgb.crop((left, top, right, bottom))


def fit_contain(image: Image.Image, max_width: int, max_height: int) -> Image.Image:
    scale = min(max_width / image.width, max_height / image.height)
    width = max(1, round(image.width * scale))
    height = max(1, round(image.height * scale))
    return image.resize((width, height), Image.Resampling.LANCZOS)


def component_png(panel: str, component: str) -> Path:
    names = {
        "a1": "Figure_01_PairedSamplingDesign",
        "a2": "Figure_02_PairedPCA_CoreEdge",
        "b": "Figure_04_MA",
        "c": "Figure_05_Volcano_PDE2A_Polished",
        "d": "Figure_06_TopDEG_to_PDE2A_CoreEdgeGroupedHeatmap",
        "e": "Figure_07_GSEA_Hallmark",
        "f1": "Figure_08_GSEA_GO_BP",
        "f2": "Figure_09_GSEA_Reactome",
        "g": "Figure_12_GeneModule_LeadingEdgeHeatmap",
        "h": "Figure_03_PDE2A_PairedExpression",
        "i": "Figure_18_LOPO_Stability_Top20",
        "j": "Figure_21_Pareto_Global",
    }
    folder = Path.cwd() / FIGURE_ROOT / names[component]
    matches = list(folder.glob("*.png"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one PNG for {panel}/{component}; found {len(matches)} in {folder}")
    return matches[0]


def source_paths(project_root: Path, panel: str, component: str) -> tuple[str, str, str, Path]:
    png = component_png(panel, component)
    relative_folder = png.parent.relative_to(project_root)
    source_data = relative_folder / "source_data.tsv"
    plotting_script = relative_folder / "plotting_script.R"
    return (
        png.relative_to(project_root).as_posix(),
        source_data.as_posix(),
        plotting_script.as_posix(),
        png,
    )


def component_layouts(cell_box: tuple[int, int, int, int], panel: str,
                       components: list[tuple[str, str]], orientation: str,
                       cell_gap: int, label_height: int,
                       inner_padding: int) -> list[dict[str, object]]:
    """Return exact panel image placements shared by PNG and editable SVG."""
    x0, y0, x1, y1 = cell_box
    content_x0 = x0 + inner_padding
    content_y0 = y0 + label_height
    content_x1 = x1 - inner_padding
    content_y1 = y1 - inner_padding
    content_width = max(1, content_x1 - content_x0)
    content_height = max(1, content_y1 - content_y0)
    images = [(component, component_png(panel, component),
               crop_to_content(Image.open(component_png(panel, component))))
              for _, component in components]
    placements: list[dict[str, object]] = []
    if orientation == "horizontal":
        item_width = max(1, (content_width - cell_gap * (len(images) - 1)) // len(images))
        for index, (component, png_path, image) in enumerate(images):
            fitted = fit_contain(image, item_width, content_height)
            x = content_x0 + index * (item_width + cell_gap) + (item_width - fitted.width) // 2
            y = content_y0 + (content_height - fitted.height) // 2
            placements.append({"panel": panel, "component": component, "png": png_path,
                               "image": fitted, "x": x, "y": y,
                               "width": fitted.width, "height": fitted.height})
    else:
        item_height = max(1, (content_height - cell_gap * (len(images) - 1)) // len(images))
        for index, (component, png_path, image) in enumerate(images):
            fitted = fit_contain(image, content_width, item_height)
            x = content_x0 + (content_width - fitted.width) // 2
            y = content_y0 + index * (item_height + cell_gap) + (item_height - fitted.height) // 2
            placements.append({"panel": panel, "component": component, "png": png_path,
                               "image": fitted, "x": x, "y": y,
                               "width": fitted.width, "height": fitted.height})
    return placements


def natural_height(components: list[tuple[str, str]], width: int, gap: int, vertical: bool = False) -> int:
    images = [crop_to_content(Image.open(component_png(panel, component))) for panel, component in components]
    if vertical:
        available = width
        return sum(round(img.height * available / img.width) for img in images) + gap * (len(images) - 1)
    available = max(1, (width - gap * (len(images) - 1)) // len(images))
    return max(round(img.height * available / img.width) for img in images)


def place_components(canvas: Image.Image, cell_box: tuple[int, int, int, int], panel: str,
                     components: list[tuple[str, str]], orientation: str, panel_label: str,
                     cell_gap: int, label_height: int, inner_padding: int) -> None:
    draw = ImageDraw.Draw(canvas)
    x0, y0, x1, y1 = cell_box
    draw.text((x0 + inner_padding, y0), panel_label, fill="#111111",
              font=load_font(64, bold=True))
    for placement in component_layouts(cell_box, panel, components, orientation,
                                       cell_gap, label_height, inner_padding):
        canvas.paste(placement["image"], (placement["x"], placement["y"]))


def png_data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def write_editable_svg(output_svg: Path, canvas_width: int, canvas_height: int,
                       labels: list[tuple[str, int, int]],
                       placements: list[dict[str, object]]) -> None:
    """Write a portable SVG with one top-level image object per component."""
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{canvas_width}px" height="{canvas_height}px" '
        f'viewBox="0 0 {canvas_width} {canvas_height}">',
        '<title>Figure 2 transcriptome overview; independently movable panel objects</title>',
        f'<rect id="figure_background" x="0" y="0" width="{canvas_width}" height="{canvas_height}" fill="#ffffff"/>',
    ]
    for label, x, y in labels:
        safe_label = html.escape(label, quote=True)
        lines.append(
            f'<text id="label_{safe_label}" x="{x}" y="{y}" fill="#111111" '
            'font-family="Arial, Helvetica, sans-serif" font-size="64px" font-weight="700" '
            'dominant-baseline="hanging">'
            f'{safe_label}</text>'
        )
    for placement in placements:
        component = str(placement["component"])
        panel = str(placement["panel"])
        png_path = Path(placement["png"])
        href = png_data_uri(png_path)
        object_id = f"panel_{panel}_{component}"
        lines.append(
            f'<image id="{object_id}" x="{placement["x"]}" y="{placement["y"]}" '
            f'width="{placement["width"]}" height="{placement["height"]}" '
            'preserveAspectRatio="none" '
            f'href="{href}" xlink:href="{href}"/>'
        )
    lines.append('</svg>')
    output_svg.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: 39_compose_figure2_transcriptome_overview.py <project_root>", file=sys.stderr)
        return 2
    project_root = Path(sys.argv[1]).resolve()
    os.chdir(project_root)
    figure_root = project_root / FIGURE_ROOT
    if not figure_root.is_dir():
        raise FileNotFoundError(f"Missing standalone figure root: {figure_root}")

    # Three-column, five-row layout. Reading order is a-j. The information-dense
    # heatmap, convergence map, GSEA pair, and Pareto plot receive double width.
    canvas_width = 6000
    margin = 120
    column_gap = 90
    row_gap = 100
    inner_padding = 42
    label_height = 92
    column_width = (canvas_width - 2 * margin - 2 * column_gap) // 3
    wide_width = column_width * 2 + column_gap

    cell_specs = [
        ("a", [("a", "a1"), ("a", "a2")], "horizontal", 0, 0, 2),
        ("b", [("b", "b")], "horizontal", 2, 0, 1),
        ("c", [("c", "c")], "horizontal", 0, 1, 1),
        ("d", [("d", "d")], "horizontal", 1, 1, 2),
        ("e", [("e", "e")], "horizontal", 0, 2, 1),
        ("f", [("f", "f1"), ("f", "f2")], "horizontal", 1, 2, 2),
        ("g", [("g", "g")], "horizontal", 0, 3, 2),
        ("h", [("h", "h")], "horizontal", 2, 3, 1),
        ("i", [("i", "i")], "horizontal", 0, 4, 1),
        ("j", [("j", "j")], "horizontal", 1, 4, 2),
    ]

    row_heights = [0, 0, 0, 0, 0]
    for panel, components, orientation, col, row, colspan in cell_specs:
        width = column_width * colspan + column_gap * (colspan - 1)
        content_width = width - 2 * inner_padding
        if orientation == "horizontal" and len(components) > 1:
            height = natural_height(components, content_width, 60, vertical=False)
        elif orientation == "vertical":
            height = natural_height(components, content_width, 60, vertical=True)
        else:
            image = crop_to_content(Image.open(component_png(panel, components[0][1])))
            height = round(image.height * content_width / image.width)
        row_heights[row] = max(row_heights[row], height + label_height + inner_padding)

    canvas_height = 2 * margin + sum(row_heights) + row_gap * (len(row_heights) - 1)
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    y_positions = []
    y = margin
    for height in row_heights:
        y_positions.append(y)
        y += height + row_gap

    for panel, components, orientation, col, row, colspan in cell_specs:
        width = column_width * colspan + column_gap * (colspan - 1)
        x = margin + col * (column_width + column_gap)
        place_components(canvas, (x, y_positions[row], x + width, y_positions[row] + row_heights[row]),
                         panel, components, orientation, panel, 60, label_height, inner_padding)

    output_dir = project_root / OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    output_png = output_dir / OUTPUT_NAME
    canvas.save(output_png, dpi=(300, 300), optimize=True)

    svg_labels: list[tuple[str, int, int]] = []
    svg_placements: list[dict[str, object]] = []
    for panel, components, orientation, col, row, colspan in cell_specs:
        width = column_width * colspan + column_gap * (colspan - 1)
        x = margin + col * (column_width + column_gap)
        y = y_positions[row]
        svg_labels.append((panel, x + inner_padding, y))
        svg_placements.extend(component_layouts(
            (x, y, x + width, y + row_heights[row]), panel, components,
            orientation, 60, label_height, inner_padding
        ))
    output_svg = output_dir / OUTPUT_SVG_NAME
    write_editable_svg(output_svg, canvas.width, canvas.height, svg_labels, svg_placements)

    manifest_rows = [
        ["2a", "2a-i", "Paired sampling design", "Figure_01_PairedSamplingDesign", "a1"],
        ["2a", "2a-ii", "Paired Core-Edge PCA", "Figure_02_PairedPCA_CoreEdge", "a2"],
        ["2b", "2b", "MA plot", "Figure_04_MA", "b"],
        ["2c", "2c", "Volcano plot", "Figure_05_Volcano_PDE2A_Polished", "c"],
        ["2d", "2d", "Top DEG heatmap through PDE2A", "Figure_06_TopDEG_to_PDE2A_CoreEdgeGroupedHeatmap", "d"],
        ["2e", "2e", "Hallmark GSEA", "Figure_07_GSEA_Hallmark", "e"],
        ["2f", "2f-i", "GO Biological Process GSEA", "Figure_08_GSEA_GO_BP", "f1"],
        ["2f", "2f-ii", "Reactome GSEA", "Figure_09_GSEA_Reactome", "f2"],
        ["2g", "2g", "Leading-edge gene-module convergence", "Figure_12_GeneModule_LeadingEdgeHeatmap", "g"],
        ["2h", "2h", "PDE2A paired expression", "Figure_03_PDE2A_PairedExpression", "h"],
        ["2i", "2i", "LOPO candidate stability", "Figure_18_LOPO_Stability_Top20", "i"],
        ["2j", "2j", "Pareto candidate context", "Figure_21_Pareto_Global", "j"],
    ]
    source_file = output_dir / "source_data.tsv"
    with source_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["panel", "component", "description", "source_figure", "source_png", "source_data", "plotting_script", "layout_rule", "svg_object_id", "svg_x", "svg_y", "svg_width", "svg_height"])
        placement_by_component = {str(item["component"]): item for item in svg_placements}
        for panel, component, description, source_figure, source_key in manifest_rows:
            source_png, source_data, plotting_script, _ = source_paths(project_root, panel, source_key)
            layout_rule = "horizontal within 2a" if panel == "2a" else "horizontal within 2f" if panel == "2f" else "single panel"
            placement = placement_by_component[source_key]
            writer.writerow([
                panel, component, description, source_figure, source_png, source_data, plotting_script,
                layout_rule, f"panel_{placement['panel']}_{placement['component']}",
                placement["x"], placement["y"], placement["width"], placement["height"]
            ])

    print(f"Created {output_png}")
    print(f"Created {output_svg}")
    print(f"Created {source_file}")
    print(f"Canvas: {canvas.width} x {canvas.height} px at 300 dpi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
