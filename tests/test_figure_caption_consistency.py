"""Guard against panel, region-name, and coordinate mismatches in captions."""

from pathlib import Path
import ast
import json


ROOT = Path(__file__).resolve().parents[1]


def text(path):
    return (ROOT / path).read_text(encoding="utf-8")


def literal_assignment(path, name):
    tree = ast.parse(text(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
                return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {path}")


def test_figure1_box_names_bounds_and_caption_are_synchronized():
    folder = ROOT / "figures/Fig01_global_trends"
    regions = json.loads((folder / "Fig01_boxes.json").read_text(encoding="utf-8"))
    caption = (folder / "Fig01_caption.txt").read_text(encoding="utf-8")

    expected = [
        (1, "Arabian Sea and northwestern Indian Ocean", 60, 80, 15, 30),
        (2, "western North Atlantic", -70, -40, 30, 45),
        (3, "southern tropical Indian Ocean", 40, 80, -20, 0),
        (4, "western North Pacific east of Japan", 140, 170, 35, 45),
    ]
    actual = [
        (
            region["box"],
            region["name"],
            region["lon_min"],
            region["lon_max"],
            region["lat_min"],
            region["lat_max"],
        )
        for region in regions
    ]
    assert actual == expected
    for region in regions:
        assert -180 <= region["lon_min"] < region["lon_max"] <= 180
        assert -90 <= region["lat_min"] < region["lat_max"] <= 90
        phrase = (
            f'{region["name"]} (box {region["box"]}; '
            f'{region["caption_coordinates"]})'
        )
        assert caption.count(phrase) == 1

    script = text("figures/Fig01_global_trends/plot_Fig01.py")
    assert 'HERE / "Fig01_boxes.json"' in script
    assert "centres =" not in script


def test_main_figure_basin_panels_match_captions():
    fig3 = literal_assignment(
        "figures/Fig03_basin_intensity_trends/plot_Fig03.py", "B"
    )
    cap3 = text("figures/Fig03_basin_intensity_trends/Fig03_caption.txt")
    fig3_caption_names = {
        "North Indian": "North Indian Ocean",
        "South Indian": "South Indian Ocean",
    }
    for _, panel, name in fig3:
        assert f"{fig3_caption_names.get(name, name)} ({panel})" in cap3

    fig4 = literal_assignment("figures/Fig04_basin_decomposition/plot_Fig04.py", "B")
    cap4 = text("figures/Fig04_basin_decomposition/Fig04_caption.txt")
    fig4_caption_names = {
        "North Indian": "North Indian Ocean",
        "South Indian": "South Indian Ocean",
    }
    for _, panel, name, _ in fig4:
        assert f"({panel}) {fig4_caption_names.get(name, name)}" in cap4


def test_hemisphere_and_record_type_panel_descriptions_match_code():
    cap2 = text("figures/Fig02_latitude_decomposition/Fig02_caption.txt")
    assert "(a) the Northern Hemisphere" in cap2
    assert "(b) the Southern Hemisphere" in cap2

    panels3 = literal_assignment("figures/FigS03_poleward_shifts/plot_FigS03.py", "P")
    caps3 = text("figures/FigS03_poleward_shifts/FigS03_caption.txt")
    labels3 = {
        ("track", "NH"): "Northern Hemisphere TC-track records",
        ("track", "SH"): "Southern Hemisphere TC-track records",
        ("event", "NH"): "Northern Hemisphere TCEP events",
        ("event", "SH"): "Southern Hemisphere TCEP events",
    }
    for kind, hemisphere, panel, _ in panels3:
        phrase = f"{labels3[(kind, hemisphere)]} ({panel})"
        assert phrase in caps3

    caps4 = text("figures/FigS04_intensity_distributions/FigS04_caption.txt")
    assert "Northern Hemisphere (a) and Southern Hemisphere (b)" in caps4

    caps5 = text("figures/FigS05_land_ocean/FigS05_caption.txt")
    assert "Panels (a,b)" in caps5
    assert "Northern and Southern Hemispheres, respectively" in caps5
    assert "Panels (c,d)" in caps5
    script5 = text("figures/FigS05_land_ocean/plot_FigS05.py")
    assert 'for j, h in enumerate(["NH", "SH"])' in script5
    assert 'heading(ax, "ab"[j]' in script5
    assert 'heading(ax, "cd"[j]' in script5
