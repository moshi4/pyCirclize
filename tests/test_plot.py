import random
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from pycirclize import Circos
from pycirclize.parser import Genbank, Gff
from pycirclize.utils import (
    ColorCycler,
    load_eukaryote_example_dataset,
    load_prokaryote_example_file,
)

np.random.seed(0)
random.seed(0)


###########################################################
# Circos Class Plot
###########################################################


def test_circos_axis_plot(fig_outfile: Path):
    """Test `circos.axis()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    circos.axis()

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_circos_text_plot(fig_outfile: Path):
    """Test `circos.text()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    circos.text("center")
    circos.text("top", r=100)
    circos.text("right", r=100, deg=90)
    circos.text("right-middle", r=50, deg=90)
    circos.text("bottom", r=100, deg=180)
    circos.text("left", r=100, deg=270)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_circos_line_plot(fig_outfile: Path):
    """Test `circos.line()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    circos.line(r=100)
    circos.line(r=80, deg_lim=(0, 270), color="red")
    circos.line(r=60, deg_lim=(90, 360), color="blue", lw=3, ls="dotted")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_circos_rect_plot(fig_outfile: Path):
    """Test `circos.rect()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    circos.rect((80, 100))
    circos.rect((60, 80), deg_lim=(0, 270), fc="tomato")
    circos.rect((30, 50), deg_lim=(90, 360), fc="lime", ec="grey", lw=2, hatch="//")
    circos.rect((30, 100), deg_lim=(0, 90), fc="orange", alpha=0.2)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_circos_link_plot(fig_outfile: Path):
    """Test `circos.link()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    name2color = {"A": "red", "B": "blue", "C": "green"}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track = sector.add_track((95, 100))
        track.axis(fc=name2color[sector.name])

    # Plot links in various styles
    circos.link(("A", 0, 1), ("A", 7, 8))
    circos.link(("A", 1, 2), ("A", 7, 6), color="skyblue")
    circos.link(("A", 9, 10), ("B", 4, 3), direction=1, color="tomato")
    circos.link(("B", 5, 7), ("C", 6, 8), direction=1, ec="black", lw=1, hatch="//")
    circos.link(("B", 18, 16), ("B", 11, 13), r1=90, r2=90, color="violet", ec="red")
    circos.link(("C", 1, 3), ("B", 2, 0), direction=1, color="limegreen")
    circos.link(("C", 11.5, 14), ("A", 4, 3), direction=2, color="olive", ec="black")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_circos_colorbar_plot(fig_outfile: Path):
    """Test `circos.colorbar()`"""
    circos = Circos(sectors=dict(data=100), start=90)
    circos.axis()

    # Plot colorbar in various style
    vmin1, vmax1 = 0, 100
    circos.colorbar(vmin=vmin1, vmax=vmax1)
    circos.colorbar(bounds=(0.7, 0.6, 0.02, 0.3), vmin=vmin1, vmax=vmax1, cmap="bwr")

    vmin2, vmax2 = -200, 200
    circos.colorbar(
        bounds=(0.8, 0.6, 0.02, 0.3), vmin=vmin2, vmax=vmax2, cmap="viridis"
    )
    circos.colorbar(
        bounds=(0.3, 0.485, 0.4, 0.03),
        vmin=vmin2,
        vmax=vmax2,
        cmap="viridis",
        orientation="horizontal",
        colorbar_kws=dict(label="Colorbar in center"),
        tick_kws=dict(labelsize=12, colors="red"),
    )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_chord_diagram_plot(fig_outfile: Path, tsv_matrix_file: pd.DataFrame):
    """Test chord diagram plot"""
    circos = Circos.initialize_from_matrix(tsv_matrix_file)
    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_cytoband_plot(fig_outfile: Path, hg38_testdata_dir: Path):
    """Test hg38 cytoband plot"""
    # Add tracks for cytoband plot
    chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset(
        "hg38", cache_dir=hg38_testdata_dir
    )
    circos = Circos.initialize_from_bed(chr_bed_file, space=2)
    circos.add_cytoband_tracks((95, 100), cytoband_file)

    # Plot and check fig file exists
    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


###########################################################
# Sector Class Plot
###########################################################


def test_sector_axis_plot(fig_outfile: Path):
    """Test `sector.axis()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    sector_a = circos.get_sector("A")
    sector_a.axis()
    sector_b = circos.get_sector("B")
    sector_b.axis(fc="lightgrey", lw=2, ls="dotted")
    sector_c = circos.get_sector("C")
    sector_c.axis(fc="tomato", ec="blue", hatch="//")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_sector_text_plot(fig_outfile: Path):
    """Test `sector.text()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    name2color = {"A": "red", "B": "blue", "C": "green"}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        sector.axis()
        sector.text(sector.name)
        color = name2color[sector.name]
        sector.text(f"Center\n{sector.name}", r=50, size=12, color=color)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_sector_line_plot(fig_outfile: Path):
    """Test `sector_line()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        sector.axis()
        sector.line(r=90)
        sector_center = (sector.start + sector.end) / 2
        sector.line(r=80, end=sector_center, color="red")
        sector.line(r=80, start=sector_center, color="blue")
        sector.line(r=60, color="green", lw=2, ls="dotted")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_sector_rect_plot(fig_outfile: Path):
    """Test `sector.rect()`"""
    ColorCycler.set_cmap("tab10")
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        sector.axis()
        sector.rect(r_lim=(90, 100), fc="tomato")
        sector_center = (sector.start + sector.end) / 2
        sector.rect(end=sector_center, r_lim=(70, 80), color="skyblue")
        sector.rect(start=sector_center, r_lim=(70, 80), color="limegreen")
        for i in range(int(sector.size)):
            sector.rect(i, i + 1, (50, 60), fc=ColorCycler(), ec="black", lw=1)
        start, end = sector.start + 3, sector.end - 3
        sector.rect(start, end, (30, 100), color="orange", alpha=0.2)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


###########################################################
# Track Class Plot
###########################################################


def test_track_axis_plot(fig_outfile: Path):
    """Test `track.axis()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track1 = sector.add_track((90, 100))
        track1.axis()
        track2 = sector.add_track((70, 80))
        track2.axis(fc="lightgrey")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_text_plot(fig_outfile: Path):
    """Test `track.text()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track1 = sector.add_track((90, 100))
        track1.axis()
        track1.text(sector.name)
        track2 = sector.add_track((70, 80))
        track2.axis(fc="lightgrey")
        track2.text(sector.name, orientation="vertical", color="red", size=15)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_rect_plot(fig_outfile: Path):
    """Test `track.rect()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track1 = sector.add_track((90, 100))
        track1.axis()
        # Plot rect & text (style1)
        for i in range(int(track1.size)):
            start, end = i, i + 1
            track1.rect(start, end, fc=ColorCycler())
            track1.text(str(end), (end + start) / 2)
        # Plot rect & text (style2)
        track2 = sector.add_track((70, 80))
        for i in range(int(track2.size)):
            start, end = i, i + 1
            track2.rect(start, end, fc=ColorCycler(), ec="white", lw=1)
            track2.text(
                str(end), (end + start) / 2, color="white", orientation="vertical"
            )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_arrow_plot(fig_outfile: Path):
    """Test `track.arrow()`"""
    ColorCycler.set_cmap("tab10")
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        sector.axis()
        # Plot forward arrow with default style
        track1 = sector.add_track((90, 100))
        for i in range(int(track1.size)):
            start, end = i, i + 1
            track1.arrow(start, end, fc=ColorCycler())
        # Plot reverse arrow with user-specified style
        track2 = sector.add_track((70, 80))
        for i in range(int(track2.size)):
            start, end = i, i + 1
            track2.arrow(end, start, head_length=4, shaft_ratio=1.0, fc=ColorCycler())

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_xticks_plot(fig_outfile: Path):
    """Test `track.xticks()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track1 = sector.add_track((90, 100))
        track1.axis()
        # Plot outer xticks
        pos_list = list(range(0, int(track1.size) + 1))
        labels = [f"{i:02d}" for i in pos_list]
        track1.xticks(pos_list, labels)
        # Plot inner xticks label
        labels = [f"Label{i:02d}" for i in pos_list]
        track1.xticks(
            pos_list,
            labels,
            outer=False,
            tick_length=0,
            label_margin=2,
            label_orientation="vertical",
        )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_xticks_by_interval_plot(fig_outfile: Path):
    """Test `track.xticks_by_interval()`"""
    sectors = {"A": 10000000, "B": 20000000, "C": 15000000}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        # Major & Minor xticks
        track1 = sector.add_track((90, 100))
        track1.axis()
        mb_size = 1000000
        track1.xticks_by_interval(mb_size, label_orientation="vertical")
        track1.xticks_by_interval(mb_size / 5, tick_length=1, show_label=False)
        # Mb formatted xticks
        track2 = sector.add_track((80, 90))
        track2.xticks_by_interval(
            mb_size,
            outer=False,
            show_bottom_line=True,
            label_orientation="vertical",
            label_formatter=lambda v: f"{v / mb_size:.1f} Mb",
        )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_yticks_plot(fig_outfile: Path):
    """Test `track.yticks()`"""
    sectors = {"A": 10000000, "B": 20000000, "C": 15000000}
    circos = Circos(sectors, space=15)
    for sector in circos.sectors:
        # Plot yticks
        track1 = sector.add_track((80, 100))
        track1.axis()
        y = [0, 5, 10, 15, 20]
        y_labels = list(map(str, y))
        track1.yticks(y, y_labels)
        # Plto yticks label
        track2 = sector.add_track((50, 70), r_pad_ratio=0.1)
        track2.axis()
        y = [10, 15, 20]
        y_labels = ["Label1", "Label2", "Label3"]
        track2.yticks(
            y,
            y_labels,
            vmin=10,
            vmax=25,
            side="left",
            line_kws=dict(color="red", lw=1),
            text_kws=dict(color="blue"),
        )

    circos.savefig(fig_outfile)
    fig_outfile.exists()


def test_track_grid_plot(fig_outfile: Path):
    """Test `track.grid()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        # Plot Y-axis grid line (Default: 6 grid line)
        track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
        track1.axis()
        track1.grid()
        # Plot X-axis grid line (interval=1)
        track2 = sector.add_track((55, 75))
        track2.axis()
        track2.grid(y_grid_num=None, x_grid_interval=1, color="red")
        # Plot both XY-axis grid line
        track3 = sector.add_track((30, 50))
        track3.grid(y_grid_num=11, x_grid_interval=0.5, color="blue", ls="dashed")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_line_plot(fig_outfile: Path):
    """Test `track.line()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track = sector.add_track((80, 100), r_pad_ratio=0.1)
        track.axis()
        track.xticks_by_interval(1)
        vmin, vmax = 0, 10
        # Line between start-end two points
        track.line([track.start, track.end], [vmin, vmax], lw=1.5, ls="dotted")
        # Line of random value points
        x = np.linspace(track.start, track.end, int(track.size) * 5 + 1)
        y = np.random.randint(vmin, vmax, len(x))
        track.line(x, y)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_scatter_plot(fig_outfile: Path):
    """Test `track.scatter()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        track = sector.add_track((80, 100), r_pad_ratio=0.1)
        track.axis()
        track.xticks_by_interval(1)
        vmin, vmax = 0, 10
        x = np.linspace(track.start, track.end, int(track.size) * 5 + 1)
        y = np.random.randint(vmin, vmax, len(x))
        track.scatter(x, y)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_bar_plot(fig_outfile: Path):
    """Test `track.bar()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        vmin, vmax = 1, 10
        x = np.linspace(sector.start + 0.5, sector.end - 0.5, int(sector.size))
        y = np.random.randint(vmin, vmax, len(x))
        # Plot bar (default)
        track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
        track1.axis()
        track1.xticks_by_interval(1)
        track1.xticks_by_interval(0.1, tick_length=1, show_label=False)
        track1.bar(x, y)
        # Plot stacked bar with user-specified params
        track2 = sector.add_track((50, 70))
        track2.axis()
        track2.xticks_by_interval(1, outer=False)

        ColorCycler.set_cmap("tab10")
        tab10_colors = [ColorCycler() for _ in range(len(x))]
        track2.bar(
            x, y, width=1.0, color=tab10_colors, ec="grey", lw=0.5, vmax=vmax * 2
        )

        ColorCycler.set_cmap("Pastel1")
        pastel_colors = [ColorCycler() for _ in range(len(x))]
        y2 = np.random.randint(vmin, vmax, len(x))
        track2.bar(
            x,
            y2,
            width=1.0,
            bottom=y,
            color=pastel_colors,
            ec="grey",
            lw=0.5,
            hatch="//",
            vmax=vmax * 2,
        )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_fill_between_plot(fig_outfile: Path):
    """Test `track.fill_between()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=5)
    for sector in circos.sectors:
        vmin, vmax = 0, 10
        # Plot fill_between with simple lines
        track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
        track1.axis()
        track1.xticks_by_interval(1)
        track1.fill_between(
            x=[track1.start, track1.end], y1=[vmin, vmax], y2=[vmin, vmax / 2]
        )
        # Plot fill_between with random points line
        track2 = sector.add_track((50, 70), r_pad_ratio=0.1)
        track2.axis()
        x = np.linspace(track2.start, track2.end, int(track2.size) * 5 + 1)
        y = np.random.randint(vmin, vmax, len(x))
        track2.fill_between(x, y, ec="black", lw=1, ls="dashed")

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_heatmap_plot(fig_outfile: Path):
    """Test `track.heatmap()`"""
    sectors = {"A": 10, "B": 20, "C": 15}
    circos = Circos(sectors, space=10)
    for sector in circos.sectors:
        # Plot heatmap
        track1 = sector.add_track((80, 100))
        track1.axis()
        track1.xticks_by_interval(1)
        data = np.random.randint(0, 10, (4, int(sector.size)))
        track1.heatmap(data, show_value=True)
        # Plot heatmap with labels
        track2 = sector.add_track((50, 70))
        track2.axis()
        x = np.linspace(1, int(track2.size), int(track2.size)) - 0.5
        xlabels = [str(int(v + 1)) for v in x]
        track2.xticks(x, xlabels, outer=False)
        track2.yticks([0.5, 1.5, 2.5, 3.5, 4.5], list("ABCDE"), vmin=0, vmax=5)
        data = np.random.randint(0, 100, (5, int(sector.size)))
        track2.heatmap(data, cmap="viridis", rect_kws=dict(ec="white", lw=1))

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_tree_plot(fig_outfile: Path):
    """Test `track.heatmap()`"""
    # Load newick tree
    tree_text = "((((A:1,B:1)100:1,(C:1,D:1)100:1)100:1,(E:2,F:2)90:1):1,G:6)100;"
    tree = Phylo.read(StringIO(tree_text), "newick")
    # Initialize circos sector by tree size
    circos = Circos(sectors={"Tree": tree.count_terminals()})
    sector = circos.sectors[0]
    # Plot tree
    track = sector.add_track((50, 100))
    track.axis(ec="lightgrey")
    track.tree(tree, leaf_label_size=12)

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_genomic_features_genbank_plot(
    fig_outfile: Path,
    prokaryote_testdata_dir: Path,
):
    """Test `track.genomic_features()` with genbank file"""
    # Load Genbank file
    gbk_file = load_prokaryote_example_file(
        "enterobacteria_phage.gbk",
        cache_dir=prokaryote_testdata_dir,
    )
    gbk = Genbank(gbk_file)
    # Initialize circos sector by genome size
    circos = Circos(sectors={gbk.name: gbk.range_size})
    circos.text("Enterobacteria phage\n(NC_000902)", size=15)
    sector = circos.sectors[0]
    # Outer track
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(5000, label_formatter=lambda v: f"{v / 1000:.0f} Kb")
    outer_track.xticks_by_interval(1000, tick_length=1, show_label=False)
    # Plot forward & reverse CDS genomic features
    cds_track = sector.add_track((90, 95))
    cds_track.genomic_features(
        gbk.extract_features("CDS", 1), plotstyle="arrow", fc="salmon"
    )
    cds_track.genomic_features(
        gbk.extract_features("CDS", -1), plotstyle="arrow", fc="skyblue"
    )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_track_genomic_features_gff_plot(
    fig_outfile: Path,
    prokaryote_testdata_dir: Path,
):
    """Test `track.genomic_features()` with gff file"""
    # Load Genbank file
    gff_file = load_prokaryote_example_file(
        "enterobacteria_phage.gff",
        cache_dir=prokaryote_testdata_dir,
    )
    gff = Gff(gff_file)
    # Initialize circos sector by genome size
    circos = Circos(sectors={gff.name: gff.range_size})
    circos.text("Enterobacteria phage\n(NC_000902)", size=15)
    sector = circos.sectors[0]
    # Outer track
    outer_track = sector.add_track((98, 100))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(5000, label_formatter=lambda v: f"{v / 1000:.0f} Kb")
    outer_track.xticks_by_interval(1000, tick_length=1, show_label=False)
    # Plot forward & reverse CDS genomic features
    cds_track = sector.add_track((90, 95))
    cds_track.genomic_features(
        gff.extract_features("CDS", 1), plotstyle="arrow", fc="salmon"
    )
    cds_track.genomic_features(
        gff.extract_features("CDS", -1), plotstyle="arrow", fc="skyblue"
    )

    circos.savefig(fig_outfile)
    assert fig_outfile.exists()
