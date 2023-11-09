from __future__ import annotations

import io
import os
from collections import Counter, defaultdict
from copy import deepcopy
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import urlparse
from urllib.request import urlopen

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from matplotlib.patches import Rectangle

from pycirclize import utils

if TYPE_CHECKING:
    from pycirclize.track import Track


class TreeViz:
    """Phylogenetic Tree Visualization Class

    Interface for changing tree properties and adding tree annotations in a track
    """

    def __init__(
        self,
        tree_data: str | Path | Tree,  # type: ignore
        *,
        format: str = "newick",
        outer: bool = True,
        align_leaf_label: bool = True,
        ignore_branch_length: bool = False,
        leaf_label_size: float = 12,
        leaf_label_rmargin: float = 2.0,
        reverse: bool = False,
        ladderize: bool = False,
        line_kws: dict[str, Any] | None = None,
        align_line_kws: dict[str, Any] | None = None,
        track: Track,
    ):
        """
        Parameters
        ----------
        tree_data : str | Path | Tree
            Tree data (`File`|`File URL`|`Tree Object`|`Tree String`)
        format : str, optional
            Tree format (`newick`|`phyloxml`|`nexus`|`nexml`|`cdao`)
        outer : bool, optional
            If True, plot tree on outer side. If False, plot tree on inner side.
        align_leaf_label: bool, optional
            If True, align leaf label.
        ignore_branch_length : bool, optional
            If True, ignore branch length for plotting tree.
        leaf_label_size : float, optional
            Leaf label size
        leaf_label_rmargin : float, optional
            Leaf label radius margin
        reverse : bool, optional
            If True, reverse tree
        ladderize : bool, optional
            If True, ladderize tree
        line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(color="red", lw=1, ls="dashed", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        align_line_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(lw=1, ls="dotted", alpha=1.0, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        track : Track
            Track for tree visualization
        """
        tree = self.load_tree(tree_data, format=format)

        # Set unique node name and branch length if not exists
        tree, _ = self._set_uniq_innode_name(tree)
        self._check_node_name_dup(tree)
        max_tree_depth = max(tree.depths().values())
        if ignore_branch_length or max_tree_depth == 0:
            tree = self._to_ultrametric_tree(tree)
        if ladderize:
            tree.ladderize()
        if reverse:
            for clade in tree.find_clades():
                clade.clades = clade.clades[::-1]
        self._tree = tree

        # Set plot parameters
        self._outer = outer
        self._align_leaf_label = align_leaf_label
        self._leaf_label_size = leaf_label_size
        self._leaf_label_rmargin = leaf_label_rmargin
        self._track = track

        line_kws = {} if line_kws is None else line_kws
        line_kws.setdefault("color", "black")
        self._line_kws = line_kws

        align_line_kws = {} if align_line_kws is None else align_line_kws
        align_line_kws.setdefault("ls", "dashed")
        align_line_kws.setdefault("alpha", 0.5)
        self._align_line_kws = align_line_kws

        self._node2label_props: dict[str, dict[str, Any]] = defaultdict(lambda: {})
        self._node2line_props: dict[str, dict[str, Any]] = defaultdict(lambda: {})

    ############################################################
    # Properties
    ############################################################

    @property
    def track(self) -> Track:
        """Track for tree visualization"""
        return self._track

    @property
    def tree(self) -> Tree:
        """BioPython's Tree Object"""
        return self._tree

    @cached_property
    def leaf_labels(self) -> list[str]:
        """Leaf labels"""
        return [str(n.name) for n in self.tree.get_terminals()]

    @cached_property
    def innode_labels(self) -> list[str]:
        """Internal node labels"""
        return [str(n.name) for n in self.tree.get_nonterminals()]

    @cached_property
    def all_node_labels(self) -> list[str]:
        """All node labels"""
        return self.leaf_labels + self.innode_labels

    @cached_property
    def max_tree_depth(self) -> float:
        """Max tree depth (root -> leaf max branch length)"""
        return max(self.tree.depths().values())

    @cached_property
    def name2xr(self) -> dict[str, tuple[float, float]]:
        """Tree node name & node xr coordinate dict"""
        return self._calc_name2xr()

    @cached_property
    def name2rect(self) -> dict[str, Rectangle]:
        """Tree node name & rectangle dict"""
        return self._calc_name2rect()

    ############################################################
    # Public Method
    ############################################################

    @staticmethod
    def load_tree(data: str | Path | Tree, format: str) -> Tree:
        """Load tree data

        Parameters
        ----------
        data : str | Path | Tree
            Tree data
        format : str
            Tree format

        Returns
        -------
        tree : Tree
            Tree object
        """
        if isinstance(data, str) and urlparse(data).scheme in ("http", "https"):
            # Load tree file from URL
            return Phylo.read(io.StringIO(urlopen(data).read().decode()), format=format)
        elif isinstance(data, (str, Path)) and os.path.isfile(data):
            # Load tree file
            return Phylo.read(data, format=format)
        elif isinstance(data, str):
            # Load tree string
            return Phylo.read(io.StringIO(data), format=format)
        elif isinstance(data, Tree):
            return data
        else:
            raise ValueError(f"{data=} is invalid input tree data!!")

    def highlight(
        self,
        query: str | list[str] | tuple[str],
        *,
        color: str,
        alpha: float = 0.5,
        **kwargs,
    ) -> None:
        """Plot highlight for target node

        Parameters
        ----------
        query : str | list[str] | tuple[str]
            Search query node name(s) for highlight. If multiple node names are set,
            MRCA(Most Recent Common Ancester) node is set.
        color : str | None, optional
            Highlight color
        alpha : float, optional
            Highlight color alpha(transparancy) value
        **kwargs : dict, optional
            Patch properties (e.g. `ec="blue", lw=1.0, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Get target rectangle for highlight
        target_node_name = self._search_target_node_name(query)
        rect = self.name2rect[target_node_name]

        # Set default kwargs
        kwargs.setdefault("zorder", 0)

        # Setup track.rect() parameters
        start, end = rect.get_x(), rect.get_x() + rect.get_width()
        r_lim = (rect.get_y(), rect.get_y() + rect.get_height())

        self.track.rect(start, end, r_lim=r_lim, color=color, alpha=alpha, **kwargs)

    def marker(
        self,
        query: str | list[str] | tuple[str],
        *,
        marker: str = "o",
        size: int = 6,
        descendent: bool = True,
        **kwargs,
    ) -> None:
        """Plot marker on target node(s)

        Parameters
        ----------
        query : str | list[str] | tuple[str]
            Search query node name(s) for plotting marker.
            If multiple node names are set,
            MRCA(Most Recent Common Ancester) node is set.
        marker : str, optional
            Marker type (e.g. `o`, `s`, `D`, `P`, `*`, `x`, `d`, `^`, `v`, `<`, `>`)
            <https://matplotlib.org/stable/api/markers_api.html>
        size : int, optional
            Marker size
        descendent : bool, optional
            If True, plot markers on target node's descendent as well.
        **kwargs : dict, optional
            Axes.scatter properties (e.g. `fc="lime", ec="black", lw=0.5, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.scatter.html>
        """
        target_node_name = self._search_target_node_name(query)

        x: list[float] = []
        r: list[float] = []
        rmin, rmax = self.track.r_plot_lim
        if descendent:
            clade: Clade = next(self.tree.find_clades(target_node_name))
            descendent_nodes: list[Clade] = list(clade.find_clades())
            for descendent_node in descendent_nodes:
                node_x, node_r = self.name2xr[str(descendent_node.name)]
                if descendent_node.is_terminal() and self._align_leaf_label:
                    node_r = rmax if self._outer else rmin
                x.append(node_x)
                r.append(node_r)
        else:
            node_x, node_r = self.name2xr[target_node_name]
            target_node: Clade = next(self.tree.find_clades(target_node_name))
            if target_node.is_terminal() and self._align_leaf_label:
                node_r = rmax if self._outer else rmin
            x.append(node_x)
            r.append(node_r)

        self.track.scatter(
            x, r, s=size**2, vmin=rmin, vmax=rmax, marker=marker, **kwargs
        )

    def set_node_label_props(self, target_node_label: str, **kwargs) -> None:
        """Set tree node label properties

        Parameters
        ----------
        target_node_label : str
            Target node label name
        kwargs : dict, optional
            Text properties (e.g. `color="red", size=15, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        self._search_target_node_name(target_node_label)
        self._node2label_props[target_node_label].update(kwargs)

    def set_node_line_props(
        self,
        query: str | list[str] | tuple[str],
        *,
        descendent: bool = True,
        apply_label_color: bool = False,
        **kwargs,
    ) -> None:
        """Set tree node line properties

        Parameters
        ----------
        query : str | list[str] | tuple[str]
            Search query node name(s) for coloring tree node line.
            If multiple node names are set,
            MRCA(Most Recent Common Ancester) node is set.
        descendent : bool, optional
            If True, set properties on target node's descendent as well.
        apply_label_color : bool, optional
            If True & `descendent=True` & kwargs contain color keyword,
            apply node line color to node label color as well.
        **kwargs : dict, optional
            Patch properties (e.g. `color="blue", lw=2.0, ls="dashed", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
        """
        target_node_name = self._search_target_node_name(query)

        clade: Clade = next(self.tree.find_clades(target_node_name))
        if descendent:
            descendent_nodes: list[Clade] = list(clade.find_clades())
            for descendent_node in descendent_nodes:
                node_name = str(descendent_node.name)
                self._node2line_props[node_name] = kwargs
                if apply_label_color and "color" in kwargs:
                    self._node2label_props[node_name].update(color=kwargs["color"])
        else:
            self._node2line_props[str(clade.name)] = kwargs

    ############################################################
    # Private Method
    ############################################################

    def _set_uniq_innode_name(self, tree: Tree) -> tuple[Tree, list[str]]:
        """Set unique internal node name (N_1, N_2, ..., N_XXX)

        Parameters
        ----------
        tree : Tree
            Tree object

        Returns
        -------
        tree, uniq_node_names: tuple[Tree, list[str]]
            Unique node name set tree object & set unique node names
        """
        tree = deepcopy(tree)
        uniq_innode_names: list[str] = []
        for idx, node in enumerate(tree.get_nonterminals(), 1):
            uniq_innode_name = f"N_{idx}"
            if node.name is None:
                node.name = uniq_innode_name
                uniq_innode_names.append(uniq_innode_name)
        return tree, uniq_innode_names

    def _to_ultrametric_tree(self, tree: Tree) -> Tree:
        """Convert to ultrametric tree

        Parameters
        ----------
        tree : Tree
            Tree

        Returns
        -------
        tree : Tree
            Ultrametric tree
        """
        tree = deepcopy(tree)
        # Get unit branch depth info
        name2depth = {str(n.name): float(d) for n, d in tree.depths(True).items()}
        name2depth = dict(sorted(name2depth.items(), key=lambda t: t[1], reverse=True))
        max_tree_depth = max(name2depth.values())
        # Reset node branch length
        for node in tree.find_clades():
            node.branch_length = None
        tree.root.branch_length = 0
        # Calculate appropriate ultrametric tree branch length
        for name, depth in name2depth.items():
            node = next(tree.find_clades(name))
            if not node.is_terminal():
                continue
            path: list[Clade] | None = tree.get_path(node)
            if path is None:
                raise ValueError(f"{name=} node not exists?")
            if depth == max_tree_depth:
                for path_node in path:
                    path_node.branch_length = 1
            else:
                # Collect nodes info which has branch length
                bl_sum, bl_exist_node_count = 0, 0
                for path_node in path:
                    if path_node.branch_length is not None:
                        bl_sum += path_node.branch_length
                        bl_exist_node_count += 1
                # Set branch length to no branch length nodes
                other_bl = (max_tree_depth - bl_sum) / (len(path) - bl_exist_node_count)
                for path_node in path:
                    if path_node.branch_length is None:
                        path_node.branch_length = other_bl
        return tree

    def _check_node_name_dup(self, tree: Tree) -> None:
        """Check node name duplication in tree

        Parameters
        ----------
        tree : Tree
            Tree object
        """
        all_node_names = [str(n.name) for n in tree.find_clades()]
        err_msg = ""
        for node_name, count in Counter(all_node_names).items():
            if count > 1:
                err_msg += f"{node_name=} is duplicated in tree ({count=}).\n"
        if err_msg != "":
            err_msg += "\nTreeViz cannot handle tree with duplicate node names!!"
            raise ValueError("\n" + err_msg)

    def _search_target_node_name(
        self,
        query: str | list[str] | tuple[str],
    ) -> str:
        """Search target node name from query

        Parameters
        ----------
        query : str | list[str] | tuple[str]
            Search query node name(s). If multiple node names are set,
            MRCA(Most Recent Common Ancester) node is set.

        Returns
        -------
        target_node_name : str
            Target node name
        """
        self._check_node_name_exist(query)
        if isinstance(query, (list, tuple)):
            target_node_name = self.tree.common_ancestor(*query).name
        else:
            target_node_name = query
        return target_node_name

    def _check_node_name_exist(
        self,
        query: str | list[str] | tuple[str],
    ) -> None:
        """Check node name exist in tree

        Parameters
        ----------
        query : str | list[str] | tuple[str]
            Query node name(s) for checking exist
        """
        if isinstance(query, str):
            query = [query]
        err_msg = ""
        for node_name in query:
            if node_name not in self.all_node_labels:
                err_msg += f"{node_name=} is not found in tree.\n"
        if err_msg != "":
            err_msg = f"\n{err_msg}\nAvailable node names:\n{self.all_node_labels}"
            raise ValueError(err_msg)

    def _calc_name2xr(self) -> dict[str, tuple[float, float]]:
        """Calculate node name & xr coordinate

        Returns
        -------
        name2xr : dict[str, tuple[float, float]]
            Tree node name & xr coordinate
        """
        track = self.track
        name2depth = {str(n.name): float(d) for n, d in self.tree.depths().items()}
        # Calculate x, r unit size of depth
        x_unit_size = track.size / self.tree.count_terminals()
        r_unit_size = track.r_plot_size / self.max_tree_depth
        # Calculate leaf node (x, r) coordinate
        name2xr: dict[str, tuple[float, float]] = {}
        node: Clade
        for idx, node in enumerate(self.tree.get_terminals()):
            x = track.start + (x_unit_size * idx) + (x_unit_size / 2)
            if self._outer:
                r = min(track.r_plot_lim) + r_unit_size * name2depth[str(node.name)]
            else:
                r = max(track.r_plot_lim) - r_unit_size * name2depth[str(node.name)]
            name2xr[str(node.name)] = (x, r)
        # Calculate internal node (x, r) coordinate
        for node in self.tree.get_nonterminals(order="postorder"):
            x = sum([name2xr[n.name][0] for n in node.clades]) / len(node.clades)
            if self._outer:
                r = min(track.r_plot_lim) + r_unit_size * name2depth[str(node.name)]
            else:
                r = max(track.r_plot_lim) - r_unit_size * name2depth[str(node.name)]
            name2xr[str(node.name)] = (x, r)
        return name2xr

    def _calc_name2rect(self) -> dict[str, Rectangle]:
        """Calculate tree node name & rectangle

        Returns
        -------
        name2rect : dict[str, Rectangle]
            Tree node name & rectangle dict
        """
        name2rect: dict[str, Rectangle] = {}
        for name, xr in self.name2xr.items():
            # Get parent node
            node: Clade = next(self.tree.find_clades(name))
            if node == self.tree.root:
                parent_node = node
            else:
                tree_path = self.tree.get_path(node.name)
                tree_path = [self.tree.root] + tree_path  # type: ignore
                parent_node: Clade = tree_path[-2]

            # Get child node xr coordinates
            child_node_names = [str(n.name) for n in node.find_clades()]
            x_list: list[float] = []
            r_list: list[float] = []
            for child_node_name in child_node_names:
                x, r = self.name2xr[child_node_name]
                x_list.append(x)
                r_list.append(r)

            # Calculate rectangle min-max xr coordinate
            x_unit_size = self.track.size / len(self.leaf_labels)
            xmin = min(x_list) - (x_unit_size / 2)
            xmax = max(x_list) + (x_unit_size / 2)

            parent_xr = self.name2xr[str(parent_node.name)]
            upper_r = (xr[1] + parent_xr[1]) / 2
            r_plot_lim = self.track.r_plot_lim
            if self._align_leaf_label:
                lower_r = max(r_plot_lim) if self._outer else min(r_plot_lim)
            else:
                lower_r = max(r_list) if self._outer else min(r_list)
            rmin, rmax = min(upper_r, lower_r), max(upper_r, lower_r)

            # Set rectangle
            rect = Rectangle(
                xy=(xmin, rmin),
                width=xmax - xmin,
                height=rmax - rmin,
            )
            name2rect[name] = rect

        return name2rect

    def _plot_tree_line(self) -> None:
        """Plot tree line"""
        # Plot tree line by node (x, r) coordinate
        for node in self.tree.get_nonterminals():
            parent_x, parent_r = self.name2xr[node.name]
            child_node: Clade
            for child_node in node.clades:
                child_x, child_r = self.name2xr[str(child_node.name)]
                # Set node color if exists
                _line_kws = deepcopy(self._line_kws)
                _line_kws.update(self._node2line_props[str(child_node.name)])
                # Plot horizontal line
                h_line_points = (parent_x, child_x), (parent_r, parent_r)
                self.track._simpleline(*h_line_points, **_line_kws)
                # Plot vertical line
                v_line_points = (child_x, child_x), (parent_r, child_r)
                self.track._simpleline(*v_line_points, **_line_kws)
                # Plot vertical line for leaf node alignment
                if self._align_leaf_label and child_node.is_terminal():
                    r_plot_lim = self.track.r_plot_lim
                    end_r = max(r_plot_lim) if self._outer else min(r_plot_lim)
                    v_align_line_points = (child_x, child_x), (child_r, end_r)
                    _align_line_kws = deepcopy(self._align_line_kws)
                    _align_line_kws.update(dict(color=_line_kws["color"]))
                    self.track._simpleline(*v_align_line_points, **_align_line_kws)

    def _plot_tree_label(self) -> None:
        """Plot tree label"""
        text_kws = dict(size=self._leaf_label_size, orientation="vertical")
        for node in self.tree.get_terminals():
            if self._leaf_label_size <= 0:
                continue
            # Set label text (x, r) position
            label = str(node.name)
            x, r = self.name2xr[label]
            if self._align_leaf_label:
                r_plot_lim = self.track.r_plot_lim
                r = max(r_plot_lim) if self._outer else min(r_plot_lim)
            rmargin = self._leaf_label_rmargin
            r = r + rmargin if self._outer else r - rmargin

            # Set label text properties
            _text_kws = deepcopy(text_kws)
            rad = self.track.x_to_rad(x)
            params = utils.plot.get_label_params_by_rad(rad, "vertical", self._outer)
            _text_kws.update(params)
            _text_kws.update(self._node2label_props[label])

            self.track.text(label, x, r, **_text_kws)  # type: ignore
