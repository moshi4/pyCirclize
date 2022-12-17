from __future__ import annotations

from collections import Counter
from copy import deepcopy
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree


class TreeUtil:
    """Tree Utility Class"""

    @staticmethod
    def load_tree(
        treedata: str | Path | StringIO | Tree,
        format: str = "newick",
    ) -> Tree:
        """Load tree

        Parameters
        ----------
        treedata : str | Path | StringIO | Tree
            Phylogenetic tree data (`File-like object` or `Tree object`)
        format : str, optional
            Tree format (e.g. `newick`, `nexus`, `phyloxml`, ...)

        Returns
        -------
        tree : Tree
            Tree object
        """
        if isinstance(treedata, Tree):
            tree = treedata
        else:
            tree: Tree = Phylo.read(file=treedata, format=format)
        return tree

    @staticmethod
    def set_unique_node_name(tree: Tree) -> Tree:
        """Set unique node name (N_1, N_2, ..., N_XXX)

        Parameters
        ----------
        tree : Tree
            Tree object

        Returns
        -------
        tree: Tree
            Unique node name set tree object
        """
        tree = deepcopy(tree)
        node: Clade
        for idx, node in enumerate(tree.get_nonterminals(), 1):
            node.name = f"N_{idx}" if node.name is None else node.name
        return tree

    @staticmethod
    def set_node_color(
        tree: Tree, node_color_list: list[tuple[list[str], str]]
    ) -> Tree:
        """Set node color

        Parameters
        ----------
        tree : Tree
            Tree
        node_color_list : list[tuple[list[str], str]]
            Node color list. If multi nodes are set, MRCA of target nodes is used.

        Returns
        -------
        Tree
            Node color set tree object
        """
        tree = deepcopy(tree)
        # Get mrca(most recent common ancestor) node & color dict
        mrca_node_name2color = {}
        for nodes, color in node_color_list:
            mrca_node = tree.common_ancestor(nodes)
            mrca_node_name2color[mrca_node.name] = color
        # Set node color
        node1: Clade
        node2: Clade
        for node1 in tree.find_clades():
            if node1.name in mrca_node_name2color:
                for node2 in tree.find_clades():
                    if node1.is_parent_of(node2):
                        # node2.color = 'color_name' cannot accept mpl color value
                        node2._color = mrca_node_name2color[node1.name]
        return tree

    @staticmethod
    def check_node_name_dup(tree: Tree):
        """Check node name duplication

        Parameters
        ----------
        tree : Tree
            Tree object
        """
        node_names = [n.name for n in tree.find_clades()]
        err_msg = ""
        for node_name, count in Counter(node_names).items():
            if count > 1:
                err_msg += f"{node_name=} is duplicated in tree ({count=})\n"
        if err_msg != "":
            raise ValueError("\n" + err_msg)

    @staticmethod
    def to_ultrametric_tree(tree: Tree) -> Tree:
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
        name2depth = {n.name: d for n, d in tree.depths(True).items()}
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
