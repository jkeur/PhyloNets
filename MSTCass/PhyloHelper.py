import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo


def _fix_node_names(tree):
    """
    Fix the node names as follows:
    - Give all anonymous nodes a name (ID = index)
    - In each node name, replace '#' by '_' (for reticulations).

    :param tree: A tree (of type BaseTree.Tree)
    :return: The tree with relabeled nodes
    """
    for i, clade in enumerate(tree.find_clades()):
        # Give all anonymous nodes a name
        if not clade.name:
            clade.name = '_' + str(i)
        else:
            clade.name = clade.name.replace('#', '_')


def read(file: str):
    """
    Read the input file and give each node a name
    :param file: A file name
    :return: The trees
    """
    trees = list(Phylo.parse(file, "newick"))
    for t in trees:
        _fix_node_names(t)
    print(f"SUCCESS: '{file}' has been read successfully.")

    return trees


def tree2newick(tree):
    """
    Get a Newick formatted string of the tree
    :param tree: A tree
    :return: The tree in the Newick format
    """
    items = []
    for clade in tree.root.clades:
        s = ''
        if len(clade.clades) > 0:
            subtree = tree2newick(clade)
            if subtree != '':
                s += subtree
        s += '' if clade.name is None else clade.name
        items.append(s)
    return '(' + ','.join(items) + ')'


def get_term_names(trees: list):
    """
    Get the sorted union of the leaf names in the given trees
    :param trees: A list of trees (of type BaseTree.Tree)
    :return: A list of terminal names
    """
    names = set()
    node_list = []
    for t in trees:
        for clade in t.find_clades():
            if clade.name not in node_list:
                node_list.append(clade.name)
                if clade.is_terminal():
                    names.add(clade.name)
    return sorted(names)


def get_terminals(tree):
    """
    Get the terminals in the tree, ordered by name
    :param tree: A tree (of type Phylo.BaseTree.Tree)
    :return: A sorted list
    """
    terms = tree.get_terminals()
    return sorted(terms, key=lambda t: t.name)


def get_clusters(trees: list):
    """
    Get the union of all clusters in the trees.
    Don't return the cluster on all taxa
    :param trees: A list of trees (of type BaseTree.Tree)
    :return: A list of clusters
    """
    clusters = []
    taxa = set()
    tree_ind = []
    for i, tree in enumerate(trees):
        t_cl, _taxa = get_tree_clusters(tree)
        taxa.update(_taxa)
        for cl in t_cl:
            if cl not in clusters:
                clusters.append(cl)
                tree_ind.append([i])
            else:
                tree_ind[clusters.index(cl)].append(i)
    return clusters, tree_ind, sorted(taxa)


def get_tree_clusters(tree):
    """ Get a list of clusters, where each cluster is a list of taxa, AND a list of the taxa
    Don't return the cluster on all taxa
    """
    clusters = []
    taxa = set()
    for clade in tree.find_clades(terminal=False):
        cl = []
        cl_to_visit = [clade]
        for subclade in cl_to_visit:
            for child in subclade:
                if child.name[0] == '_':
                    cl_to_visit.append(child)
                else:
                    cl.append(child.name)
        cl.sort()
        clusters.append(cl)
        taxa.update(cl)
    return clusters, sorted(taxa)


def draw_trees(trees, title=None, filename=None):
    """
    Draw the trees in one figure
    :param trees: A list of trees
    :param title: A title
    :param filename: A file name to save the figure
    """
    n_figs = len(trees)  # #figures
    if n_figs == 0:  # If there is nothing to show, then return
        print("INFO: There is nothing to show.")
        return
    plt.clf()
    n_cols = int(np.ceil(np.sqrt(n_figs)))
    n_rows = n_cols - 1 if n_figs <= n_cols * (n_cols - 1) else n_cols
    fig, ax = plt.subplots(n_rows, n_cols, num=1)
    for i, t in enumerate(trees):
        if n_cols == 1:
            axes = ax
        else:
            ix = np.unravel_index(i, ax.shape)
            axes = ax[ix]
        Phylo.draw(t, do_show=False, axes=axes)
    if title is not None:
        fig.suptitle(title)
    if filename is not None:
        import os
        if not os.path.isdir('output/'):
            os.makedirs('output/')
        plt.savefig('output/' + filename, bbox_inches='tight')  # Save the figure
    plt.show()
