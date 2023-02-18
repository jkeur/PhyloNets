import string

import networkx as nx
from matplotlib import pyplot as plt

from ClusterSet import ClusterSet, Cluster
from PhyloNet import PhyloNet


class IncGraph(nx.Graph):
    """ Incompatibility Graph """

    def __init__(self, cs: ClusterSet = None):
        """
        Create the incompatibility graph IG (using NetworkX).
        Edge (c1, c2) is in IG iff cluster c1 and c2 are incompatible.

        Example
        -------
        >>> cs = ClusterSet([['A', 'B'], ['B', 'C']])
        >>> IG = IncGraph(cs)
        >>> print(IG)
        Graph with 2 nodes and 1 edges
        :return: The incompatibility graph
        """
        super(IncGraph, self).__init__()
        if cs is None:
            self.X = None
            return
        self.X = list(cs.X)
        clusters = list(cs.C)
        # Print.var(cs.get_mst_sets())
        # TODO: Simplify the labels
        cl_labels = []
        for c in clusters:
            cl_labels.append(','.join(c))
        self.add_nodes_from(cl_labels)
        for i, c1 in enumerate(clusters):
            c1 = Cluster(c1)
            for j in range(i + 1, len(clusters)):
                c2 = clusters[j]

                # If clusters c1 and c2 are incompatible, add an edge between them
                if not c1.iscompatible(c2):
                    self.add_edge(','.join(c1), ','.join(c2))

    def subgraph(self, nodes):
        nodes = list(nodes)
        node_labels = nodes
        if nodes and isinstance(nodes[0], list):
            # Convert the clusters to labels
            node_labels = []
            for c in nodes:
                node_labels.append(','.join(c))
        g_sub = super(IncGraph, self).subgraph(node_labels)
        g_sub.X = self.X
        return g_sub

    def get_nontrivial_cc(self):
        """
        Get a list of the non-trivial connected components, i.e. the non-singleton components
        :return: A list (of type [ClusterSet])
        """
        ccs = []
        for cc in nx.connected_components(self):
            if len(cc) <= 1:
                continue
            C = []
            for c in cc:
                C.append(c.split(','))
            cs = ClusterSet(C, self.X)
            cs._remove_unused_leaves()
            ccs.append(cs)
        return ccs

    def draw(self, title='Incompatibility Graph', axes=None, filename=None, exclude_singletons=True, node_color=None):
        """
        Draw the incompatibility graph
        :param title: A title
        :param axes: axes
        :param filename: A filename to store the figure
        :param exclude_singletons: Exclude singletons T/F
        """

        if axes is None:
            plt.clf()

        # If singletons should be excluded, remove them
        if exclude_singletons:
            iso = list(nx.isolates(self))
            if len(iso) >= 1:
                self.remove_nodes_from(iso)

        # To color the edges corresponding to the incompatibilities, get a list of incompatibilities
        clusters = [c.split(',') for c in self.nodes]
        cs = ClusterSet(clusters, self.X)
        incs = cs.get_incompatibilities()
        e_color_map = []
        for u, v in self.edges:
            inc = list(set(u.split(',')).intersection(v.split(',')))
            inc.sort()
            e_color_map.append(incs.index(inc))

        labels = None
        if node_color:
            # Node colors
            n_colors = []
            for node in self.nodes():
                if node in node_color.keys():
                    n_colors.append(node_color[node])
                else:
                    n_colors.append('lightGrey')
            node_color = n_colors
        else:
            node_color = 'lightGrey'
        # Create a nice network layout
        nx.draw_circular(self, with_labels=True, edge_color=e_color_map, node_color=node_color, ax=axes, labels=labels)
        e_labels = nx.get_edge_attributes(self, 'weight')
        if len(e_labels) > 0:
            # nx.draw_networkx_edge_labels(net, pos, e_labels, ax=axes)
            nx.draw_networkx_edge_labels(self, e_labels, ax=axes)
        # Set the title
        if title is not None:
            plt.gcf().suptitle(title)
        if axes is None:
            plt.show()
        if filename is not None:
            import os
            if not os.path.isdir('output/'):
                os.makedirs('output/')
            plt.savefig('output/' + filename, bbox_inches='tight')
