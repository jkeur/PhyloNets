import os
import string
from inspect import stack
from itertools import combinations

import networkx as nx
from matplotlib import pyplot as plt
from networkx.drawing.nx_pydot import graphviz_layout

LIVE_DRAW_ENABLED = False  # Disable the function draw for live renderings
DRAW_ENABLED = True  # Disable the function draw, but not draw_multi


def config(enable_live_draw=False, enable_draw=True):
    global LIVE_DRAW_ENABLED, DRAW_ENABLED
    LIVE_DRAW_ENABLED = enable_live_draw
    DRAW_ENABLED = enable_draw


class PhyloNet(nx.DiGraph):
    """
    Helper class for a phylogenetic network (of instance NetworkX.DiGraph)
    """

    def __new__(cls, di_graph=nx.DiGraph()):
        """ Create a PhyloNet instance. The first argument can be used to define either a DiGraph or MultiDiGraph. """
        if not isinstance(di_graph, nx.DiGraph):
            print(di_graph, type(di_graph), isinstance(di_graph, nx.DiGraph))
        #     raise TypeError('PhyloNet should be of type DiGraph (or MultiDiGraph).')
        return di_graph.__new__(cls)

    def __init__(self, type1=None):
        super().__init__(type1)
        self._z_edges = None

    def copy(self, as_view=False):
        return super().copy(as_view)

    def get_root(self):
        """
        Get the root of the network
        :return: The unique root
        """
        for n in self.nodes:
            if self.in_degree(n) == 0:
                return n

    def get_ret_num(self):
        """
        Get the number of binary hybrid nodes/reticulations
        :return: The number of reticulations
        """
        return sum(self.in_degree(n) - 1 for n in self.nodes if self.in_degree(n) >= 2)

    def get_level(self):
        """ Get the level of the network, i.e. the maximal number of reticulations per bi-connected component. """
        g = self.to_undirected()
        level = 0
        for cc in nx.biconnected_components(g):
            if len(cc) <= 2:
                continue
            cc_ret_num = len([n for n in cc if self.in_degree(n) >= 2])
            level = max(level, cc_ret_num)
        return level

    def get_leaves(self):
        """
        For the given acyclic digraph, get the leaf nodes
        :return: A sorted list of leaves
        """
        return sorted(n for n in self.nodes() if self.out_degree(n) == 0)

    def get_dummy_leaves(self):
        """
        Get the dummy leaves
        :return: A list of names
        """
        return sorted(n for n in self.nodes() if len(n) >= 2 and n[0:2] == "d_")

    def get_hardwired_cluster(self, source: str, filter_name=False):
        """
        Get the hardwired cluster below the given source node
        :param source: A node
        :return: The hardwired cluster below the given source/node
        """
        if source not in self.nodes:
            return None
        # If the source is a leaf
        if source[0] != '_':
            return [source]
        # assert source[0] == '_':
        # The source is an internal node
        if filter_name:
            return sorted([n for n in nx.descendants(self, source) if n[0] != '_'])
        return sorted([n for n in nx.descendants(self, source) if self.out_degree(n) == 0])

    # NICE TODO: Implement a more efficient method
    def _get_soft_wired_clusters(self, source):
        """
        Get all soft-wired clusters
        :param source: A node
        :return: A set of the soft-wired clusters below the given source/node
        """

        def _select_ret_inputs(ret_ins: dict):
            """
            In the network net, select one incoming edge per reticulation to leave a tree
            :param ret_ins: A dictionary:
                Each (r: i) determines that the i-th (starting from 0) edge is selected for reticulation r.
            :return: The resulting tree
            """
            _tree = self.copy()
            for r, i in ret_ins.items():
                E_in = list(_tree.in_edges(r))
                E_in.__delitem__(i)
                _tree.remove_edges_from(E_in)
            for n in list(_tree.nodes):
                if _tree.out_degree(n) == 0 and n[0] == '_':
                    _tree.remove_node(n)
            # Remove singletons
            _tree.remove_nodes_from(list(nx.isolates(_tree)))
            return _tree

        # Find all reticulations below the source
        rets = [n for n in nx.descendants(self, source) if self.in_degree(n) >= 2]

        # If there are no reticulations below, then return the hardwired cluster
        if len(rets) == 0:
            c = self.get_hardwired_cluster(source)
            if c:
                return [c]
            return set()

        # Now find all soft-wired clusters
        c_strings = set()
        done = False
        ret_inputs = {r: 0 for r in rets}
        while not done:
            # For each reticulation, let 1 incoming edge remain
            tree = _select_ret_inputs(ret_inputs)
            c = tree.get_hardwired_cluster(source, filter_name=True)
            if c:
                c_strings.add(','.join(c))
            # Select the next tree configuration
            for r in rets:
                if ret_inputs[r] + 1 < self.in_degree(r):
                    ret_inputs[r] += 1
                    break
                elif r == rets[-1]:
                    done = True
                else:
                    ret_inputs[r] = 0
        clusters = []
        for c_s in c_strings:
            clusters.append(c_s.split(','))
        return clusters

    def get_st_sets(self, source: str):
        """
        Given that this network is a tree, return all ST-sets from a given source node and below it
        :param source: A source node
        :return: A list of lists (ST-sets)
        """
        # If the source is a leaf
        if source[0] != '_':
            return source
        sts = [self.get_hardwired_cluster(source)]
        childs = list(self.successors(source))
        # For each combination of children, add the hardwired cluster as ST-set
        for k in range(2, len(childs)):
            # Choose k children
            for combi in combinations(childs, k):
                # Get the hardwired cluster for this combination
                cl = []
                for n in combi:
                    cl.extend(self.get_hardwired_cluster(n))  # Singletons will be included
                sts.append(sorted(cl))
        # Recursively search for the ST-sets below the child split nodes
        for n in childs:
            if n[0] == '_':  # If n is a split node
                sts.extend(self.get_st_sets(n))
            else:
                sts.append([n])
        return sts

    def represents_old(self, C: set):
        """
        Check if the given network represents all clusters in C
        :param C: A set of clusters
        :param print_instr: Print instructions T/F
        :return: Does the PN represent the clusters [T/F]?
        """
        # Make a copy of the clusters C and remove each cluster that is represented. Note: This is inefficient!
        _C = C.copy()
        for node in self.nodes:
            if self.out_degree(node) == 0:
                continue
            clusters = self._get_soft_wired_clusters(node)
            for c in clusters:
                if c in _C:
                    _C.remove(c)

        return len(_C) == 0

    def represents(self, cs, x_added, lcas, print_instr=False):
        """
        Does this network represent the cluster defined in ClusterSet cs? T/F
        :param cs: A ClusterSet
        :param x_added:
        :param print_instr: Print instructions T/F
        :return: T/F
        """
        return self.represents_old(cs.C)

    def represents_cl(self, cluster, lcas):
        """ Check if the given cluster is represented """
        cl = list(cluster)
        if len(cl) == 1:
            return True
        c_s = ','.join(cluster)
        if c_s not in lcas.keys():
            # Check if the cluster is in the network
            print(f'ERROR: I do not see cluster {c_s}. Is it there?')
            self.draw(True, '?', all_labels=True)
        for ca in lcas[c_s]:  # For each common ancestor
            # Get the hardwired cluster
            hw_cl = self.get_hardwired_cluster(ca)
            if hw_cl == cluster:
                return True
        # Create a copy of the network & remove the leaves that are not in the cluster
        net = self.copy()
        # Remove each (parent) node that is not part of the cluster
        X_not = set(self.get_leaves()).difference(cluster)
        X_not_p = {list(net.predecessors(x))[0] for x in X_not if net.in_degree(x) == 1}
        net.remove_nodes_from(X_not)
        net.clean()
        while X_not_p:
            X_not = X_not_p
            X_not_1p = [x for x in X_not if
                        self.in_degree(x) == 1 and self.out_degree(x) <= 2]  # All x in X_not having 1 parent
            X_not_p = [list(self.predecessors(x))[0] for x in X_not_1p]
            net.remove_nodes_from(X_not)
            net.clean()
        for ca in lcas[c_s]:  # For each common ancestor
            # Get the hardwired cluster
            hw_cl = net.get_hardwired_cluster(ca)
            if not hw_cl:
                continue
            if hw_cl == cl:
                return True

        return False

    def get_lcas_2(self, n1: str, n2: str):
        """ Get the LCA of the 2 given nodes """
        if n1[0] != '_':
            # n1 is a leaf
            n1 = self.predecessors(n1).__next__()
        if n2[0] != '_':
            # n2 is a leaf
            n2 = self.predecessors(n2).__next__()
        if n2 in nx.ancestors(self, n1):
            lcas = set()
            lcas.add(n2)
            return lcas
        if n1 in nx.ancestors(self, n2):
            lcas = set()
            lcas.add(n1)
            return lcas
        # n1 is not an ancestor of n2, nor vice versa
        P1 = [n1]
        P2 = [n2]
        lcas = set()
        net = self.copy()
        while True:
            P1_new = set()
            for n in P1:
                if n not in lcas:
                    P1_new.update(net.predecessors(n))
            if not P1_new:
                break
            P1 = P1_new
            for p in P1:
                # If node 2 is reachable from p, then p is an LCA
                if n2 in nx.descendants(net, p):
                    lcas.add(p)
                    net.remove_node(p)
            P2_new = set()
            for n in P2:
                if n not in lcas:
                    P2_new.update(net.predecessors(n))
            if not P2_new:
                break
            P2 = P2_new
            for p in P2:
                # If node 1 is reachable from p, then p is an LCA
                if n1 in nx.descendants(net, p):
                    lcas.add(p)
                    net.remove_node(p)
        if lcas:
            return lcas
        self.draw(True, f'LCA({n1}, {n2}) = ?', all_labels=True)

    def _subdivide_edge(self, e):
        """
        Insert node v into edge e in the given network
        :param e: An edge
        :return: The added node (name)
        """
        # Let v be the new node
        v = self.new_node_name()
        # Remove the old edge
        self.remove_edge(e[0], e[1])
        # Add the new edges
        self.add_edge(e[0], v)
        self.add_edge(v, e[1])
        return v

    def new_node_name(self) -> str:
        """
        Get a name for a new (internal) node (splitting or reticulation) to insert in the given network
        :return: A new node name
        """
        node_nrs = [int(n[1:]) for n in self.nodes if n[1:].isnumeric() and n[0] == '_']
        if not node_nrs:
            return '_1'
        return '_' + str(max(node_nrs) + 1)

    def add_new_root(self):
        """
        Add a new root as parent of the current root and return the labels of the new and old root, respectively.

        :return: The labels of the new and old root, respectively
        """
        root = self.get_root()
        label = self.new_node_name()
        self.add_edge(label, root)
        return label, root

    def add_ret_below_2(self, parents: list, x: str):
        """
        Add a reticulation below two of the given parent nodes
        If no parents are given, then insert a reticulation instead of two dummy leaves, and a dummy leaf below it
        :return: A list of networks having a new reticulation
        """
        # If dummy leaves are given, then search their parents
        if parents[0][0:2] == 'd_':
            # Search for the dummy leaves & their parents
            p = {}
            for dum in parents:
                p[dum] = list(self.predecessors(dum))[0]
            parents = p
        else:
            parents = dict.fromkeys(parents, None)
        parent_combis = combinations(parents, 2)
        nets = []
        # t becomes the new reticulation/hybrid node
        t = self.new_node_name()
        for ret_parents in parent_combis:
            net = self.copy()
            # If dummy leaves are selected, then hang the reticulation below its parent & remove the dummy leaf
            p = parents[ret_parents[0]]
            if p:
                # Preserve one dummy root. Subdivide its incoming edge and make it a reticulation.
                net.remove_node(ret_parents[1])
                # If there is another dummy leaf in the network, then add a dummy leaf below the new reticulation
                if len(parents) >= 3:
                    v = net._subdivide_edge((p, ret_parents[0]))
                    net.add_edge(parents[ret_parents[1]], v)
                else:
                    # #parents = 2 => Add leaf x below reticulation t
                    net.remove_node(ret_parents[0])
                    net.add_edge(p, t)
                    net.add_edge(parents[ret_parents[1]], t)
                    net.add_edge(t, x)
            else:
                # Add the edges going to t
                net.add_edges_from((v, t) for v in ret_parents)
                # Add dummy leaves
                parents_new = list(parents)
                parents_new.remove(ret_parents[0])
                parents_new.remove(ret_parents[1])
                parents_new.append(t)
                net.add_edges_from((d, f'd_{i}') for i, d in enumerate(parents_new))
            nets.append(net)
        return nets

    def add_leaf_below_ret(self, x: str, edges, cs, stay_binary: bool = False, print_info=False):
        """
        Add a leaf below a reticulation below the given edges as follows.
        Create a reticulation t, a leaf l labeled x and an edge from t to l. Then, for each edge e_i, i = 1,2,... in the
        list, insert a node v_i into e_i and add an edge from v_i to t
        :param x: A leaf
        :param cs: A set of clusters (ClusterSet)
        :param edges: A list of edges
        :param stay_binary: Let the returned network(s) be binary
        :param print_info:
        :return: The new PN
        """
        parent_nodes = []
        E_done = []
        for e in edges:
            # Note: it can happen that an edge occurs multiple times in E
            if e in E_done:
                parent_nodes.append(parent_nodes[E_done.index(e)])
                continue
            # If e represents a node, i.e. e = (u, u)
            if e[0] == e[1]:
                parent_nodes.append(e[0])
            else:
                parent_nodes.append(self._subdivide_edge(e))
            E_done.append(e)
        # If the new leaf (or ST-set) is to be hung below exactly 1 edge (or node)
        if len(edges) == 1:
            # Hang the leaf below this edge (or node)
            self.add_edge(parent_nodes[0], x)
        else:
            # If leaf x should be hung below 3 edges or more and if the network should stay binary, then add 1 ret. now
            if len(parent_nodes) > 2 and stay_binary:
                # Print.var(parent_nodes, print_info)
                # Choose two nodes to create a reticulation below, and a dummy leaf below it
                return self.add_ret_below_2(parent_nodes, x)

            # t becomes the new reticulation/hybrid node
            t = self.new_node_name()
            # Add the edges going to t
            self.add_edges_from((v, t) for v in parent_nodes)
            # Add edge (t, x)
            self.add_edge(t, x)

        return [self]

    def get_e_out(self, cs_new, cs_in, x: str, lcas):
        """
        Get the lowest edge below which a new reticulation should be hung, such that all clusters below the edge will be
        represented in the newly formed network
        :param cs_new: The new ClusterSet
        :param cs_in: The minimal clusters containing the ST-set
        :param x: A leaf
        :return:
        """
        # If cs_in is empty, then return None
        if not cs_in:
            raise RuntimeError('cs_in should not be empty')

        # If the cluster X\{x} is in C, then E_out = the root edge
        Xm = list(cs_new.X)
        Xm.remove(x)
        if Xm in cs_new.C:
            r = self.get_root()
            return cs_new.X, list(self.out_edges(r))[0]

        # Check if the cluster X_in\st exists as (not necessarily strict) subset of any cluster in C.
        # Make a list (set) of all such edges.
        c_out = set(cs_in.get_used_taxa())
        c_out.remove(x)
        E_out = set()
        cl_out = set()
        # for c in cs_new.maximal().C:
        for c in cs_new.C:
            # If cluster c contains c_out and not leaf x
            if c_out.issubset(c) and x not in c:
                cl_out.update(c)

                # Cluster c in C contains c_out. Hence, the ST-set should be hung below an edge that is above cluster c.
                c_s = ','.join(c)
                if c_s in lcas.keys():
                    lca = lcas[c_s]
                    e = list(self.in_edges(lca))[0]
                else:
                    e = self.get_cluster_edge(c)
                if not e:
                    # Cluster c is not represented
                    print(f'ERROR: Cluster {c} is not represented')
                    return None, None
                E_out.add(e)

        # Check the length of E_out and return if possible
        _l = len(E_out)
        if _l == 1:
            return cl_out, list(E_out)[0]
        if _l == 0:
            return None, None

        # Find the highest E
        cl_nodes = {v for u, v in E_out}
        lca = self.get_comm_anc(cl_nodes)
        # if len(lcas) > 2:
        #     Print.var(cl_nodes)
        #     Print.var(lcas)
        #     self.draw(True, f'LCA({cl_nodes}) ?= {lcas}', all_labels=True)
        # lca = lcas[0]
        e_in = self.in_edges(lca)
        if len(e_in) == 0:
            # This node has no incoming edge. It must be the root
            e_out = self.out_edges(lca)
            return cl_out, list(e_out)[0]
        return cl_out, list(e_in)[0]

    def get_cluster_edge(self, c: list):
        """
        Find the lowest edge that represents the cluster
        :param c: A cluster
        :return: An edge
        """
        # Get the cluster leaves
        if len(c) == 1:
            comm_anc = c[0]
        else:
            # Quickly check if the root edge is present
            r = self.get_root()
            assert r not in c

            # Get the (lowest) common ancestor and then its incoming edge
            comm_anc = self.get_comm_anc(c)
            if not comm_anc:
                # This cluster is not represented
                return None
        e_in = self.in_edges(comm_anc)
        if len(e_in) == 0:
            # This node has no incoming edge. It must be the root
            return None, comm_anc
        return list(e_in)[0]

    def get_comm_anc(self, nodes):
        """
        Get the common ancestor of the given nodes
        :param nodes: A list of nodes
        :return: A node
        """

        assert nodes
        nodes = list(nodes)
        # Print.var(nodes, True)
        # If all leaves in the network are given, then LCA = the child of the root
        if nodes == self.get_leaves():
            r = self.get_root()
            return list(self.out_edges(r))[0]
        net = self.copy()
        # If the given nodes are leaves
        if net.out_degree(nodes[0]) == 0:
            # Remove each (parent) node that is not part of the cluster
            X_not = set(self.get_leaves()).difference(nodes)
            X_not_p = {list(net.predecessors(x))[0] for x in X_not}
            while X_not_p:
                X_not = X_not_p
                X_not_1p = [x for x in X_not if
                            self.in_degree(x) == 1 and self.out_degree(x) <= 2]  # All x in X_not having 1 parent
                X_not_p = [list(self.predecessors(x))[0] for x in X_not_1p]
                net.remove_nodes_from(X_not)
                net.clean()
            # net.draw(True, 'OK?', all_labels=True)
        comm_anc = None
        while comm_anc is None:  # TODO: will this always work well?
            # Search the common ancestor for each pair of nodes
            new_comm_ancs = set()
            for combi in combinations(nodes, 2):
                lcas = net.get_lcas_2(combi[0], combi[1])
                if not lcas:
                    print(f'ERROR: LCA({combi[0]}, {combi[1]}) = ?')
                    self.draw(True, f'LCA({combi[0]}, {combi[1]}) = ?', all_labels=True)
                # lcas = sorted(lcas)
                new_comm_ancs.update(lcas)
            if not new_comm_ancs:
                # Check if for some combination, a wrong CA is taken
                print(comm_ancs)
                self.draw(True, f'Check CA({nodes})', all_labels=True)
            comm_ancs = new_comm_ancs
            if len(comm_ancs) == 0:
                return None
            if len(comm_ancs) == 1:
                comm_anc = comm_ancs.pop()
            else:
                # If they are not all the same, then search for the common ancestor of these nodes!
                nodes = comm_ancs
        return comm_anc

    def _get_ret_ret_in_edges(self, node, incl_intermediate=True):
        """
        Suppose reticulation [node] has a reticulation v as parent, then return the union of the non-ret. parents
        of both [node] and v.
        """
        in_edges = []
        rets = [node]
        for r in rets:
            for e in self.in_edges(r):
                if self.in_degree(e[0]) > 1:
                    rets.append(e[0])
                    if incl_intermediate:
                        in_edges.append(e)
                else:
                    in_edges.append(e)
        return in_edges

    def get_hang_edges(self, cs_new, x: str, lcas, filter_edges: bool = True, guess_edges: bool = True,
                       stay_bin=False, full_cx=False, print_instr=False):
        """
        Get a list of lists of edges per cluster in Cin, to hang the ST-set below
        :param cs_new: A set of clusters (of type ClusterSet)
        :param x: A leaf to add
        :param lcas: The LCAs per cluster
        :param filter_edges: Filter the edges T/F
        :param guess_edges:
        :param stay_bin: Stay binary T/F
        :param full_cx: Use Cx as containing all clusters containing x
        :param print_instr: T/F
        :return: A list of edges
        """
        # Init
        leaves = self.get_leaves()
        if len(cs_new.X) != len(leaves) + 1:
            self.draw(True, 'Check #leaves', all_labels=True)
            return {}
        # st_str = BitString.from_bool(x in st for x in X)
        cl_edges = {}
        cl_parents = {}
        meta_candidates = {}  # Candidates for being a meta-node to be split
        cs_in = cs_new.get_cs_in(x, full_cx)

        # Find the edge e_out
        if full_cx:
            _cs_in = cs_in.minimal()
            cl_out, e_out = self.get_e_out(cs_new, _cs_in, x, lcas)
        else:
            cl_out, e_out = self.get_e_out(cs_new, cs_in, x, lcas)
        if e_out:
            E_cl = set()
            # E_cl.add(e_out)
            # Add the edges above e_out
            E_cl.update(self.out_edges(e_out[0]))
            c_s_out = ','.join(cl_out)  # Cluster as string
            cl_edges[c_s_out] = E_cl
            P = set()
            P.add(e_out[0])
            cl_parents[c_s_out] = P

        # For each cluster c in cs_in, get the lowest edge in the network that represents it.
        # Then, find all edges E_cl below which x may be hung.
        self.update_z_edges()
        for c in cs_in.C:
            if len(c) == len(cs_new.X):
                continue
            # cl = cl \ {x}
            c.remove(x)
            if filter_edges:
                E_cl, meta_candidate = self._get_cluster_ext_edges(c, cs_new, lcas, full_cx)
                # assert len(set(E_cl)) == len(E_cl)  # There should be no duplicate edges
                if not E_cl:  # If this cluster c is not represented, then return
                    return {}
            else:
                E_cl = list(self.edges)
                meta_candidate = None

            c_s = ','.join(c)  # Cluster as string
            cl_edges[c_s] = E_cl
            cl_parents[c_s] = set(e[0] for e in E_cl)
            if meta_candidate:
                meta_candidates[c_s] = meta_candidate
        # Do not use only minimal clusters in Cx, but all clusters containing x
        if full_cx:
            # Filter the edge sets: only preserve the minimal edge sets
            cl_edges_new = {}
            # Print.var(cl_edges, True)
            for c_s1, E1 in cl_edges.items():
                E1 = set(E1)
                c1 = set(c_s1.split(','))
                issupset = False
                for c_s2, E2 in reversed(cl_edges.items()):
                    if c_s2 == c_s1:
                        break
                    if E1.issuperset(E2):
                        issupset = True
                        break
                if 2 <= len(c1) < len(cs_new.X) - 1:
                    for _x in c1:
                        E1.update(self.in_edges(_x))
                if not issupset:
                    cl_edges_new[c_s1] = E1
            cl_edges = cl_edges_new
            # Update cl_parents
            for c_s in reversed(list(cl_parents.keys())):
                if c_s not in cl_edges.keys():
                    del cl_parents[c_s]
        if not guess_edges:
            return cl_edges

        for c_s1, e in meta_candidates.items():
            p = e[0]
            # Is this parent p a parent of an edge in another hang-edge group?
            # has_overlap = False
            for c_s2, P in cl_parents.items():
                if p in P:
                    E = list(filter(lambda _e: _e[0] == p, cl_edges[c_s2]))
                    assert E
                    e_div = E[0]
                    if e_div == e:
                        continue
                    # If v in e = (u, v) is not a reticulation, then continue
                    if self.in_degree(e_div[1]) == 1:
                        continue
                    # Regard p as meta node & split it
                    n_new = self._subdivide_edge(e_div)
                    # In order to minimize the number of reticulations (according to our Greedy approach, we should
                    # subdivide edge e_div and hang cluster (c_s1) below the new node n_new
                    self.remove_edge(e[0], e[1])
                    self.add_edge(n_new, e[1])
                    # Set the edges E_cl below which the new leaf should be hung
                    E_cl = set()
                    E_cl.add((e[0], n_new))
                    E_cl.add((n_new, e[1]))
                    E_cl.add((n_new, e_div[1]))
                    cl_edges[c_s1] = E_cl

                    # Add 0-edges
                    ret = e_div[1]
                    while ret:
                        succ = list(self.successors(ret))
                        if len(succ) == 0:
                            self.draw(True, f'Check ret {ret}', all_labels=True)
                            raise RuntimeError(f'Reticulation {ret} does not have any childs.')
                        desc = succ[0]
                        # If its child is a reticulation too, then add the 0-edge between the two rets.
                        if self.in_degree(desc) >= 2:
                            cl_edges[c_s1].add((ret, desc))
                            ret = desc
                        else:
                            ret = None

                    # Clean cl_edges
                    for _c_s in cl_edges.keys():
                        if e_div in cl_edges[_c_s]:
                            cl_edges[_c_s].remove(e_div)
                        if e in cl_edges[_c_s]:
                            cl_edges[_c_s].remove(e)

                    # For cluster 2, change edge e_div to (n_new, e_div[1])
                    if e_div not in cl_edges[c_s2]:
                        # In this case, there are multiple possibilities to split the metanode!
                        #  Which one to choose?
                        continue
                    cl_edges[c_s2].add((n_new, e_div[1]))
                    break
        return cl_edges

    def update_z_edges(self):
        """
        Get the zero-edges
        :return:
        """
        self._z_edges = set()
        for n in self.nodes:
            if self.in_degree(n) >= 2:  # If n is a ret.
                for e in self.in_edges(n):
                    self._z_edges.add(e)
                    p = e[0]
                    # If p has only child reticulations
                    onlyRets = True
                    for c in self.successors(p):
                        if self.in_degree(c) == 1:
                            onlyRets = False
                            break
                    if onlyRets:
                        self._z_edges.add(list(self.in_edges(p))[0])

    def get_z_edges_below(self, node: str):
        z_edges = set()
        for e in self.out_edges(node):
            if e in self._z_edges:
                z_edges.add(e)
                z_edges.update(self.get_z_edges_below(e[1]))
        return z_edges

    def _get_cluster_ext_edges(self, c: list, cs_new, lcas, full_cx=False):
        """
        Get the edges below which another cluster cl2 can be hung, such that cluster (cl1 union cl2) will be represented
        :param c: The cluster to search the 'extension edges' for
        :param cs_new: The set of clusters to be represented after hanging cluster cl back
        :return: A list of 'extension edges'
        """
        # Get the cluster edge
        if len(c) == 1:
            lca = c[0]
        else:
            c_str = ','.join(c)
            lca = lcas[c_str]
        e = list(self.in_edges(lca))[0]
        if not e:
            # Cluster c is not represented
            return [], None
        E_cl = []

        # If edge e is the root edge, then return
        if len(c) >= 2:
            # Add 0-edges connected to v
            E_cl.extend(self.get_z_edges_below(e[1]))
        # Add 0-edges connecting the subtree below the LCA with a reticulation
        for x in c:
            p = list(self.predecessors(x))[0]  # Parent of leaf x
            rest_only_z_outs = True
            for _e in self.out_edges(p):
                if _e[1] != x and _e not in self._z_edges:
                    rest_only_z_outs = False
                    break
            if rest_only_z_outs:
                for _e in self.out_edges(p):
                    if _e in self._z_edges:
                        if _e not in E_cl:
                            E_cl.append(_e)
                        for _ez in self.get_z_edges_below(_e[1]):
                            if _ez not in E_cl:
                                E_cl.append(_ez)
            else:
                E_cl.append((p, x))

        if self.out_degree(e[0]) == 1:  # If u has out-degree 1, then (u is a ret, so) add edge e
            E_cl.append(e)
            # Add all in-edges of the reticulation (and of its parent reticulations, etc.)
            E_cl.extend(self._get_ret_ret_in_edges(e[0]))
            return list(reversed(E_cl)), None
        # u of e = (u, v) is a split node.
        # If the other children of 'u' are only reticulations, then add (w, u)
        neighbors = list(self.neighbors(e[0]))

        # Consider the other children of u. Are they all reticulations?
        neighbors.remove(e[1])
        only_rets = True
        has_ret = False
        for n in neighbors:
            if self.in_degree(n) == 1:  # If node n is a splitting
                only_rets = False
                if has_ret:
                    break
            else:  # If node n is a reticulation
                has_ret = True
                if not only_rets:
                    break
        if only_rets and self.in_degree(e[0]) == 1:
            # # Add (w, u)
            E_cl.extend(self.in_edges(e[0]))
            # Add (u, u); let the new node be hung below node u!
            E_cl.append((e[0], e[0]))
            # The other children of node e[0] are all reticulations
            return E_cl, None
        # u is a splitting and there is a splitting among the neighbours of v
        if len(c) == 1:
            # If this cluster cl is a single leaf, then surely, x should be hung below this edge,
            # unless there is a ret. below e[0]
            if e not in E_cl:
                E_cl.append(e)
            if has_ret:
                # Node u has a child ret.
                return E_cl, e  # Node u of edge e=(u,v) is a meta-node candidate
            return E_cl, None
        # Otherwise, add e
        assert e[1][0] == '_', self.draw(True, f'Edge {e}')
        if e not in E_cl:
            if c in cs_new.C or full_cx:
                E_cl.append(e)
        if self.in_degree(e[1]) == 1 and not E_cl:
            E_cl.append((e[1], e[1]))
        return E_cl, None

    def clean(self, rem_dummy_root=False):
        """
        Clean the network from redundant edges, etc.
        :param rem_dummy_root: Remove the dummy root T/F
        :return: The cleaned PN
        """
        # If multiple edges go from v1 to v2, then remove redundant edges among them
        E = []
        for i, e in reversed(list(enumerate(self.edges()))):
            if not e in E:
                E.append(e)

        # If a node has both in- and out-degree 1, then remove the node and connect its ancestor with its descendant
        for node in list(self.nodes):
            if self.in_degree(node) == 1 and self.out_degree(node) == 1:
                e1 = list(self.in_edges(node))[0]
                e2 = list(self.out_edges(node))[0]
                self.remove_node(e1[1])
                self.add_edge(e1[0], e2[1])

        if rem_dummy_root:
            self.remove_dummy_root()

    def remove_dummy_root(self):
        r = self.get_root()
        if self.out_degree(r) == 1:
            self.remove_node(r)

    def merge(self, other: "PhyloNet"):
        X_other = other.get_leaves()
        # Search the node in this network to be replaced by the other network
        for n in self.nodes:
            if self.out_degree(n) <= 1:  # Only consider splittings
                continue
            # Node n is an internal node
            c = self.get_hardwired_cluster(n)
            if X_other == c:
                # Remove node n and its descendants
                nodes2rem = [n]
                nodes2rem.extend(nx.descendants(self, n))
                preds = list(self.predecessors(n))
                if len(preds) == 0:  # If n does not have a parent, then return the other tree!
                    self.remove_nodes_from(list(self.nodes))
                    self.add_edges_from(other.edges)
                    return
                p_n = preds[0]  # The parent of n
                # For each node containing 3 or more child leaves in the other network, replace it by the
                # pendant subtree in this network
                self.remove_nodes_from(nodes2rem)
                r_other = other.get_root()
                # This is the i-th merged subnetwork
                subnet_prefixes = {n[0:3] for n in self.nodes if n[0:2] == '__'}
                i = len(subnet_prefixes)
                # Add the subnetwork below the parent of node n
                prefix = f"__{string.ascii_uppercase[i]}"
                # Relabel the internal nodes to prevent duplicate labels
                mapping = {}
                for n in other.nodes:
                    if n[0] == '_':  # If node n is an internal node
                        mapping[n] = prefix + n
                other = nx.relabel_nodes(other, mapping)
                self.add_edges_from(other.edges)
                self.add_edge(p_n, prefix + r_other)
                return
            if set(X_other).issubset(c):
                # Remove the children of node n that are not reticulations
                child_leaves = sorted(_n for _n in self.successors(n) if self.in_degree(_n) == 1 and _n[0] != '_')
                if not child_leaves:
                    continue
                if set(X_other).issubset(child_leaves):
                    # Merge the networks
                    self.remove_nodes_from(X_other)
                    r_other = other.get_root()
                    # This is the i-th merged subnetwork
                    subnet_prefixes = {n[0:3] for n in self.nodes if n[0:2] == '__'}
                    i = len(subnet_prefixes)
                    # Add the subnetwork below the parent of node n
                    prefix = f"__{string.ascii_uppercase[i]}"
                    # Relabel the internal nodes to prevent duplicate labels
                    mapping = {}
                    for _n in other.nodes:
                        if _n[0] == '_':  # If node n is an internal node
                            mapping[_n] = prefix + _n
                    other = nx.relabel_nodes(other, mapping)
                    self.add_edges_from(other.edges)
                    self.add_edge(n, prefix + r_other)
                    return
        print('ERROR: How should I merge these networks?')
        # self.draw(True, 'Final Net', all_labels=True)
        # other.draw(True, 'Merge It', all_labels=True)

    def _get_newick(self, _root=None, _hybrids=None):
        """
        Convert the phylogenetic network to a string in (e-)Newick format.

        :param _root: A root
        :param _hybrids: A list of hybrid nodes
        :return: A Newick-formatted string
        """
        if _hybrids is None:
            _hybrids = {}
        if _root is None:
            _root = self.get_root()
        sub_str = []
        for child in self[_root]:
            # Is this a hybrid node?
            if self.in_degree(child) > 1:
                if child in _hybrids.keys():
                    hybrid_label = _hybrids[child]
                    sub_str.append(hybrid_label)
                else:
                    hybrid_label = '#H' + str(len(_hybrids) + 1)
                    _hybrids[child] = hybrid_label
                    sub_str.append(self._get_newick(_root=child, _hybrids=_hybrids) + hybrid_label)
            elif len(self[child]) > 0:
                sub_str.append(self._get_newick(_root=child, _hybrids=_hybrids))
            else:
                sub_str.append(child)
        return "(" + ','.join(sub_str) + ")"

    def save_net(self, filename):
        """
        Save the given network in the Newick (or e-Newick) format.

        :param filename: A file name
        """
        nwk = self._get_newick()
        # Create the path if it does not exist
        path = os.path.dirname(filename)
        if len(path) and not os.path.exists(path):
            os.makedirs(path)

        with open(filename, 'w') as f:
            f.write(nwk + ';')
        print(f'SUCCESS: The network is saved as {filename}.')

    def draw(self, hold=False, title=None, axes=None, filename=None, all_labels=False, labels=None, edge_color=None,
             node_color=None):
        """
        Draw the given NetworkX digraph
        :param hold: Show the figure until the user closes it
        :param title: A title
        :param axes: Figure axes
        :param filename: A filename to store the figure
        :param edge_color:
        :param labels:
        :param all_labels:
        :param node_color:
        """
        if not LIVE_DRAW_ENABLED and not hold:
            print('WARNING: Drawing live is disabled.')
            return
        if not DRAW_ENABLED and stack()[1].function != 'draw_multi':  # TODO: remove stack
            print('WARNING: Drawing is disabled.')
            return
        if axes is None:
            plt.clf()  # HERE the program halts often :( WHY?
            plt.cla()

        def _decrease_leaf_labels():
            ys = set(y for (x, y) in pos.values())
            y_shift = (max(ys) - min(ys)) / (2 * len(ys))
            return {_n: (x, y - y_shift * (self.out_degree(_n) == 0)) for _n, (x, y) in pos.items()}

        # If self is a digraph, then nicely color its nodes
        if isinstance(self, nx.DiGraph):
            # Color the nodes: distinguish leaf nodes, hybrid nodes/reticulations and splits
            color_leaf = 'green'
            color_ret = 'red'
            color_splt = 'black' if not all_labels else 'lightGrey'
            color_lookup = {}
            node_size = []
            for node in self.nodes():
                if ',' in node:
                    raise RuntimeError(f'A node may not include a comma (I see node {node})')
                if self.out_degree(node) == 0:  # If this is a leaf
                    color_lookup[node] = color_leaf
                    node_size.append(100)
                elif self.in_degree(node) >= 2:  # If this is a hybrid node
                    color_lookup[node] = color_ret
                    node_size.append(60)
                else:
                    color_lookup[node] = color_splt
                    node_size.append(60)
                if node_color and node in node_color.keys():
                    color_lookup[node] = node_color[node]

            # Create a nice network layout
            pos = graphviz_layout(self, prog="dot")
            # Draw the network
            nx.draw_networkx(self, pos, True, with_labels=False, node_color=list(color_lookup.values()),
                             node_size=node_size, ax=axes, edge_color=edge_color)
            pos_nodes = _decrease_leaf_labels()
            # Only show the labels of the taxa/leaves
            if labels is None:
                labels = {}
                for n in self.nodes:
                    # If n is a leaf, then show the name
                    if self.out_degree(n) == 0 or all_labels:
                        labels[n] = n
            else:
                # Filter the labels for the leaves in this network
                nodes = self.get_leaves()
                labels = {n: l for n, l in labels.items() if n in nodes}
            nx.draw_networkx_labels(self, pos=pos_nodes, labels=labels, ax=axes)
        else:
            if labels is not None:
                # Filter the labels for the leaves in this network
                nodes = list(self.nodes)
                labels = {n: l for n, l in labels.items() if n in nodes}
            # Create a nice network layout
            nx.draw_circular(self, with_labels=True, edge_color=edge_color, node_color='lightGrey', ax=axes,
                             labels=labels)
            # If there are edge labels, show them
            e_labels = nx.get_edge_attributes(self, 'weight')
            if len(e_labels) > 0:
                nx.draw_networkx_edge_labels(self, nx.circular_layout(self), e_labels, ax=axes)

        # Set the title
        if title is not None:
            plt.gcf().suptitle(title)
        plt.subplots_adjust(left=0, bottom=0, right=1)
        if filename is not None:
            plt.pause(1E-9)
            plt.savefig(f'output/{filename}', bbox_inches='tight')
        if axes is None:
            if hold:
                plt.show()
            else:
                plt.pause(1E-9)
