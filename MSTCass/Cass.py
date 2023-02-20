# In the given set of clusters, find all ST-sets
import multiprocessing as mp
import queue
import string
import time
from collections import deque
from math import ceil
from multiprocessing.managers import BaseManager

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from networkx.algorithms.approximation import max_clique

from ClusterSet import ClusterSet
from IncGraph import IncGraph
from PhyloNet import PhyloNet

# Global program variables
DEBUG = False
FIND_ALL = False  # Find only one solution or all solutions
USE_MULTICORES = False
NUM_PROCESSES = 4  # mp.cpu_count()  # #Processes
MAX_IDLE_TIME = .1  # Maximum idle time [seconds]
ENABLE_DRAW = True

# Define the different CASS solving strategies
ORIGINAL = 0  # CASS: brute force approach
FILTER_EDGES = 1  # Filter the edges below which an ST-set could/should be hung. Needed
FILTER_X = 2  # Filter the possible ST-sets to choose from to remove in an iteration. Not used yet!
GUESS_EDGES = 4  # Guess an edge combination to hang the ST-set below. Not used yet!
GUESS_X = 8  # Guess which ST-set to remove next
STAY_BIN = 16  # Let the network stay binary
ALL_STSETS = 32  # Try removing any ST-set insted of only MST-sets
MIN_HANG_ONLY = 64  # Only use the smallest hang-edge combinations
FULL_CX = 128  # Let Cx comprise the minimal clusters only (False) or all clusters containing x (True)
STRATEGY = FILTER_EDGES | FILTER_X | GUESS_EDGES | GUESS_X | MIN_HANG_ONLY


class Combinations:
    """ Given a list of lists, generate al combinations of one value per list """

    def __init__(self, lst):
        if not lst:
            raise RuntimeError('Combinations should be initialized with a non-empty list')
        self.lst = lst
        self.sizes = [len(_l) for _l in lst.values()]
        self.counters = [0] * len(lst)
        self.i = 0
        self.stop = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop:
            raise StopIteration
        lst = []
        for i, cnt in enumerate(self.counters):
            lst.append(list(self.lst[list(self.lst)[i]])[cnt])
        # Update the counters
        c = 1  # Carry
        i = 0
        while c == 1 and i < len(self.counters):
            c = 0
            if self.counters[i] + 1 < self.sizes[i]:
                self.counters[i] += 1
            else:
                self.counters[i] = 0
                c = 1
            i += 1
        if c == 1:
            self.stop = True
        else:
            self.i += 1
        return lst


class TQueue(queue.LifoQueue):  # A LIFO queue results in a depth-first search :) & = faster than BFS
    def __init__(self):
        super(TQueue, self).__init__()
        self.id = 1

    def put(self, cn):
        cn.id = self.id
        self.id += 1
        super(TQueue, self).put(cn)


class PrioQueue(queue.PriorityQueue):  # A LIFO priority queue results in a depth-first search :) = faster than BFS
    def __init__(self):
        super(PrioQueue, self).__init__()
        self.id = 1

    def put(self, cn):
        cn.id = self.id
        self.id += 1
        # Give prio to continue building + give more prio to more complete networks
        super(PrioQueue, self).put((cn.net is None, cn))

    def get(self, block: bool = True, timeout: float = None):
        _, cn = super(PrioQueue, self).get(block, timeout)
        return cn


class RQueue(queue.Queue):
    """
    A result queue
    It checks whether a newly inserted network has a topology that was not found before
    """

    def __init__(self):
        super(RQueue, self).__init__()
        self.topos = deque()

    def put(self, net: PhyloNet):
        # Get the topology (ID)
        nodes = nx.topological_sort(net)
        topo = ''.join('_' if n[0] == '_' else n for n in nodes)
        # If a network having a new topology has been found, then save it
        if topo in self.topos:
            # This network has a topology that has been seen before
            return False
        self.topos.append(topo)
        super(RQueue, self).put(net)
        ret_num = net.get_ret_num()
        print(f"SUCCESS: We have found a new (sub)network of level {ret_num}!")
        return True


class ClusterSets:
    def __init__(self):
        self._lst = {}

    def save(self, cs: ClusterSet):
        c_str = ','.join(cs.X)
        self._lst[c_str] = cs.copy()

    def get_on(self, X):
        c_str = ','.join(X)
        if c_str not in self._lst.keys():
            raise RuntimeError(f'The leaf combination {c_str} has not been seen before')
        return self._lst[c_str].copy()


# Global variable
cs_seen = ClusterSets()


# Define an instance of a clusters-to-network instance
class CNInstance:
    def __init__(self, cs: ClusterSet, k: int, k1: int, X_removed=None, net: PhyloNet = None, comm_ancs=None,
                 parent_id=None):
        """
        Initialize
        :param cs: A set of clusters
        :param k: The desired level
        :param X_removed: A list of triples (b, l, X)
         If b is true, then the leaves X are implicitly removed when removing leaves l.
         If b is false, then the leaves X are collapsed to leaf l.
        """
        if not isinstance(cs, ClusterSet):
            raise TypeError(f"C should be a ClusterSet (instead of a {type(cs)}).")
        if X_removed is None:
            X_removed = deque()
        if comm_ancs is None:
            comm_ancs = {}
        self.cs = cs.copy()
        self.k = k
        self.k1 = k1
        self.X_removed = X_removed.copy()  # A dict of the removed leaves and values are the leaves that were removed afterwards
        self.net = net.copy() if net else None  # A PN
        self.comm_ancs = dict(comm_ancs)  # Common ancestor per cluster
        self.id = None
        self.parent_id = parent_id

    def __repr__(self):
        return f'CNInstance(id={self.id}, parent_id={self.parent_id}, cs={self.cs}, k={self.k}, k1={self.k1}, ' \
               f'X_removed={self.X_removed}, net={self.net}, LCA={self.comm_ancs})'

    def __lt__(self, other):
        if self.net:
            return self.k1 < other.k1
        return self.k1 > other.k1

    def collapse_mst_sets(self, print_instr=False):
        """
        Collapse all leaves that always occur together in the clusters
        :return: A collapsed set of clusters and a reduced set of taxa
        """
        cs_seen.save(self.cs)
        # Get the non-singleton MST-sets
        mst_sets = self.cs.get_mst_sets()
        if not mst_sets:
            return False  # Nothing collapsed
        if print_instr:
            print('INSTRUCTION: We have found some maximal ST-sets, which we are collapsing:')
        # Get the offset for the number of removed sets. Which letter may we use?
        letters = [_x for _x in self.cs.X if len(_x) == 1 and 'A' <= _x <= 'Z']
        if not letters:
            x_last = chr(ord('A') - 1)
        else:
            x_last = max(letters)
        for _, x, _ in self.X_removed:
            if len(x) == 1 and 'A' <= x <= 'Z' and x > x_last:
                x_last = x
        offset = ord(x_last) - ord('A') + 1
        for i, mst in enumerate(mst_sets):
            # Get a new meta-leaf name
            label = string.ascii_uppercase[i + offset]
            self.collapse_mst_set(mst, label, print_instr)
        return True

    def collapse_mst_set(self, mst, label, print_instr=False):
        # Get the clusters on this MST-set & save them
        cs_mst = self.cs.get_on_x(mst)
        cs_seen.save(cs_mst)
        if print_instr:
            print("INSTRUCTION: We collapse the leaves '" + "', '".join(mst) +
                  f"' to meta-taxon '{label}'")
        self.X_removed.append((False, label, mst))
        # Collapse the nodes in cs
        self.cs.collapse(mst, label)
        cs_seen.save(self.cs)

    def build_tree(self):
        """
        Build a NetworkX tree from this set of clusters
        :return:
        """
        self.net = self.cs.to_nxtree()
        if self.net:
            self._init_lca()
        return self.net

    def _init_lca(self):
        """
        For each cluster in C|X, determine the lowest common ancestor (LCA)
        :return:
        """
        for n in self.net.nodes():
            if n[0] != '_':
                continue
            # n is an internal node
            c = self.net.get_hardwired_cluster(n)
            # If this hard-wired cluster is in C, then LCA(C) = n
            if c in self.cs.C:
                _c_str = ','.join(c)
                lcas = set()
                lcas.add(n)
                self.comm_ancs[_c_str] = lcas

    def update_lca(self, x_added: str):
        lcas_new = {}
        for c in self.cs.C:
            c = c.copy()
            c_str = ','.join(c)
            # If x is not in this cluster, then LCA remains the same
            if x_added in c:
                # Find LCA(x_added, LCA(c\{x}))
                c.remove(x_added)
                if len(c) == 1:
                    lcas1 = c
                else:
                    c_str_old = ','.join(c)
                    if c_str_old not in self.comm_ancs.keys():
                        print(f'WARNING: LCA({c_str_old}) is unknown')
                        lcas1 = []
                    else:
                        lcas1 = self.comm_ancs[c_str_old]
                        assert isinstance(lcas1, set)
                        assert lcas1, print(self.comm_ancs)
                lcas = set()
                for lca1 in lcas1:
                    lcas2 = self.net.get_lcas_2(x_added, lca1)
                    lcas.update(lcas2)
                if not lcas:
                    print(self.comm_ancs)
                    print(lcas_new)
                    print(f'ERROR: The LCAs of {x_added} and {lcas1} could not be determined.')
                    # self.net.draw(True, f'Check LCAs of {x_added} and {lcas1}', all_labels=True)
                lcas_new[c_str] = lcas
            else:
                if c_str not in self.comm_ancs.keys():
                    print(f'ERROR: LCA({c_str}) does not exist')
                else:
                    lcas_new[c_str] = self.comm_ancs[c_str]
        # self.net.draw(True, 'Check LCAs', all_labels=True)
        self.comm_ancs = lcas_new

    def update_lca_after_dec(self, x: str, leaves):
        """
        Update the LCAs after decollapsing leaf x to the given leaves
        :param x:
        :param leaves:
        :return:
        """
        leaves = list(leaves)
        # For each cluster containing x, replace x by the leaves
        for c_s, lcas in list(self.comm_ancs.items()):
            c = c_s.split(',')
            if x in c:
                # Remove x by the leaves
                c.remove(x)
                c.extend(leaves)
                del self.comm_ancs[c_s]
                self.comm_ancs[','.join(sorted(c))] = lcas

        # Update the clusters on the leaves
        for _x in leaves:
            p = _x
            while True:
                preds = list(self.net.predecessors(p))
                if not preds:
                    print(f'ERROR: Check the predecessors of node {p}')
                    self.net.draw(True, f'Check Preds({p})')
                p = preds[0]
                descs = sorted(x for x in nx.descendants(self.net, p) if x[0] != '_')
                lcas = set()
                lcas.add(p)
                self.comm_ancs[','.join(descs)] = lcas
                if descs == leaves:
                    break
                if set(leaves).issubset(descs):
                    self.comm_ancs[','.join(leaves)] = lcas
                    break

    def get_hang_edges(self, x: str, print_instr=False):
        return self.net.get_hang_edges(self.cs, x, self.comm_ancs,
                                       filter_edges=(STRATEGY & FILTER_EDGES == FILTER_EDGES),
                                       full_cx=(STRATEGY & FULL_CX == FULL_CX), print_instr=print_instr)

    def add_leaves(self, X_add):
        """
        Let X_t be the target set of leaves in the network after adding the ST-set. Add the leaves that are in X_t,
        though not in the network, as leaves below a new root
        :param X_add: the leaves to add before adding x (BitString)
        :return: The network in which the necessary leaves are added
        """
        if not X_add:
            return

        # Add these nodes below a new root, the parent of the current root
        root, old_root = self.net.add_new_root()
        self.net.add_edges_from([(old_root, n) for n in X_add])

        # Update C, X
        X_new = list(self.cs.X)
        X_new.extend(X_add)
        X_new.sort()
        self.cs = cs_seen.get_on(X_new)

    def remove_leaf(self, x: str, print_instr=False):
        X_unused = self.cs.remove_leaf(x)
        self.X_removed.append((True, x, X_unused))
        if print_instr:
            print(f'INSTRUCTION: We remove leaf {x}.')
            if X_unused:
                print(f'INSTRUCTION: Thereby, leaves {X_unused} are removed implicitly.')
        if X_unused:
            cs = self.cs.copy()
            cs.X.extend(X_unused)
            cs.X.sort()
            cs_seen.save(cs)

    def copy(self):
        """ Deep copy """
        return CNInstance(self.cs, self.k, self.k1, self.X_removed, self.net, self.comm_ancs, self.parent_id)


def draw_multi(nets: list[PhyloNet], title: str = None, hold=True, file=None, full_screen=False, all_labels=False,
               labels=None):
    """
    Draw all networks that are given, in a single figure
    :param nets: A list of networks
    :param title: A title
    :param hold: Hold the drawing T/F
    :param file: A filename to store the figure as image
    :param full_screen: Make the figure screen filling T/F
    :param all_labels: Show all labels T/F
    :param labels: A list of node labels, which can be used to replace the labels defined in the network
    """
    if not ENABLE_DRAW:
        return
    n_figs = len(nets)  # #figures
    if n_figs == 0:  # If there is nothing to show, then return
        print("WARNING: There is nothing to show.")
        return
    N_MAX = 12
    if n_figs > N_MAX:
        nets = nets[:N_MAX]  # Take the first 9 networks
        title += f' ({N_MAX}/{n_figs})'
        n_figs = N_MAX

    plt.clf()  # NOTE: In the debug mode, this might cause a return in the main process. I don't know why...
    n_cols = int(np.ceil(np.sqrt(n_figs)))
    n_rows = n_cols - 1 if n_figs <= n_cols * (n_cols - 1) else n_cols
    fig, ax = plt.subplots(n_rows, n_cols, num=1)
    for i, net in enumerate(nets):
        if n_cols == 1:
            axes = ax
        else:
            ix = np.unravel_index(i, ax.shape)
            axes = ax[ix]
        if isinstance(net, nx.DiGraph):
            # If this is a PN
            net.draw(hold=hold, axes=axes, all_labels=all_labels, labels=labels)
        else:
            # If this is an IG
            net.draw(axes=axes)

    # Set the title, etc.
    if title is not None:
        fig.suptitle(title)
    if full_screen:
        mgr = plt.get_current_fig_manager()
        if not getattr(mgr, "flag_is_max", None):
            mgr.full_screen_toggle()
            mgr.flag_is_max = True
    if file is not None:
        plt.pause(1E-9)
        plt.savefig(f'output/{file}', bbox_inches='tight')  # Save the figure
    if hold:
        plt.show()
    else:
        plt.pause(1E-9)


def get_ret_bounds(cc: ClusterSet, IG=None, ig_is_min=False, report_bounds=False):  # TODO: finish
    """
    Get bounds on the reticulation number
    :return: A lower and upper bound, respectively
    """
    cs = cc.copy()
    # Create the IG on the given set of clusters
    if IG is None:
        IG = IncGraph(cs)
        if not IG.edges:
            if report_bounds:
                return {','.join(cs.X): {
                    'min_ints': 0,
                    'min_rests': 0,
                    'gl_min_ints': 0,
                    'gl_min_rests': 0,
                    'n_min_rests': 0,
                    # 'frac_min_diffs': 0,
                    'clique': 0,
                    'max': 0,
                    # 'max2': 0
                }}
            return 0, 0
        # Remove singletons
        IG.remove_nodes_from(list(nx.isolates(IG)))
    else:
        # Note: Collapse MST-sets before
        IG = IG.subgraph(cs.C)

    bounds_report = {}
    # Get the connected components in IG
    ccs = nx.connected_components(IG)
    # Analyse each conn. comp.
    r_min = r_max = 0
    for cc in ccs:
        C = [c.split(',') for c in cc]
        cs = ClusterSet(C)
        bounds = _get_cc_ret_bounds(cs, IG.subgraph(cc), ig_is_min, report_bounds=report_bounds)
        if report_bounds:
            # bounds['n_min_rests'] = bounds['n_min_rests'] / bounds['']
            bounds_report[','.join(cs.X)] = bounds
        else:
            r_min += bounds[0]
            r_max += bounds[1]
    if report_bounds:
        return bounds_report
    return r_min, r_max


def _get_cc_ret_bounds(cs: ClusterSet, IG: IncGraph, ig_is_min=False, report_bounds=False):
    """
    Get the bounds on r(G) for the given connected component G in IG(C)
    :param cs:
    :param IG:
    :param ig_is_min:     IG is of the minimal clusters T/F
    :param report_bounds:
    :return:
    """
    # Init
    r_max = min(len(cs.X), len(cs.C)) - 1
    bounds = {
        'min_ints': 0,
        'min_rests': 0,
        'gl_min_ints': 0,
        'gl_min_rests': 0,
        'n_min_rests': 0,
        'frac_min_diffs': 0,
        'clique': 0,
        'max': r_max,
        # 'max2': r_max,
    }
    if len(cs.C) == 1:
        if report_bounds:
            return bounds
        return 0, 0
    r_min = 1
    nds = list(IG.nodes)

    # Create a digraph G, with the nodes in IG, indicating clusters that contain others
    # G = IG.get_contain_graph()
    # G = G.to_undirected()
    all_ints = set()  # All intersections
    all_diffs = set()  # All rests
    for i, c1 in enumerate(nds):
        cl1 = set(c1.split(','))

        # Minimal differences degree bound. This is only applicable to IG_min
        diffs = []
        for e in IG.edges(c1):  # For each neighbouring edge e of node c1
            _c1 = e[0].split(',')
            _c2 = e[1].split(',')
            diff = list(set(_c2).difference(_c1))
            diff.sort()
            if diff not in diffs:
                diffs.append(diff)
                # print(f'INSTRUCTION: {_c2} \\ {_c1} = {diff}')
            all_diffs.add(','.join(diff))
        min_diffs = ClusterSet(diffs, cs.X).minimal()
        r = len(min_diffs.C) - 1
        if r > r_min:
            r_min = r
            # print(f'INFO: Min. diff. degree of node {c1} = {r + 1} => r >= {r}')
        if r > bounds['min_rests']:
            bounds['min_rests'] = r

        # Triangle => min. 2 reticulations needed
        nb1 = [n for n in IG.neighbors(c1)]
        for c2 in nb1:
            # If this edge is in a triangle, #reticulations >= 2
            comm_nbs = set(nx.neighbors(IG, c1)).intersection(nx.neighbors(IG, c2))
            if comm_nbs:
                if r_min < 2:
                    r_min = 2
                    # print(f'INFO: Triangle in IG => r >= 2')
                break

        # Special K(2,n) bound
        # Minimal incompatibility degree bound
        intersections1 = {n: list(cl1.intersection(n.split(','))) for n in nb1}
        # Get the minimal clusters among the intersections with neighbours
        int_unique = set(','.join(n) for n in intersections1.values())
        all_ints.update(int_unique)
        cs_int_min = ClusterSet(list(c.split(',') for c in int_unique)).minimal()
        # Count the number of neighbours having the same minimal intersection
        int_cnt = {}
        for n, inters in intersections1.items():
            if sorted(inters) in cs_int_min.C:
                key = ','.join(inters)
                if key not in int_cnt.keys():
                    int_cnt[key] = 1
                else:
                    int_cnt[key] += 1
        X_ints_len = len(cs_int_min.get_used_taxa())
        num_ints = len(cs_int_min.C)  # #minimal intersections
        if len(cl1) > X_ints_len:
            if r_min < num_ints:
                # print(f'INFO: A) By a minimal neighbour intersection set of {c1} of size {num_ints}, '
                #       f'r >= {num_ints}.')
                r_min = num_ints
            if bounds['min_ints'] < r_min:
                bounds['min_ints'] = r_min
        else:
            # print(f'INFO: B| X_ints_len = {X_ints_len}, num_ints = {num_ints}, r_min = {r_min}')
            if r_min < num_ints - 1:
                # print(f'INFO: B) By a minimal neighbour intersection set of {c1} of size {num_ints}, '
                #       f'r >= {num_ints - 1}.')
                r_min = num_ints - 1
            if bounds['min_ints'] < r_min:
                bounds['min_ints'] = r_min
    # Minimal intersections globally
    intersections = [c.split(',') for c in all_ints]
    cs_min_ints = ClusterSet(intersections).minimal()
    r = int(ceil(len(cs_min_ints.C) / 3))
    bounds['gl_min_ints'] = r
    if r > r_min:
        r_min = r
        print(f"INFO: Globally, there are {len(cs_min_ints.C)} different minimal intersections => r >= {r}")

    # Minimal differences globally
    diffs = [c.split(',') for c in all_diffs]
    cs_min_diffs = ClusterSet(diffs).minimal()
    amt = len(cs_min_diffs.C)
    r = int(ceil(amt / 5))
    bounds['gl_min_rests'] = r
    bounds['n_min_rests'] = amt
    if r > r_min:
        r_min = r
        # print(f"INFO: Globally, there are {amt} different minimal differences => r >= {r}")

    # Max clique bound
    c = max_clique(IG)
    r = len(c) - 1
    if r > r_min:
        r_min = len(c) - 1
        # print(f"INFO: I've found a max. clique of size {len(c)} => r >= {r_min}")
    bounds['clique'] = r

    # Calculate an upper bound on the number of reticulations
    # For each leaf, determine the number of minimal incompatible clusters. This is a measure for an UB.
    L = []
    for x in list(cs.X):
        cs_in = cs.get_cs_in(x)
        # cs_out = cs.get_cs_out(x, cs_in)
        l = len(cs_in.C)  # r(G) <= (|C_x| - 1) + 1 (to e.g. the root)
        L.append(l)
    r = sum(L[2:])
    if r < r_max:
        # print(f"INFO: By Conjecture 7.2, r <= {r}")
        r_max = r
    if r_min > r_max:
        raise RuntimeError(f'r_min = {r_min} > r_max = {r_max}')
    bounds['max'] = r_max

    if report_bounds:
        return bounds
    return r_min, r_max


def get_level_and_ret_bounds(cs: ClusterSet):
    """
    Get bounds on the level and number of reticulations of a network that represents the clusters
    :return: A lower and upper bound, respectively
    """
    # Collapse the maximal ST-sets
    cs = cs.copy()
    cs.collapse_mst_sets()

    IG = IncGraph(cs)
    if not IG.edges:
        return 0, 0, 0, 0
    IG_min = IncGraph(cs.minimal())

    # Remove singletons
    IG.remove_nodes_from(list(nx.isolates(IG)))
    IG_min.remove_nodes_from(list(nx.isolates(IG_min)))

    # Init
    r_min = 0  # Min. ret. num.
    r_max = 0
    l_min = 0  # Min. level
    l_max = 0  # Max. level
    # For each connected component
    ccs = IG.get_nontrivial_cc()
    for cc in ccs:
        rl, ru = get_ret_bounds(cc, IG)
        # print(f'INFO: {rl} <= r <= {ru} reticulations are needed for connected component {cc}.')
        r_min += rl
        r_max += ru
        l_min = max(l_min, rl)
        l_max = max(l_max, ru)

    # Repeat the above using IG_min
    r_min2 = 0  # Min. ret. num.
    r_max2 = 0
    l_min2 = 0  # Min. level
    l_max2 = 0  # Max. level
    # For each connected component
    ccs = IG_min.get_nontrivial_cc()
    for cc in ccs:
        rl, ru = get_ret_bounds(cc, IG_min)
        # print(f'INFO: In IG_min, {rl} <= r <= {ru} reticulations are needed for conn. comp. {cc}.')
        r_min2 += rl
        r_max2 += ru
        l_min2 = max(l_min2, rl)
        l_max2 = min(l_max2, ru)
    l_min = max(l_min, l_min2)
    l_max = min(l_max, l_max2)
    r_min = max(r_min, r_min2)
    r_max = min(r_max, r_max2)
    return l_min, l_max, r_min, r_max


def optcass(cs: ClusterSet, strategy: int = ORIGINAL, timeout=None):
    """
    Find a phylogenetic network for the given set of clusters C with minimum level
    :param cs: A set of clusters
    :param strategy: The search strategy, e.g. FILTER_ST | FILTER_EDGES | GUESS_EDGES | GUESS_ST
    :param timeout: A timeout in seconds
    :return: A max k-level network
    """
    # Set the strategy
    global STRATEGY
    STRATEGY = strategy

    print('INFO: OptCass starts searching a network for you...')
    strats = []
    if STRATEGY & FILTER_EDGES == FILTER_EDGES:
        strats.append('filter the edges')
    if STRATEGY & FILTER_X == FILTER_X:
        strats.append('filter the leaves to remove/add')
    if STRATEGY & GUESS_EDGES == GUESS_EDGES:
        strats.append('guess good edges')
    if STRATEGY & GUESS_X == GUESS_X:
        strats.append('guess good leaves to remove/add')
    if STRATEGY & ALL_STSETS == ALL_STSETS:
        strats.append('try removing/adding all ST-sets (i.e. not only MST-sets)')
    # print('INFO: The following strategies are used: ' + ', '.join(strats))
    if timeout:
        print(f'INFO: The program will run at most {timeout / 60} minutes.')

    cn = CNInstance(cs, 0, 0)

    t0 = time.time()  # Record the elapsed time
    q_res = mp.Queue()  # Record the results
    p = mp.Process(target=_optcass, name='OptCass', args=(q_res, cn, timeout))
    p.start()
    p.join(timeout)
    p.terminate()
    t_tot = time.time() - t0
    print('INFO: OPT_CASS needed {:.3f} s'.format(t_tot))
    # Get the results
    nets = []
    k = 0
    while not q_res.empty():
        cn = q_res.get()
        if cn.net:
            nets.append(cn.net)
        k = max(k, cn.k)
    return nets, t_tot, k


def _optcass(q_cn_results, cn: CNInstance, timeout=None) -> None:
    # Can a tree be formed? If yes, return it!
    if cn.build_tree():
        cn.k = 0
        q_cn_results.put(cn)
        return

    # For each connected component in IG(C), build a subnetwork that represents teh clusters therein.
    # Create the incompatibility graph
    IG = IncGraph(cn.cs)

    # The isolates determine the global network structure
    isos = list(nx.isolates(IG))
    C_iso = [x.split(',') for x in isos]
    cs_iso = ClusterSet(C_iso, cn.cs.X)
    final_net = cs_iso.to_nxtree()
    # final_net.draw(True, 'The Overall Structure')

    # Determine the minimum level l_min
    l_min, l_max, _, _ = get_level_and_ret_bounds(cn.cs)
    print(f"INFO: The level will lie between {l_min} <= k <= {l_max}.")

    # Init
    level = 0
    r_tot = 0
    # Get the non-singleton connected components
    ccs = IG.get_nontrivial_cc()
    ccs = sorted(ccs, key=lambda cs: len(cs.X), reverse=True)  # Handle the largest connected component firstly
    # Add each isolate to the smallest CC containing the leaves of the iso. plus at least one more
    for iso in isos:
        iso = set(iso.split(','))
        for i, cc in reversed(list(enumerate(ccs))):
            # If iso is a strict subset of X(CC)
            if iso.issubset(cc.X) and len(iso) < len(cc.X):
                c = sorted(iso)
                if c not in cc.C:
                    cc.C.append(c)
                break

    for i, cc in enumerate(ccs):  # For each connected component
        # Search a network with the lowest number of reticulations
        found = False
        rl, ru = get_ret_bounds(cc, IG)
        k = rl  # The level
        while not found:
            found = True
            # print(f"INFO: I start searching a level-{k} network for the connected component on X={cc.X}")
            # Apply the algorithm
            t0 = time.time()  # Record the starting time to determine the runtime
            _cn = CNInstance(cc, k, k)
            q_cn_results.put(_cn)  # Save this instance to be able to derive that the level is >= k, on a timeout
            # Run the actual algorithm
            nets = optcass_k(_cn, timeout)
            # Record the runtime
            t = time.time() - t0
            print('INFO: OPT_CASS(k={:}) needed {:.3f} s'.format(k, t))
            for n in nets:
                n.remove_dummy_root()
            if len(nets) >= 1:
                r = nets[0].get_ret_num()
                r_tot += r
                print(f'k={max(level, k)},r={r_tot}')  # System output: level k and ret.num. so far
                final_net.merge(nets[0])
            elif k == l_max:
                raise RuntimeError(f'I could not find any network for all {l_min} <= k <= {l_max}. '
                                   'There might be a bug in your code.')
            else:
                found = False
                k += 1  # Search for a higher level network
                if not q_cn_results.empty():
                    q_cn_results.get(False)
        level = max(level, k)
    # Decollapse the remained super/meta nodes
    cn.net = final_net
    cn.k = level
    while cn.X_removed:
        removed, x, X_add = cn.X_removed.pop()
        assert not removed  # The leaves in X_add had been collapsed
        _decollapse(cn, x, X_add)
    q_cn_results.put(cn)


def _add_leaf(q, q_res, cn: CNInstance, x: str, X_add, print_instr=False):
    if X_add:
        if print_instr:
            print(f"INSTRUCTION: Since the leaves {X_add} were removed at removing leaf '{x}', we add them now.")
        cn.add_leaves(X_add)

    X_new = list(cn.cs.X)
    assert x not in X_new, cn.net.draw(True, f'CHECK: ADD {x}?', all_labels=True)
    # Update the set of clusters (i.e. clusters C and taxa X)
    X_new.append(x)
    X_new.sort()
    cn.cs = cs_seen.get_on(X_new)
    E = cn.get_hang_edges(x, print_instr=DEBUG)
    _MIN_LEN_ONLY = STRATEGY & MIN_HANG_ONLY == MIN_HANG_ONLY
    if _MIN_LEN_ONLY:
        # Analyse the edge combinations. Only use those that result in the lowest number of reticulations in the resulting N
        min_len = len(cn.cs.X) - 1
        for edges in Combinations(E):
            r = len(set(edges))
            if r < min_len:
                min_len = r
    ret_num0 = cn.net.get_ret_num()
    found = False
    for edges in Combinations(E):
        edges = set(edges)
        if _MIN_LEN_ONLY and len(edges) > min_len:
            continue
        # Note: #edges - 1 = #reticulations that will be added.
        # If the number of reticulations in the new network would exceed the maximum allowed number, then return
        new_ret_num = cn.net.get_ret_num() + len(edges) - 1
        if new_ret_num > cn.k:
            continue
        net = cn.net.copy()
        nets = net.add_leaf_below_ret(x, edges, cn.cs, STRATEGY & STAY_BIN == STAY_BIN)
        if not nets:
            if print_instr:
                print(
                    f'INSTRUCTION: I could not find any level-{cn.k} network after inserting leaf {x} below the edges '
                    f'{E}.')
                continue
        for net in nets:
            ret_num = net.get_ret_num()
            num_rets_added = ret_num - ret_num0
            k1_new = cn.k1 + num_rets_added
            # If the network contains dummy leaves, then save the network to handle later. No LCAs update
            if net.get_dummy_leaves():
                cn1 = CNInstance(cn.cs, cn.k, k1_new, cn.X_removed, net, cn.comm_ancs, cn.id)
                q.put(cn1)
                continue
            # Does the new network represent all clusters?
            # cn.update_lca(x)
            if not net.represents(cn.cs, x, cn.comm_ancs):  # This may be the case for eluFig9
                continue

            # cn1 = CNInstance(cn.cs, cn.k, k1_new, cn.X_removed, net, cn.comm_ancs, cn.id)
            # if not cn1.net_represents_cl(x, edges):  # TODO: let it work; it would be faster than net.represents!
            #     continue
            found = True
            # TODO
            # comm_ancs_new = _update_lca(cn.net, net, cs_new, st[0], cn.comm_ancs)
            # cn1 = CNInstance(cs_new, cn.k, k1_new, cn.X_removed, net, comm_ancs_new, cn.id)
            cn1 = CNInstance(cn.cs, cn.k, k1_new, cn.X_removed, net, cn.comm_ancs, cn.id)
            # If the next step is a decollapse of leaves, then do it immediately
            if cn1.X_removed:
                toAdd, y, leaves = cn1.X_removed.__getitem__(len(cn1.X_removed) - 1)
                if not toAdd:
                    cn1.X_removed.pop()
                    cn1.update_lca(x)  # TODO: use LCAs
                    _decollapse(cn1, y, leaves, print_instr=print_instr)
            # If all original leaves are present in this network, then we have found a final network!
            if not cn1.X_removed:
                cn1.net.clean(True)
                q_res.put(cn1.net)
            else:
                # A new network has been found, containing leaf x.
                cn1.update_lca(x)  # TODO: use LCAs
                q.put(cn1)
    if not found and print_instr and new_ret_num > cn.k:
        # new_ret_num = cn.net.get_ret_num() + len(edges) - 1
        print(f'INSTRUCTION: There would be too many reticulations ({new_ret_num}, while k = {cn.k}) '
              'in the resulting network.')


def _decollapse(cn, x, leaves, print_instr=False):
    if print_instr:
        print(f"INSTRUCTION: We decollapse meta-leaf '{x}' to {leaves}")
    cs_mst = cs_seen.get_on(leaves)
    subtree = cs_mst.to_nxtree()
    # Replace the meta-node by the pendant subtree
    preds = list(cn.net.predecessors(x))
    if len(preds) == 0:
        # ERROR
        print(f'ERROR: Check the predecessors of node {x}')
        cn.net.draw(True, f'Predecessors {x}?', all_labels=True)
    p_x = preds[0]  # parent of x
    cn.net.remove_node(x)
    # Merge the network & subtree
    r = subtree.get_root()
    # This is the i-th pendant subtree that is decollapsed
    dec_leaf_prefixes = {n[0:3] for n in cn.net.nodes if n[0:2] == '__'}
    i = len(dec_leaf_prefixes)
    # Add the subnetwork below the parent of node n
    prefix = f"__{string.ascii_uppercase[i]}"
    cn.net = nx.union(cn.net, subtree, rename=(None, prefix))
    # Relabel the leaves/taxa to their original names
    mapping = {}
    for leaf in subtree.get_leaves():
        mapping[prefix + leaf] = leaf
    cn.net = nx.relabel_nodes(cn.net, mapping)
    childs = list(cn.net.successors(p_x))
    # If node p_x has childs and if they are all reticulations, then merge p_x and r
    non_rets = [n for n in cn.net.successors(p_x) if cn.net.in_degree(n) == 1]
    if childs and len(non_rets) == 0:
        # Merge p_x and r
        childs = cn.net.successors(prefix + r)
        cn.net.remove_node(prefix + r)
        cn.net.add_edges_from((p_x, n) for n in childs)
        # Node x has been decollapsed now, and merged with its former parent
    else:
        # Connect the parent p_x of meta-leaf x to the root of the tree
        cn.net.add_edge(p_x, prefix + r)

    # Update C and X
    cn.cs.X.remove(x)
    cn.cs.X.extend(leaves)
    cn.cs = cs_seen.get_on(sorted(cn.cs.X))
    # Update the LCAs
    cn.update_lca_after_dec(x, leaves)  # TODO: use LCAs


def process_cn(q, q_res, cn: CNInstance, print_instr=False):
    # Is this instance in phase 1, i.e. should a leaf be removed?
    if cn.net is None:
        # Try to build a tree from the remaining clusters. (tree = None if no such tree exists.)
        cn.build_tree()

        # If k1 == 0, then return the tree, if it exists
        if cn.k1 == 0 and cn.net is None:
            return
        if cn.k1 == 0 or cn.net:
            # Add a new dummy root
            cn.net.add_new_root()
    if cn.net is None:
        cs_seen.save(cn.cs)
        # Remove an MST-set
        if cn.X_removed:
            msts_only = not (STRATEGY & ALL_STSETS == ALL_STSETS)
        else:
            # In the first iteration, always remove an MST-set (otherwise the network will not be optimal)
            msts_only = True
        for m in cn.cs.get_st_sets2rem(msts_only, STRATEGY & GUESS_X == GUESS_X):
            # Remove an (M)ST-set
            cn1 = cn.copy()
            cn1.k1 -= 1
            if len(m) == 1:
                cn1.remove_leaf(m[0])
            else:
                cn1.collapse_mst_set(m, 'MST', DEBUG)
                cn1.remove_leaf('MST', DEBUG)
            cn1.parent_id = cn.id
            q.put(cn1)
    else:  # Is this instance in phase 2, i.e. should a leaf be added?
        # Continue building a network
        # If there are dummy leaves in the network, then add a reticulation
        dummies = cn.net.get_dummy_leaves()
        if dummies:
            # leaves = [n for n in cn.cs.X if n not in cn.net.get_leaves()]
            leaves = set(cn.cs.X).difference(cn.net.get_leaves())
            if len(leaves) >= 2:
                print(leaves)
                print('ERROR: Check this!')
            x = leaves.pop()
            # Choose two dummy leaves and add a reticulation
            nets = cn.net.add_ret_below_2(dummies, x)
            for net in nets:
                cn1 = cn.copy()
                cn1.k1 += 1
                cn1.parent_id = cn.id
                cn1.net = net
                # If x has been added, update the LCAs
                if x in net.get_leaves():
                    # Check if we have found a solution
                    if not cn1.X_removed:
                        q_res.put(cn1.net)
                        continue
                    cn1.update_lca(x)  # TODO: use LCAs
                q.put(cn1)
            return
        added = False
        if not cn.X_removed:
            represents_clusters = cn.net.represents(cn.cs, None, cn.comm_ancs, True)
            print(represents_clusters)
            print('ERROR: cn.X_removed is empty')
            # cn.net.draw(True, 'Nothing to remove. OK?', all_labels=True)
            # raise RuntimeError('cn.X_removed is empty.')
        else:
            while not added and cn.X_removed:  # While decollapsed
                added, x, X_add = cn.X_removed.pop()
                if added:
                    _add_leaf(q, q_res, cn, x, X_add, print_instr)
                else:
                    _decollapse(cn, x, X_add, print_instr)
        # It may be that only decollapse steps are done. Then a solution is found!
        if not added and not cn.X_removed:
            # A solution has been found
            cn.net.remove_dummy_root()
            q_res.put(cn.net)


def queue_handler(proc_busy, q, q_res, _cs_seen):
    global cs_seen
    cs_seen = _cs_seen
    worker_name = mp.current_process().name

    while True:
        # If a solution has been found and no more are needed, then stop.
        if not FIND_ALL and not q_res.empty():
            break

        cn = None
        try:
            cn = q.get(False)
        except queue.Empty:
            # This process is idle now
            if worker_name in proc_busy:
                proc_busy.remove(worker_name)
            # If all processes are done, then stop the process
            if not proc_busy:
                break
            # There is some other process running. Wait.
            time.sleep(MAX_IDLE_TIME)

        if cn:
            # This process starts handling a problem instance
            if not worker_name in proc_busy:
                proc_busy.append(worker_name)
            process_cn(q, q_res, cn)

    exit(1)


def optcass_k(cn: CNInstance, timeout=None):
    """
    Run the OPT_CASS (Optimal Cass) algorithm
    :param cn: A CNInstance
    :param timeout: A timeout [s]
    :return: A list of CNInstances
    """

    cn.collapse_mst_sets()
    if USE_MULTICORES:
        q_res = mp.Queue()  # A queue to store the results
        proc_busy = mp.Manager().list()
        BaseManager.register('TQueue', TQueue)
        # BaseManager.register('TQueue', PrioQueue)  # Does not work fast
        BaseManager.register('ClusterSets', ClusterSets)
        with BaseManager() as mgr:
            q = mgr.TQueue()  # The queue for the instances to process
            _cs_seen = mgr.ClusterSets()
            q.put(cn)

            # Create the processes
            procs = []
            for i in range(NUM_PROCESSES):
                procs.append(mp.Process(target=queue_handler, args=(proc_busy, q, q_res, _cs_seen)))

            # Run the processes
            for proc in procs:
                proc.start()

            # Wait (block) until all processes are done
            for proc in procs:
                proc.join(timeout)

            for p in procs:
                p.terminate()
    else:
        q = TQueue()  # The queue for the instances to process
        # q = queue.Queue()  # This queue works too SLOW
        # q = PrioQueue()  # The queue for the instances to process
        q_res = RQueue()  # A queue to store the results
        process_cn(q, q_res, cn)
        while not q.empty() and (FIND_ALL or q_res.empty()):
            process_cn(q, q_res, q.get())

    # Read the results
    nets = []
    while not q_res.empty():
        nets.append(q_res.get())
    print(f'INFO: I finished with {len(nets)} result(s).')
    return nets
