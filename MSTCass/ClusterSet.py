from collections import Counter

import networkx as nx

from PhyloNet import PhyloNet


class Cluster(list):
    def __new__(cls, c=None):
        # if not isinstance(c, (list, set)):
        #     raise RuntimeError(f"A Cluster should be initialized with a list or set.\nc = {c}")
        return list.__new__(cls, c)

    def issubset(self, other):
        return set(self).issubset(other)

    def intersection(self, other):
        """ Get the intersection of two clusters, assumed that they are sorted """
        return list(set(self).intersection(other))

    def iscompatible(self, other):
        """ Check if this cluster is compatible with another one

        At least one of the following conditions should hold:
         1. The one is contained in the other, or
         2. The clusters are independent
        """
        return self.issubset(other) or set(other).issubset(self) or not self.intersection(other)


class ClusterSet:
    """
    Helper class for a set of clusters, where each cluster is a list of taxa
    """

    def __init__(self, C, X=None):
        if X is None:
            if isinstance(C, ClusterSet):
                X = C.X
                C = C.C
            elif isinstance(C, list):
                X = set()
                for c in C:
                    X.update(c)
            else:
                raise RuntimeError('What to do with C?')
        elif not isinstance(X, (list, set)):
            print(X)
            raise RuntimeError('X should be a list or set of taxa')
        else:
            C = list(C)
        # Create a list of Clusters
        self.C = []
        for c in C:
            self.C.append(Cluster(c))
        self.X = list(X)

    def __repr__(self):
        return f"ClusterSet({len(self.C)} clusters {self.C} on X={self.X})"

    def copy(self):
        return ClusterSet(self.C, self.X)

    def get_incompatibilities(self):
        """
        Get a set of incompatibilities.

        Example
        -------
        >>> cs = ClusterSet([['A', 'B'], ['B', 'C']])
        >>> incs = cs.get_incompatibilities()
        >>> incs
        [BitString('010')]

        :return: A sorted set of incompatibilities
        """
        incs = []
        for i, cl1 in enumerate(self.C):
            for j in range(i + 1, len(self.C)):
                cl2 = self.C[j]
                if not cl1.iscompatible(cl2):
                    inc = list(set(cl1).intersection(cl2))
                    inc.sort()
                    if inc not in incs:
                        incs.append(inc)
        return sorted(incs, key=lambda inc: len(inc))

    def get_cinc_freqs(self):  # TODO: Needed? Then clarify
        """
        Get the cluster in cluster frequencies.
        Be careful: this is an expensive function
        :return: The cluster in cluster frequencies
        """
        cinc = {}
        for c1 in self.C:
            cnt = 0
            for c2 in self.C:
                if c1 != c2 and set(c2).issuperset(c1):
                    cnt += 1
            cinc[','.join(c1)] = cnt
        return cinc

    def is_separating(self):
        """ Is this set of clusters separating, i.e. is there no pair of conflicting clusters? """
        for n1 in self.C:
            c1 = Cluster(n1)
            for n2 in self.C:
                if not c1.iscompatible(n2):
                    return False
        return True

    def get_used_taxa(self):
        """
        Get a BitString that represents the taxa that are in at least one cluster
        :return: A BitString
        """
        X = set()
        for c in self.C:
            X.update(c)
        X = list(X)
        X.sort()
        return X

    def is_subset(self, c: list):
        """
        Check if cluster c is a strict subset of another one in C
        :param c: A cluster
        :return: c is a strict subset T/F
        """
        is_sub = False
        for c2 in self.C:
            is_sub |= set(c2).issuperset(c) and c != c2
        return is_sub

    def is_supset(self, c):
        """
        Check if cluster c is a strict superset of another one in C
        :param c: A cluster
        :return: c is a strict superset T/F
        """
        is_super = False
        for c2 in self.C:
            is_super |= c2.issubset(c) and c != c2
        return is_super

    def minimal(self):
        """
        Get the minimal clusters
        :return: A set of minimal clusters
        """
        return ClusterSet(filter(lambda c: not self.is_supset(c), self.C), self.X)

    def maximal(self):
        """
        Get the maximal clusters
        :return: A set of minimal clusters
        """
        return ClusterSet(filter(lambda c: not self.is_subset(c), self.C), self.X)

    def get_cs_in(self, x: str, full_cx=False):
        # Get the clusters containing leaf x, making the set of clusters cs_in
        cs_in = ClusterSet(filter(lambda c: x in c, self.C), self.X)
        if full_cx:
            return cs_in
        # Let the minimal clusters remain; these determine where leaf x should be hung
        return cs_in.minimal()

    def get_cs_out(self, x: str, cs_in: 'ClusterSet'):
        """
        Given that I am cs_in, check if the cluster X_in\{x} exists as (not necessarily strict) subset of any cluster
        in C. Make a list (set) of all such edges.
        :param x: The leaf being added
        :param cs_in: cs_in
        :return: A list of edges
        """
        c_out = set(cs_in.get_used_taxa())
        c_out.remove(x)
        cs_out = ClusterSet([], self.X)
        if len(c_out) == 1:  # If a singleton is left, it will always be represented
            return cs_out
        for c in self.C:
            # If cluster c contains c_out and not leaf x
            if c_out.issubset(c) and x not in c:
                cs_out.C.append(c)
                # Cluster c in C contains c_out. Hence, the ST-set should be hung below an edge that is above cluster c.
        return cs_out

    def remove_trivial(self):
        """ Remove the trivial clusters, i.e. the singletons """
        c_rem = [c for c in self.C if c.count('1') <= 1]
        for c in reversed(c_rem):
            self.C.remove(c)

    def collapse(self, leaves: list, label: str, convert2cl=False):
        """ Collapse the given leaves and give the new meta-leaf a new label """
        # Collapse X
        leaves = list(leaves)
        self.X = list(set(self.X).difference(leaves))
        self.X.append(label)
        self.X.sort()
        C_new = []
        for c in self.C:
            cl_new = [x for x in c if x not in leaves]
            # If some leaves in this cluster should be collapsed, then add the new label
            if set(c).intersection(leaves):
                cl_new.append(label)
                cl_new.sort()
            if len(cl_new) >= 2 and cl_new not in C_new:
                if convert2cl:
                    cl_new = Cluster(cl_new)
                C_new.append(cl_new)
        self.C = C_new

    def collapse_mst_sets(self):
        """
        Collapse all leaves that always occur together in the clusters
        :return:
        """
        # Get the non-singleton MST-sets
        mst_sets = self.get_mst_sets()
        if not mst_sets:
            return False  # Nothing collapsed
        for i, mst in enumerate(mst_sets):
            # Get a new meta-leaf label
            label = mst[0]
            # Collapse the nodes in cs
            self.collapse(mst, label, True)
        return True

    def remove_leaves(self, X):
        X = set(X)
        C_new = []
        for i, c in reversed(list(enumerate(self.C))):
            if X.intersection(c):
                c = set(c).difference(X)  # Remove X from cluster c
                if len(c) >= 2:  # Singletons are implicitly removed
                    c = sorted(c)
                    if c not in C_new:
                        C_new.append(c)
            elif c not in C_new:
                C_new.append(c)
        self.C = C_new
        self.X = sorted(set(self.X).difference(X))
        # There might be unused leaves now
        return self._remove_unused_leaves()

    def remove_leaf(self, x: str):
        C_new = []
        for c in self.C:
            # Remove the leaf from each cluster
            if x in c:
                if len(c) > 2:  # Neglect singeletons, i.e. effectively remove them
                    c.remove(x)
                    if c not in C_new:
                        C_new.append(c)
            elif c not in C_new:
                C_new.append(c)
        self.C = C_new
        self.X.remove(x)
        # There might be unused leaves now
        if not self.C:
            return []
        return self._remove_unused_leaves()

    def _remove_unused_leaves(self):
        """ Remove the unused nodes from X """
        X_unused = set(self.X)
        X_used = set()
        for c in self.C:
            X_unused = X_unused.difference(c)
            X_used.update(c)
        X_used = list(X_used)
        X_used.sort()
        self.X = X_used
        return X_unused

    def get_on_x(self, X):
        X_rem = set(self.X).difference(X)
        cs = self.copy()
        cs.remove_leaves(X_rem)
        cs.X = X
        return cs

    def _is_st_set(self, c):
        """
        Is the given cluster an ST-set?
        :param c: A cluster
        :return: T/F
        """
        # print(c)
        c = Cluster(c)
        for _c in self.C:
            if not c.iscompatible(_c):
                # Print.debug(f'{c} is incompatible with {_c}')
                return False
        # The ST-set S is compatible. Is C|S separating?
        return self.get_on_x(c).is_separating()

    def get_mst_sets(self, non_trivial_only=True):
        """
        Get the maximal ST-sets (which partition the taxa)
        :return: A list of the MST-sets
        """
        msts = [[x] for x in self.X]
        i = 0
        while i < len(msts):
            st1 = msts[i]
            j = i + 1
            while j < len(msts):
                st2 = msts[j]
                c = set(st1).union(st2)
                # c = st1 | st2  # Verify if this (cluster) is an ST-set
                if self._is_st_set(c):
                    # Replace the two ST-sets by their union
                    msts[i] = sorted(c)
                    del msts[j]
                    i = -1
                    break
                else:
                    # Print.debug(f'{c} is not an MST-set')
                    j += 1
            i += 1
        if non_trivial_only:
            return [s for s in msts if len(s) >= 2]
        return msts

    def get_st_sets2rem(self, msts_only=True, guess_mst=True, print_info=False):
        """
        Get all MST-sets (to remove).
        Use a heuristic to sort them.
        Given all MST-sets M, If there is an MST-set S such that:
         1. There is an S' in M: S -> S' and
         2. S" -/-> S, for all S" in M,
        then put it as far as possible in the beginning of the list. Using a depth-first search, using a LIFO-queue,
        these items would be investigated the latest.
        :return: A list of MST-sets
        """
        msts = self.get_mst_sets(False)
        if msts_only:
            if guess_mst:
                impls = self.get_implications()
                low_prio_msts = set(impls.keys())
                for to_msts in impls.values():
                    low_prio_msts = low_prio_msts.difference(to_msts)
                msts = sorted(msts, key=lambda mst: len(self.get_cs_in(mst[0]).C))  # For a LIFO queue
                # First sorting heuristic: sort MST-sets based on low prio flag
                if low_prio_msts:
                    return sorted(msts, key=lambda mst: not low_prio_msts.intersection(mst))  # if LIFO
                # Second sorting heuristic: sorting based on |C_x|
            return msts
        else:
            # Get the ST-sets
            sts = []
            for m in msts:
                if len(m) <= 2:
                    sts.append(m)
                else:
                    # Find the tree for this MST-set
                    cs_m = self.get_on_x(m)
                    if len(cs_m.C) == 1:
                        sts.append(m)
                        continue
                    # cs_m.len >= 2
                    tree = cs_m.to_nxtree()
                    r = tree.get_root()
                    ST = tree.get_st_sets(r)
                    sts.extend(ST)
            if guess_mst:
                impls = self.get_implications()
                low_prio_sts = set(impls.keys())
                for to_msts in impls.values():
                    low_prio_sts = low_prio_sts.difference(to_msts)
                return sorted(sts, key=lambda st: not low_prio_sts.intersection(st))  # If LIFO
                # return sorted(sts, key=lambda st: low_prio_sts.intersection(st))  # If FIFO
            return sts

    def get_implications(self):
        """
        If for each cluster containing leaf x, y is present too, then x -> y.
        :return: A dict of implications
        """
        impl = {}
        for x in self.X:
            s = set(self.X)
            s.remove(x)
            for c in self.C:
                if x in c:
                    s = s.intersection(c)
            if s:
                impl[x] = s
        return impl

    def to_nxtree(self):
        """
        Get the tree on X that represents the clusters C, as NetworkX multi digraph
        :return: A NetworkX digraph
        """
        # If some clusters are incompatible with each other, then no tree exists
        if not self.is_separating():
            return None

        # Sort the clusters from large to small
        C = list(self.C)
        C.sort(key=lambda c: len(c))
        tree = PhyloNet(nx.MultiDiGraph())
        tree.add_nodes_from(self.X)
        leaf_root = {x: x for x in self.X}  # Key: leaf. Value: the root of the connected component containing the leaf
        for c in C:
            r = tree.new_node_name()  # A new root
            heads = {leaf_root[x] for x in c}
            tree.add_edges_from((r, x) for x in heads)
            # Update leaf_root
            for x in c:
                leaf_root[x] = r
        roots = set(leaf_root.values())
        # If there are multiple roots, connect them to a new root
        if len(roots) >= 2:
            r = tree.new_node_name()  # A new root
            tree.add_edges_from((r, x) for x in roots)
        return tree

    def print(self, title='clusters', marked_leaf=None):
        """
        Print the clusters C by the leaf names
        :param title: A phrase to name the clusters
        :param marked_leaf: A leaf to highlight (optional)
        """
        # Init
        IWHITE = '\033[97m'
        ENDC = '\033[0m'
        X1 = list(self.X)
        if marked_leaf is not None:
            X1[X1.index(marked_leaf)] = IWHITE + marked_leaf + ENDC
        print(f'The {len(self.C)} {title} are:')
        i = 1
        for c in self.C:
            if marked_leaf and marked_leaf in c:
                c = c.copy()
                c[c.index(marked_leaf)] = IWHITE + marked_leaf + ENDC
            print('{:3d}. ('.format(i) + ', '.join(c) + ')')
            i += 1
