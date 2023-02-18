"""
This script is created to search for a level-k phylogenetic network with the lowest level possible and with a low
number of reticulations (hybrid nodes).

Required packages:
- biopython
- networkx-v2.8.3 (it does not work well with v2.8.4)
"""
import argparse
from pathlib import Path

import networkx as nx

from Cass import draw_multi, FILTER_X, GUESS_X, GUESS_EDGES, FILTER_EDGES, optcass, \
    STAY_BIN, MIN_HANG_ONLY
from ClusterSet import ClusterSet
from PhyloHelper import read, get_clusters, get_term_names

proc_times = {}
start_times = {}
TRACE = 1
BENCHMARK = 2


def _check_reqs():
    """
    Check the package requirements
    """
    if nx.__version__ == '2.8.4':
        ImportError('Package NetworkX v2.8.4 contains a bug. I recommend to install v2.8.3.\n'
                    "Hint: use pip3 install networkx==2.8.3'")


def _print_line(*args):
    print('{:12} {:>8} {:>4} {:>4} {:>4}'.format(*args))


def _main():
    """
    The main program
    """
    _check_reqs()

    # Give the input parameters: filename, k_max
    filename = "../data/eluFig9.nwk"  # k=3 Works for k=4 20220706. Works for k=3 20221003!

    # Parse the input arguments
    parser = argparse.ArgumentParser(
        prog="STCASS Algorithm",
        description="Determine a lowest-level network for the clusters in the given trees."
    )
    parser.add_argument('-file', type=str, nargs='?', default=filename,
                        help="A filename. The file should contain trees formatted in the Newick format.")
    parser.add_argument('-timeout', type=int, nargs='?', default=5,
                        help="The maximal runtime in minutes.")
    parser.add_argument('-sys_out', type=bool, nargs='?', default=False,
                        help="Give the results in an easily interpretable format for a computer.")
    args = parser.parse_args()
    FILE_NAME = args.file
    TIMEOUT = args.timeout * 60  # Timeout in seconds
    SYS_OUT = args.sys_out

    # Read the input file
    trees = read(FILE_NAME)
    # short_name = Path(filename).stem
    # draw_trees(trees, 'Input Trees', filename=(short_name + '_trees'))

    clusters, _, _ = get_clusters(trees)
    leaves = get_term_names(trees)
    cs = ClusterSet(clusters, leaves)
    strategies = {
        # 'STCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | ALL_STSETS | MIN_HANG_ONLY | FULL_CX,
        # 'STCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | MIN_HANG_ONLY | FULL_CX,
        'SuperCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | MIN_HANG_ONLY,
    }
    results = []
    nets = []
    for label, strategy in strategies.items():
        _nets, t, k = optcass(cs, strategy, TIMEOUT)
        nets.extend(_nets)
        if _nets:
            net0 = _nets[0]
            r = net0.get_ret_num()
        else:
            r = 0
        n = len(_nets)
        keypoints = {
            'strategy': strategy,
            'alg': label,
            't': '{:.3f}'.format(t),  # Runtime [s]
            'k': str(k),  # Level
            'r': str(r),  # #Reticulations
            'nSols': str(n),  # #Solutions
            'clAmt': len(clusters),  # #Clusters
            'leafAmt': len(leaves)  # #Leaves
        }
        # If a timeout occurred, log it
        if t > TIMEOUT:
            keypoints['timeout'] = True
        results.append(keypoints)

    print("\nRESULTS")
    _print_line('Algorithm', 't [s]', 'k', 'r', 'n')
    print('-' * 36)
    for result in results:
        if SYS_OUT:
            print(','.join(f'{k}={v}' for k, v in result.items()))
        else:
            _print_line(*list(result.values())[1:])

    short_name = Path(filename).stem
    filename_part = short_name.split('.')[0]
    # Show the result
    # draw_multi(nets, f'A Level-{results[0]["k"]} Network for {filename_part}', full_screen=True, all_labels=False)
    # file=f'{filename_part}_net')

    # Save the network in Newick format
    # nets[0].save_net(short_name + '_net.nwk')


if __name__ == '__main__':
    _main()
