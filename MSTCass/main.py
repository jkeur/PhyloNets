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
from PhyloHelper import read


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

    # Parse the input arguments
    parser = argparse.ArgumentParser(
        prog="STCASS Algorithm",
        description="Determine a lowest-level network for the clusters in the given trees."
    )
    parser.add_argument('-infile', type=str, nargs='?', default='ElusivenessFig9.nwk',
                        help="The filename of the input file. The file should contain trees formatted in the Newick "
                             "format or it should contain a set of clusters.")
    parser.add_argument('-outfile', type=str, nargs='?', default='',
                        help="Save the result (in the Newick format).")
    parser.add_argument('-fig_file', type=str, nargs='?', default='',
                        help="Save the resulting network as figure.")
    parser.add_argument('-timeout', type=int, nargs='?', default=5,
                        help="The maximal runtime in minutes.")
    parser.add_argument('-sys_out', type=bool, nargs='?', default=False,
                        help="Give the results in an easily interpretable format for a computer.")
    parser.add_argument('-show', type=bool, nargs='?', default=False,
                        help="Show the result in a figure.")
    args = parser.parse_args()
    IN_FNAME = args.infile
    OUT_FNAME = args.outfile
    FIG_FNAME = args.fig_file
    TIMEOUT = args.timeout * 60  # Timeout in seconds
    SYS_OUT = args.sys_out
    SHOW = args.show

    # Read the input file
    clusters, leaves = read(IN_FNAME)
    cs = ClusterSet(clusters, leaves)

    strategies = {
        # 'STCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | ALL_STSETS | MIN_HANG_ONLY | FULL_CX,
        # 'STCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | MIN_HANG_ONLY | FULL_CX,
        'MSTCass': FILTER_X | GUESS_X | FILTER_EDGES | GUESS_EDGES | STAY_BIN | MIN_HANG_ONLY,
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

    if SHOW or len(FIG_FNAME):
        # Show the result
        short_name = Path(IN_FNAME).stem
        filename_part = short_name.split('.')[0]
        draw_multi(nets, f'A Level-{results[0]["k"]} Network for {filename_part}', file=FIG_FNAME, full_screen=True,
                   all_labels=False, show=False)

    if OUT_FNAME and len(OUT_FNAME) > 0:
        # Save the network in Newick format
        nets[0].save_net(OUT_FNAME)


if __name__ == '__main__':
    _main()
