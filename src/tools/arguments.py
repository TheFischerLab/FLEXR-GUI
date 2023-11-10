import argparse

def create_parser():
    """Arguments for FLEXR"""

    CLI = argparse.ArgumentParser()

    CLI.add_argument(
        '-m',
        '--mode',
        nargs="?",
        type=str,
        default='FLEXR',
        help='Mode: FLEXR, DELTA, or FLEXRSCORE. Default: FLEXR'
    )

    CLI.add_argument(
        '-f',
        '--filename',
        nargs="?",
        type=str,
        default=None,
        help='Input ringer.csv file for FLEXR/FLEXRSCORE'
    )

    CLI.add_argument(
        '-f2',
        '--filename2',
        nargs="?",
        type=str,
        default=None,
        help='Input ringer.csv file for FLIPPER'
    )

    CLI.add_argument(
        '-g',
        '--geotolerance',
        type=int,
        default=30,
        help='Tolerance for match between measured chi and ideal \
        chi in library for FLEXR/FLEXRSCORE. Default = 30'
    )

    CLI.add_argument(
        '-t',
        '--sigmathreshold',
        nargs="?",
        type=float,
        default=0.3,
        help='Sigma threshold for peak finding in FLEXR/FLEXRSCORE. Default = 0.3'
    )

    CLI.add_argument(
        '-ph',
        '--peakheight',
        type=float,
        default=0.03,
        help='Required height of peaks in FLEXR peak finding. Default = 0.03'
    )

    CLI.add_argument(
        '-pp',
        '--peakprominence',
        type=float,
        default=0.05,
        help='Required prominence of peaks in FLEXR peak finding. Default = 0.05'
    )

    CLI.add_argument(
        '-pw',
        '--peakwidth',
        type=float,
        default=1,
        help='Required width of peaks in FLEXR peak finding. Default = 1'
    )

    CLI.add_argument(
        '-pd',
        '--peakdistance',
        type=int,
        default=5,
        help='Required minimal horizontal distance between neighboring peaks \
        in FLEXR peak finding. Default = 5'
    )

    CLI.add_argument(
        '-p',
        '--plot',
        type=str,
        default=False,
        help='Save individual plots showing peak finding results? This is slow. \
        Default = False'
    )

    CLI.add_argument(
        '-s',
        '--step',
        nargs=1,
        type=int,
        default=5,
        help='Step sized used for sigma measurements in ringer.csv file.'
    )

#    CLI.add_argument(
#        '-fs',
#        '--flexrscore',
#        nargs=1,
#        type=bool,
#        default=False,
#        help='CRYO-EM FLEXR: run FLEXR-SCORE to find optimal map value threshold. \
#        Default: False'
#    )

    CLI.add_argument(
        '-l',
        '--build_limit',
        type=int,
        default=3,
        help='Limit the number of rotamers FLEXR will build per residue. Default: 3'
    )

    CLI.add_argument(
        '-li',
        '--ligand',
        type=str,
        default=None,
        help='Define ligand to restrict FLEXR building around. Example: "39 A"'
    )

    CLI.add_argument(
        '-d',
        '--distance',
        nargs=1,
        type=float,
        default=5.0,
        help='Distance for ligand restricted building in FLEXR. Default: 5.0'
    )

    CLI.add_argument(
        '-pdb',
        '--pdb_in',
        type=str,
        default=None,
        help='PDB file that corresponds with ringer.csv file.'
    )

    CLI.add_argument(
        '-sconfs',
        '--singleconfs',
        type=str,
        default=False,
        help='Output single conformers detected in FLEXR. Default = False'
    )

### Ringer-delta options

    CLI.add_argument(
        '-chi',
        '--chi',
        nargs="?",
        type=str,
        default='chi1',
        help='What dihedral angle do you want to study with DELTA? chi1, chi2, chi3,etc? \
        Example: -chi chi1 . Default: chi1 '
    )

    CLI.add_argument(
        '-colors',
        '--colors',
        nargs="?",
        type=str,
        default='blue',
        help='Colors to plot each model in DELTA - files are read in alphabetically. \
        Example: -colors blue,red \
        Default: random assignment'
    )

    CLI.add_argument(
        '-safety',
        '--safety',
        nargs=1,
        type=str,
        default=True,
        help='Safety: compare only matched residues that are either branched or unbranched \
        in DELTA. Default: True'
    )

    CLI.add_argument(
        '-reload',
        '--reload',
        nargs=1,
        type=str,
        default=False,
        help='Use manually tweaked alignment file in DELTA? Default: False'
    )

    CLI.add_argument(
        '-pearson',
        '--pearson',
        nargs=1,
        type=str,
        default=False,
        help='Calculate Pearson CC values in DELTA. Default: False'
    )

    CLI.add_argument(
        '-render',
        '--render',
        nargs=1,
        type=str,
        default=False,
        help='Add median CC value to B-factor column in first file to render in PyMOL \
        in DELTA. Default: False'
    )

## flexr_building options

    CLI.add_argument(
        '-build',
        '--build',
        type=str,
        default=True,
        help='add alt confs to model'
    )

    CLI.add_argument(
        '-branching',
        '--branching',
        type=str,
        default='CA',
        help='where do you want alt confs to branch in model? "CA" or "ALL". Default: "CA" '
    )

    CLI.add_argument(
        '-exitcoot',
        '--exitcoot',
        type=str,
        default=True,
        help='Close Coot when finished buildng?'
    )

    CLI.add_argument(
        '-cootmolnum',
        '--cootmolnum',
        type=int,
        default=0,
        help='Mol num in Coot for building. Default: 0'
    )

    return CLI
