""" run python -m pyADAPT
usage:
=====
- pyADAPT analyze --sbml glycolysis.xml
- pyADAPT convert --sbml trehalose.xml
"""
import argparse
import click
from pyADAPT.bio.sbml_model import SBMLModel
from pyADAPT.adapt import ADAPT


def Main():
    parser = argparse.ArgumentParser(
        prog="pyADAPT",
        description="Run analysis of dynamic adaptations of parameter trajectories on a given SBML model.",
        allow_abbrev=True,
        add_help=True,
        epilog="For bug report or suggestions: <https://github.com/tonyxwz/pyADAPT>",
    )

    parser.add_argument(
        "job", choices=["analyze", "convert"], help="analyze or convert a model"
    )
    parser.add_argument(
        "-f",
        "--sbml",
        metavar="SBML",
        required=True,
        dest="sbml",
        help="a sbml file that contains the model to by analyzed",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        metavar=("p1", "p2"),
        action="append",
        nargs="+",
        required=True,
        help='the parameters to "ADAPT"',
    )
    parser.add_argument(
        "--gui", action="store_true", help="use command line interface only"
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="don't print result to console"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="results.txt",
        action="store",
        help="file to write analysis",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    args = parser.parse_args()
    print(args)
    # pass argument to the analysis or convert routine


if __name__ == "__main__":
    Main()
