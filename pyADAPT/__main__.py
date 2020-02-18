""" run python -m pyADAPT
usage:
=====
- pyADAPT analyze --sbml glycolysis.xml
- pyADAPT convert --sbml trehalose.xml
"""
import sys
import argparse
import click
from pyADAPT.bio.sbml_model import SBMLModel
from pyADAPT.adapt import ADAPT
from pyADAPT.optimize import optimize


def Main():
    parser = argparse.ArgumentParser(
        prog="pyADAPT",
        description=
        "Run analysis of dynamic adaptations of parameter trajectories on a given SBML model.",
        allow_abbrev=True,
        add_help=True,
        epilog=
        "For bug report or suggestions: <https://github.com/tonyxwz/pyADAPT>",
    )

    parser.add_argument("job",
                        choices=["analyze", "convert", "gui"],
                        help="analyze/convert a model or start the GUI")
    parser.add_argument(
        "-f",
        "--sbml",
        metavar="SBML",
        dest="sbml",
        help="a sbml file that contains the model to by analyzed",
    )
    parser.add_argument(
        "--parameters",
        metavar=("p1", "p2"),
        action="append",
        nargs="+",
        help='the parameters to "ADAPT"',
    )
    parser.add_argument("--quiet",
                        action="store_true",
                        help="don't print result to console")
    parser.add_argument(
        "--output",
        metavar="results.txt",
        action="store",
        help="file to write analysis",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    args = parser.parse_args()

    if (args.job != "gui"):
        if ((args.parameters is None) or (args.sbml is None)):
            print(
                "error:",
                f"{parser.prog} {args.job} requires at least a model file AND a parameter list"
            )
            parser.print_help()
            sys.exit()
    else:
        # gui will simply prepare the paremters for minimize or convert
        print("start gui")
    if args.job == "convert":
        print("converting")
    elif args.job == "analysis":
        print("analyzing")
    else:
        raise Exception(f"error: pyADAPT: unknown task {args.job}")
    # prepare arguments and call optimize here
    print(args)
    # pass argument to the analysis or convert routine


if __name__ == "__main__":
    Main()
