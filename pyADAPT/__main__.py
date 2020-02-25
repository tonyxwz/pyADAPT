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
        description="Run analysis of dynamic adaptations of parameter trajectories on a given SBML model.",
        allow_abbrev=True,
        add_help=True,
        epilog="For bug report or suggestions: <https://github.com/tonyxwz/pyADAPT>",
    )

    parser.add_argument(
        "job",
        choices=["analyze", "convert", "new", "gui"],
        help="analyze/convert a model or start the GUI",
    )  # TODO write a model cookiecutter
    parser.add_argument(
        "-f",
        "--sbml",
        metavar="SBML",
        dest="sbml",
        help="a sbml file that contains the model to by analyzed",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        metavar=("p1", "p2"),
        nargs="+",
        help='the parameters to "ADAPT"',
    )
    parser.add_argument(
        "--quiet", action="store_true", help="don't print result to console"
    )
    parser.add_argument(
        "--output",
        metavar="results.txt",
        action="store",
        help="file to write analysis",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    print(sys.argv)
    args = parser.parse_args()
    print(args)

    if args.job != "gui":
        if (args.parameters is None) or (args.sbml is None):
            print(
                "error:",
                f"{parser.prog} {args.job} requires at least a model file AND a parameter list",
            )
            parser.print_help()
            sys.exit()
        # pass argument to the analysis or convert routine
        elif args.job == "convert":
            print("converting")
        elif args.job == "analyze":
            print(f"analyzing {args.parameters} of {args.sbml}")
        else:
            raise Exception(f"error: pyADAPT: unknown task {args.job}")
    else:
        print("start gui")


if __name__ == "__main__":
    Main()
