import os
from keyword import iskeyword
import importlib
import re
import jinja2
from jinja2 import PackageLoader, Environment
from pyADAPT.bio.sbml_model import SBMLModel


class Converter(object):
    """ pyADAPT convert handler
    - pyADAPT convert --sbml XXX.sbml --parameters p1 p2 --output mymodel.py
    """

    def __init__(self, sbml, params=[], output=""):
        # TODO: consider Jinja objects:  http://zetcode.com/python/jinja/
        # since I already have the SBMLModel implemented
        # self.fsloader = jinja2.FileSystemLoader(searchpath="./templates")
        self.pkgloader = jinja2.PackageLoader(
            package_name="pyADAPT", package_path="templates"
        )
        self.sbml_model = SBMLModel(sbml)
        self.environ = jinja2.Environment(loader=self.pkgloader)

        self.params = params

        if len(output):
            if output[-3:] != ".py":
                output = output + ".py"
        else:
            output = sbml[::-1].replace(".xml"[::-1], ".py"[::-1], 1)[::-1]
        self.output = output
        self.change_to_keyword = lambda varStr: "Model" + re.sub(
            r"\W|^(?=\d)", "_", varStr
        )

    def convert(self):
        tmplt = self.environ.get_template("model.py.j2")

        if iskeyword(self.sbml_model.name):
            model_name = self.sbml_model.name
        else:
            model_name = self.change_to_keyword(self.sbml_model.name)

        model_specs = dict()
        model_specs["model_name"] = model_name
        model_specs["model_description"] = self.sbml_model.notes
        model_specs["predictor"] = None

        model_string = tmplt.render(
            {
                "model_name": model_name,
                "model_description": self.sbml_model.notes,
                "predictors": [
                    {"name": "g2", "value": 63},
                    {"name": "t4", "value": 56},
                ],
                "constants": [{"name": "u1", "value": 1}, {"name": "u2", "value": 3}],
                "parameters": [{"name": "k1", "value": 4, "vary": True, "lb": 0}],
                "states": [{"name": "s1", "init": 1}, {"name": "s2", "init": 2}],
            }
        )

        with open(self.output, "w") as f:
            f.write(model_string)

        spec = importlib.util.spec_from_file_location(model_name, self.output)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # should I return the imported python module?
        return getattr(module, model_name)


if __name__ == "__main__":
    converter = Converter(r"data\trehalose\smallbone.xml")
    print(converter.convert())
