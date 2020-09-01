# -*- coding: utf-8 -*-
import os
from keyword import iskeyword
import importlib
import re
import jinja2
from jinja2 import PackageLoader, Environment
from pyADAPT.sbml.sbml_model import SBMLModel
import numpy as np


class Converter(object):
    """ pyADAPT convert handler
    - pyADAPT convert --sbml XXX.sbml --parameters p1 p2 --output mymodel.py
    """
    def __init__(self, sbml, params=[], output=""):
        # since I already have the SBMLModel implemented
        # self.fsloader = jinja2.FileSystemLoader(searchpath="./templates")
        self.pkgloader = jinja2.PackageLoader(package_name="pyADAPT",
                                                package_path="templates")
        self.sbml_model = SBMLModel(sbml)
        self.environ = jinja2.Environment(loader=self.pkgloader)

        self.params = params
        self.change_to_keyword = lambda varStr: "Model" + re.sub(
            r"\W|^(?=\d)", "_", varStr)

        if len(output):
            if output[-3:] != ".py":
                output = output + ".py"
        else:
            # only replace the last instance
            output = sbml[::-1].replace(".xml"[::-1], ".py"[::-1], 1)[::-1]
        self.output = output

    def model_expr(self):
        pass

    def convert(self):
        tmplt = self.environ.get_template("model.py.j2")

        if iskeyword(self.sbml_model.name):
            model_name = self.sbml_model.name
        else:
            model_name = self.change_to_keyword(self.sbml_model.name)

        model_specs = dict()
        model_specs["model_name"] = model_name
        model_specs["model_description"] = self.sbml_model.notes
        # model_specs["predictor"] = None

        model_string = tmplt.render({
            "model_name": model_name,
            "model_description": self.sbml_model.notes,
            "model": self.sbml_model,
            "sm": np.array2string(self.sbml_model.stoich_matrix, separator=","),
            "initial_states": np.array2string(self.sbml_model.initial_states, separator=",")
        })

        with open(self.output, "w") as f:
            f.write(model_string)

        spec = importlib.util.spec_from_file_location(model_name, self.output)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        # should I return the imported python module?
        return getattr(module, model_name)


if __name__ == "__main__":
    converter = Converter(r"data/trehalose/smallbone.xml")
    print(converter.convert())
