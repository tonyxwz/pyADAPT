# The Python GUI for ADAPT

This GUI should first show a file pick window for user to choose a SBML model
file. After the path to the SBML file is given, it call `pyADAPT` to analyze
the model and give a list of all the parameters in the model in the next
windows.

## Implementation

Framework: PyQt5

## Components

MainWindow:

- `__init__` should first show a dialogue, asking user to choose a SBML model
	(call QFileDialogue)
- Process the SBML file and show the parameters in a list (maybe in a simular
	way to Pascal's MATLAB GUI), this is the MainWindow.
- Two buttons

	- Convert: with a check box: save to the same directory as model file
	- Analysis: call pyADAPT.analyze.analyze


