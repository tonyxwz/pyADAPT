# pyADAPT

Analysis of Dynamic Adaptations in Parameter Trajectories

[http://bmi.bmt.tue.nl/sysbio/software/adapt.html](http://bmi.bmt.tue.nl/sysbio/software/adapt.html)

Implementation in Python

## Getting Started

Clone this repo, checkout to branch stable from where you would be able to make changes suitable for your own model.

Install pyADAPT to the current python environment: `pip install -e .`

### Define model

Convert sbml file to model definition script: `python -m pyADAPT convert -f yourmodel.xml -o yourmodel.py`

```python
from yourmodel import YourModel
model = YourModel()
```
Alternatively you could also define the models directly:

- Inherit from `pyADAPT.basemodel.BaseModel`
- `add_parameter` according to your model in method `__init__`
- define methods `fluxes` and `state_ode`

### Define dataset

The dataset is required to be several time-dependent pairs of mean values and standard deviations.

```python
data = pyADAPT.dataset.Dataset(raw_data_path="toyData.mat", data_specs_path="toyData.yaml")
```

### Run the simulation

Analyze parameter 'k1' of the model using 4 processes, ODE solver "LSODA":

```python
from pyADAPT.optimize import optimize

(parameter_trajectories,
  state_trajectories,
  flux_trajectories,
  time) = optimize(model, dataset, 'k1',
                   n_iter=256,
                   delta_t=0.2,
                   odesolver="LSODA", 
                   n_core=4)
```

## Analysis of the parameter trajectories

Please refer to https://research.tue.nl/en/publications/computational-analysis-of-adaptations-during-disease-and-interven
