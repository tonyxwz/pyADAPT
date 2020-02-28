""" This module provides some pre-defined regularization functions for ADAPT
TODO a consistent calling signature for regularization functions
the objective have no idea what a regularization function needs, but it can
offer all the knowledge it have.
    - the parameter trajectory so far, `parameter_trajectory`
    - the state trajectory so far, `state_trajectory`
    - which iteration number, i_iter
    - which time step is it, i_ts
    - errors between data and prediction, error
    - time span, if regularization function want to be time dependent, time_span
    - T.B.C
"""
# TODO integrate with ADAPT

TIEMANN = 0


def tiemann():
    pass
