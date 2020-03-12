""" This module provides some pre-defined regularization functions for ADAPT.

A consistent calling signature for regularization functions should be defined.
The objective function have no idea what a regularization function needs, but
it can offer all the knowledge it have.
    - the parameter trajectory so far, `parameter_trajectory`
    - the state trajectory so far, `state_trajectory`
    - which iteration number, i_iter
    - which time step is it, i_ts
    - errors between data and prediction, error
    - time span, if regularization function want to be time dependent, time_span
    - T.B.C
"""

TIEMANN = 0


def tiemann(params=None,
            parameter_trajectory=None,
            state_trajectory=None,
            time_span=None,
            i_iter=None,
            i_ts=None,
            errors=None,
            **kw):
    pass
