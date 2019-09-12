function MODEL = toyModelModDavid(~)

MODEL.DESCRIPTION = 'Toy model created by N.A.W van Riel.';

MODEL.PREDICTOR = {
    't'  [0 10] {'time' 'days' 'Time'}
};

MODEL.CONSTANTS = {
    'u1' 1 [] {}
    'u2' 1 [] {}
};

MODEL.PARAMETERS = {
    'k1' 1  [] {}
    'k2' 1  [] {}
    'k3' .1 [] {}
    'k4' .5 [] {}
    'k5' 1  [] {}
};

MODEL.STATES = {
    's1' 1.03 's1' 'v1 - v3 - v4' {}
    's2' 0.38 's2' '-v1 + v2'     {}
    's3' 0.62 's3' 'ds3dt'        {}
    's4' 0.52 's4' 'v4 - v5'      {}
};

MODEL.REACTIONS = {
    'v1' 'k1 * u1 * s2' {}
    'v2' 'k2 * u2 * s3' {}
    'v3' 'k3 * s1'      {}
    'v4' 'k4 * s1'      {}
    'v5' 'k5 * s4'      {}
    
    'ds3dt' 'v1 - v2'   {}
};