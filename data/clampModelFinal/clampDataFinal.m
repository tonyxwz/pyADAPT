function DATASET = Data()

DATASET.DESCRIPTION = 'Clamp data.';

% DATASET.FILE = 'Data2test';
DATASET.FILE = 'DataFinal2';


DATASET.GROUPS = {
    'BC003pre'
    'BC003post'
    'BC004pre'
    'BC004post'
};

DATASET.FIELDS = {
    % OUTPUTS
    'glucose'           1       'glu_time' 'glu_conc' 'glu_std' 1   [] % Measured plasma glucose concentration
    'glucose_TTR'       1       'glu_time2' 'glu_ttr' 'glu_ttr_std' 1   [] %
    
    'glycerol'          1       'gly_time' 'gly_conc' 'gly_std' 1   [] %
    'glycerol_TTR'      1      'gly_time2' 'gly_ttr' 'gly_ttr_std' 1   [] %
    
    'insulin'           1       'ins_time' 'ins_conc' 'ins_std' 1   [] %
    'FFA'               1       'FFA_time' 'FFA_conc' 'FFA_std' 1   [] %
%     'cpeptide'  1       'cpe_time' 'cpe_conc' 'cpe_std' 1   [] %
%     'cho'       1       'cho_time' 'cho_conc' 'cho_std' 1   [] %
%     'ldl'       1       'ldl_time' 'ldl_conc' 'ldl_std' 1   [] %
%     'hdl'       1       'hdl_time' 'hdl_conc' 'hdl_std' 1   [] %
%     'tri'       1       'tri_time' 'tri_conc' 'tri_std' 1   [] %    
    
    'Gb'                0       []          'Gb'        'Gb'    1   [] %
    'Yb'                0       []          'Yb'        'Yb'    1   [] %
    'Ib'                0       []          'Ib'        'Ib'    1   [] %
    'FFAb'              0       []          'FFAb'      'FFAb'  1   [] %
    
    'TTR_VAR'           0       []          'glu1_ttr'  'glu1_ttr' 1 []
    'TTR_FIX'           0       []          'glu2_ttr'  'glu2_ttr' 1 []
    
    'BW'                0       []          'BW'        'BW'       1 [] %        
    
    % INPUTS
    'ins1'              1       'ins1_time' 'ins1_flux' 'ins1_std' 1 [] %
    'glu1'              1       'glu1_time' 'glu1_flux' 'glu1_std' 1 [] %
    'glu2'              1       'glu2_time' 'glu2_flux' 'glu2_std' 1 [] %
    'glu3'              1       'glu3_time2' 'glu3_flux2' 'glu3_std2' 1 [] %
    'gly1'              1       'gly1_time' 'gly1_flux' 'gly1_std' 1 [] %
    'gly2'              1       'gly2_time2' 'gly2_flux2' 'gly2_std2' 1 [] %
};
