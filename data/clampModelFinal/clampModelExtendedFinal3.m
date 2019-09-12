function MODEL = Model(~)

MODEL.DESCRIPTION = 'Model, clamp';

MODEL.PREDICTOR = {
    't'  [0 385] {'time' 'min' 'Time'}
%     't'  [0 1000] {'time' 'min' 'Time'}

};

MODEL.CONSTANTS = {
    'Gb'                3.9233      'Gb'                {}
    'Ib'                23.3333     'Ib'                {}
    'MM_glucose'        180.182     []                  {} % Molecular Weight Glucose
    'TTR_fix'           1.0000      'TTR_FIX'           {} % Tracer to tracee ratio of the bolus (so 100% isotope)
    'TTR_var'           0.0087      'TTR_VAR'           {} % Tracer to tracee ratio of the variable glucose infusion
    'BW'                121.0000    'BW'                {} % Body Weight
    'FFAb'              920.0000    'FFAb'              {} % Basal (?) level of FFA
    'v_g'               0.150       []                  {} % Liter / kg (anders per substanties, anders voor gluc en insuline) % study shows 200ml/kg, so 0.2. Other studies show somewhere between 0.1 and 0.16
%     'v_g'               0.04       []                  {} % Liter / kg (anders per substanties, anders voor gluc en insuline) % study shows 200ml/kg, so 0.2. Other studies show somewhere between 0.1 and 0.16
    'v_i'               0.10        []                  {} % 15.6 / 67.4 , or according to Rang, H. P. (2003). Pharmacology. Edinburgh: Churchill Livingstone. 0.05 - 0.1
    'v_y'               0.340       []                  {} % https://www.ncbi.nlm.nih.gov/pubmed/7316993 says 34%, https://www.ncbi.nlm.nih.gov/pubmed/3585175 31 %
    'v_f'               0.10        []                  {} % http://diabetes.diabetesjournals.org/content/52/7/1641 says 55 ml / fat free kg, i'll say that it is 0.1 for most of the patients (rought estimate)
};

MODEL.PARAMETERS = {
    'p0'            0.035                       []  {}
% %     'p0'            'p0_data'                       'p0_data'  {}
%     'p1'            0.5                     []  {}
    'p2'            0.0001                  []  {}
%     'p3'            0.1                  []  {}
    'p4'            0.5                  []  {}
%     'p5'            0.25                    []  {}
%     'p6'            0.001                  []  {}
%     'p7'            0.04                    []  {}
    'p8'            350                     []  {}
%     'p9'            1                       []  {}
    'p10'           0.07                  []  {}
%     'p11'           0.02                    []  {}
    'p12'           4                     []  {}
    'p13'           0.09                     []  {}
    'p14'           0.4                     []  {}
    'p15'           0.4                     []  {}
%     'p16'           0.4                     []  {}
    'k1'            0.004                       []  {}
    'k2'            0.0005                       []  {}
    'k3'            22                       []  {}
    'k4'            1.3                       []  {}
    };

% MODEL.PARAMETERS = {
%     'p0'            0.03                       []  {}
%     'p1'            0.5                     []  {}
%     'p2'            0.0001                  []  {}
%     'p3'            0.1                  []  {}
%     'p4'            20                  []  {}
%     'p5'            0.25                    []  {}
%     'p6'            0.001                  []  {}
%     'p7'            0.04                    []  {}
%     'p8'            14                      []  {}
%     'p9'            1                       []  {}
%     'p10'           0.0025                  []  {}
%     'p11'           0.02                    []  {}
%     'p12'           5                     []  {}
%     'p13'           0.085                      []  {}
%     'p14'           0.05                     []  {}
%     'k1'            0.015                       []  {}
%     'k2'            0.01                       []  {}
%     'k3'            23                       []  {}
%     'k4'            5                       []  {}
%     };

MODEL.INPUTS={ % the glucose infusions are in mg/kg/min for glu1 and glu2, and mg/kg for glu3. 
    'glu1'  'Data'  {[0 125-1e-9 125 148-1e-9 148 158-1e-9 158 168-1e-9 168 178-1e-9 178 207-1e-9 207 228-1e-9 228 245-1e-9 245 255-1e-9 255 272-1e-9 272 288-1e-9 288 298-1e-9 298 307-1e-9 307 385-1e-9 385],   [0 0 1.0006 1.0006 1.2639 1.2639 2.0011 2.0011 2.5014 2.5014 3.5019 3.5019 3.0017 3.0017 2.5014 2.5014 2.7647 2.7647 3.5019 3.5019 4.0022 4.0022 5.0028 5.0028 6.0297 6.0297 6.5299 6.5299 0]}              {'glu1'}     'LINEAR' {'conc.' '?mol/min' 'Glucose'}
    'glu2'  'Data'  {[0 385-1e-9 385],                                                      [0.0200 0.0200 0]}                                      {'glu2'}     'LINEAR' {'conc.' '?mol/min' 'Glucose'} 
    'glu3'  'Data'  {[0 1-1e-9 1 385],                                                      [2.0020 2.0020 0 0]}                                    {'glu3'}     'LINEAR' {'conc.' '?mol/min' 'Glucose'}  
    'gly1'  'Data'  {[0 385-1e-9 385],                                                      [0.1100 0.1100 0]}                                      {'gly1'}     'LINEAR' {'conc.' '?mol/min' 'Glycerol'}  
    'gly2'  'Data'  {[0 1-1e-9 1 385],                                                      [1.6000 1.6000 0 0]}                                    {'gly2'}     'LINEAR' {'conc.' '?mol/min' 'Glycerol'}   
    'ins1'  'Data'  {[0 125-1e-9 125 265-1e-95 265 385-1e-9 385],                           [0 0 280.0000 280.0000 830.0000 830.0000 0]}            {'ins1'}     'LINEAR' {'conc.' '?mol/min' 'Insulin'}
};

MODEL.STATES={
    'glucose_tracer'        0           []      '(- glu_upt_Id - glu_upt_nId)*glucose_TTR + glu_inf_tracer'                              {} % molarity? -> milimol / liter
%     'glucose_tracer'        1000           []      '(- glu_upt_Id - glu_upt_nId)*glucose_TTR + glu_inf_tracer'                              {} % molarity? -> milimol / liter

    'glucose_tracee'        3.9233      'Gb'    '(- glu_upt_Id - glu_upt_nId)*(1-glucose_TTR) + glu_prod_nId + glu_inf_tracee'               {} % molarity? -> milimol / liter
    'insulin'               23.3333     'Ib'    '-I_deg + I_prod + I_inf'                                                               {} % molarity? -> picomol / liter
    'insulinAct'            23.3333     'Ib'    'p14 *(insulin - insulinAct)'                                              {}
    'insulinAct2'           23.3333     'Ib'    'p15 *(insulin - insulinAct2)'                                              {}
%     'insulinAct3'           23.3333     'Ib'    'p16 *(insulin - insulinAct3)'                                              {}
    'FFA'                   920.0000    'FFAb'  'FFA_prod_Id - FFA_upt_Id + FFA_prod_nId'                                {} % molarity? -> micromol / liter                                  
%     'FFA'                   0    'FFAb'  'FFA_prod_Id - FFA_upt_Id + FFA_prod_nId'                                {} % molarity? -> micromol / liter                                  

    'glycerol_tracer'       0           []      '(-gly_upt_nId)*glycerol_TTR + gly_inf_tracer'                                              {} % molarity? -> micromol / liter
    'glycerol_tracee'       151         'Yb'    '(-gly_upt_nId)*(1-glycerol_TTR) + gly_prod_Id + gly_prod_nId + gly_inf_tracee'                            {} % molarity? -> micromol / liter
    'bolus_total'       	0           []      'glu3*BW/((MM_glucose+2)*V_G)'                            {} % molarity? -> micromol / liter
%     'bolus_total'       	1000        []      'glu3*BW/((MM_glucose+2)*V_G)'                            {} % molarity? -> micromol / liter

};

MODEL.REACTIONS = {
    % Distribution volumes
    'V_G'               'v_g*BW'    {}
    'V_I'               'v_i*BW'    {}
    'V_Y'               'v_y*BW'    {}
    'V_F'               'v_f*BW'    {}
    
    % Exogenous inputs
    'glucose_VAR'       'glu1'                                                                                              {}
    'glucose_FIX'       'glu2'                                                                                              {}
    'glucose_BOLUS'     'glu3'                                                                                              {}
    'glycerol_FIX'      'gly1'                                                                                              {}
    'glycerol_BOLUS'    'gly2'                                                                                              {}
    
    % TTR
    'glucose_TTR'       'glucose_tracer     / (glucose_tracer  + glucose_tracee)'                                           {}
    'glycerol_TTR'      'glycerol_tracer    / (glycerol_tracer + glycerol_tracee)'                                          {}
    
    % Total concentrations
    'glucose'           'glucose_tracee + glucose_tracer'                                               {} % molarity? -> milimol / liter
    'glycerol'          'glycerol_tracee + glycerol_tracer'                                             {} % molarity? -> micromol / liter
    
    'glu_prod_nId'      'p0'                                                             {} % milimol / liter % gluconeogenesis and glycogenolysis                                                                          {} % milimol / liter % gluconeogenesis and glycogenolysis
    'glu_upt_Id'        'p2*insulinAct*glucose'                                                         {} % milimol / liter   
    'glu_upt_nId'       'k1*glucose'    {}
    
    'glu_inf_tracer'    'BW/((MM_glucose+2)*V_G)*(glu1*(TTR_var)+glu2*(TTR_fix)+glu3*(TTR_fix))'            {} % milimol / liter
    'glu_inf_tracee'    'BW/(MM_glucose*V_G)*(glu1*(1-TTR_var)+glu2*(1-TTR_fix)+glu3*(1-TTR_fix))'      {} % milimol / liter
    
    'FFA_prod_nId'      'p4'                                                             {} % micromol / liter   
    'FFA_prod_Id'       '3 * p8/(k3+insulinAct2)'                                                             {} % micromol / liter   
%     'FFA_upt_Id'        'p6*FFA*insulinAct2'                                                             {} % micromol / liter   
%     'FFA_upt_nId'       'k2*FFA'    {}
    'FFA_upt_Id'       'k2*FFA*insulinAct2'    {}
    
    
    'gly_prod_nId'      'k4'                                                             {} % micromol / liter
    'gly_prod_Id'       'p8/(k3+insulinAct2)'                                                             {} % micromol / liter
    'gly_upt_nId'       'p10*glycerol'                                                       {} % micromol / liter

    'gly_inf_tracer'    'BW/V_Y * (gly1*TTR_fix + gly2*TTR_fix)'                                        {} % micromol / liter
    'gly_inf_tracee'    'BW/V_Y * (gly1*(1-TTR_fix) + gly2*(1-TTR_fix))'                                {} % micromol / liter
    
    'I_prod'            'p12'                                                                   {} % picomol / liter
    'I_inf'             '1/V_I * ins1'                                                                  {} % picomol / liter
    'I_deg'             'p13*insulin'                                                                   {} % picomol / liter
    
    'glu_RA'            'glu_prod_nId'  {}
    'glu_RD'            'glu_upt_Id + glu_upt_nId'  {}
 
    'gly_RA'            'gly_prod_nId + gly_prod_Id'  {}
    'gly_RD'            'gly_upt_nId'  {}
    
    'IS'                'glu_RD / insulin'  {}
    
    };