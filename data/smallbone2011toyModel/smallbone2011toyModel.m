function MODEL = smallbone2011toyModel(~)

MODEL.DESCRIPTION = 'smallbone2011 Toy model created by D. Lao-Martil';

MODEL.PREDICTOR = {
    't'  [0 10] {'time' 'minutes' 'Time'}
};

MODEL.CONSTANTS = {
    'u1' 1 [] {}
    'u2' 1 [] {}
};

MODEL.PARAMETERS = {
    'p_extra_cell' 1  [] {}
    'p_extra_medium' 1  [] {}
    % kinetic paramters: G6P isomerase (PGI)
    'p_PGI_Vmax' 1071   [] {}
    'p_PGI_Kg6p' 1.4   [] {} %Mm
    'p_PGI_Kf6p' 0.29   [] {} %Mm
    'p_PGI_Keq' 0.3   [] {}
    'p_PGI_shock' 1   [] {}
    'p_PGI_heat' 0   [] {}
    % kinetic paramters: glucose transport (GLT)
    'p_GLT_Vmax' 97.24   [] {}
    'p_GLT_Kglc' 1.1918   [] {} %mM
    'p_GLT_Ki' 0.91   [] {}
    'p_GLT_shock' 8   [] {}
    'p_GLT_heat' 0   [] {}
    % kinetic paramters: hexokinase (HK)
    'p_HXK_Vmax' 289.6   [] {}
    'p_HXK_Kglc' 0.08   [] {} %mM
    'p_HXK_Katp' 0.15   [] {} %mM
    'p_HXK_Kg6p' 30   [] {} %mM
    'p_HXK_Kadp' 0.23   [] {} %mM
    'p_HXK_Keq' 2000   [] {}
    'p_HXK_Kit6p' 0.04   [] {} %mM
    'p_HXK_shock' 8   [] {}
    'p_HXK_heat' 0   [] {}
    % kinetic paramters: phosphoglucomutase (PGM1)
    'p_PGM1_Vmax' 0.3545   [] {}
    'p_PGM1_Kg6p' 0.05   [] {} %mM
    'p_PGM1_Kg1p' 0.023   [] {} %mM
    'p_PGM1_Keq' 0.1667   [] {}
    'p_PGM1_shock' 16   [] {}
    'p_PGM1_heat' 0   [] {}
    % kinetic paramters: T6P phosphatase (TPS2)
    'p_TPS2_Vmax' 6.5   [] {}
    'p_TPS2_Kt6p' 0.5   [] {} %mM
    'p_TPS2_shock' 18   [] {}
    'p_TPS2_heat' 0   [] {}
    % kinetic paramters: T6p synthase (TPS1)
    'p_TPS1_Vmax' 1.371   [] {}
    'p_TPS1_Kg6p' 3.8   [] {} %mM
    'p_TPS1_Kudg' 0.886   [] {} %mM
    'p_TPS1_shock' 12   [] {}
    'p_TPS1_activity' 1   [] {}
    'p_TPS1_heat' 0   [] {}
    % kinetic paramters: trehalase (NTH1)
    'p_NTH1_Vmax' 15.2   [] {}
    'p_NTH1_Ktrh' 2.99   [] {} %mM
    'p_NTH1_shock' 6   [] {}
    'p_NTH1_heat' 0   [] {}
    % kinetic paramters: UDP glucose phosphorylase (UPG)
    'p_UPG_Vmax' 36.82   [] {}
    'p_UPG_Kutp' 0.11   [] {}
    'p_UPG_Kiutp' 0.11   [] {}
    'p_UPG_Kg1p' 0.32   [] {}
    'p_UPG_Kiudg' 0.0035   [] {}
    'p_UPG_shock' 16   [] {}
    'p_UPG_heat' 0   [] {}
    % fixed metabolite concentrations
    'ADP'      1.282   [] {}
    'ATP'      2.525   [] {}
    'diphosphate' 1   [] {}
    'F6P'      0.625   [] {}
    'Hplus'        1   [] {}
    'phosphate'    1   [] {}
%     'UDP'      0.2815   [] {}
    'UTP'      0.6491   [] {}
    'water'         1   [] {}
    'glucose_ext' 100   [] {} % main metabolite to change
};

MODEL.STATES = {
    's1' 0.09765    's1' 'v1 - v2 + 2 * v7'  {} %xdot(1) = + v_GLT - v_HXK + 2 * v_NTH1;        %Glc
    's2' 0.1        's2' 'v3 - v4'           {} %xdot(2) = + v_PGM1 - v_UDP;                     %G1P
    's3' 2.675      's3' 'v2 - v3 - v5 - v8' {} %xdot(3) = + v_HXK - v_PGM1 - v_TPS1 - v_PGI;    %G6P
    's4' 0.05       's4' 'v6 - v7'           {} %xdot(4) = + v_TPS2 - v_NTH1;                    %Trehalose
    's5' 0.02       's5' 'v5 - v6'           {} %xdot(5) = + v_TPS1 - v_TPS2;                    %T6P
    's6' 0.7        's6' 'v4 - v5'           {} %xdot(6) = + v_UDP - v_TPS1;                     %UDP_glucose
};

MODEL.REACTIONS = {
    'v1' 'p_GLT_Vmax * (glucose_ext - s1) / p_GLT_Kglc / (1 + (glucose_ext + s1) / p_GLT_Kglc + (p_GLT_Ki * s1 * glucose_ext) / (p_GLT_Kglc * p_GLT_Kglc))' {}
    'v2' 'p_HXK_Vmax * 1 / (p_HXK_Kglc * p_HXK_Katp) * (s1 * ATP - s3 * ADP / p_HXK_Keq) / (1 + s1 / p_HXK_Kglc + s3 / p_HXK_Kg6p + s5 / p_HXK_Kit6p) / (1 + ATP / p_HXK_Katp + ADP / p_HXK_Kadp)' {}
    'v3' 'p_PGM1_Vmax * 1 / p_PGM1_Kg6p * (s3 - s2 / p_PGM1_Keq) / (1 + s3 / p_PGM1_Kg6p + s2 / p_PGM1_Kg1p)'      {}
    'v4' 'p_UPG_Vmax * UTP * s2 / (p_UPG_Kutp * p_UPG_Kg1p) / (p_UPG_Kiutp / p_UPG_Kutp + UTP / p_UPG_Kutp + s2 / p_UPG_Kg1p + UTP* s2 / (p_UPG_Kutp * p_UPG_Kg1p) + p_UPG_Kiutp / p_UPG_Kutp * s6 / p_UPG_Kiudg + s2 / p_UPG_Kg1p * s6 / p_UPG_Kiudg)'      {}
    'v5' 'p_TPS1_Vmax * s3 / p_TPS1_Kg6p * s6 / p_TPS1_Kudg / (1 + s3 / p_TPS1_Kg6p) / (1 + s6 / p_TPS1_Kudg)'      {}
    'v6' 'p_TPS2_Vmax * s5 / p_TPS2_Kt6p / (1 + s5 / p_TPS2_Kt6p)'      {}
    'v7' 'p_NTH1_Vmax * s4 / p_NTH1_Ktrh / (1 + s4 / p_NTH1_Ktrh)'      {}
    'v8' 'p_PGI_Vmax * 1 / p_PGI_Kg6p * (s3 - F6P / p_PGI_Keq) / (1 + s3 / p_PGI_Kg6p + F6P / p_PGI_Kf6p)'      {}
%     'v1' 'p_GLT_Vmax * (glucose_ext - Glc) / p_GLT_Kglc / (1 + (glucose_ext + Glc) / p_GLT_Kglc + (p_GLT_Ki * Glc * glucose_ext) / (p_GLT_Kglc * p_GLT_Kglc))' {}
%     'v2' 'p_HXK_Vmax * 1 / (p_HXK_Kglc * p_HXK_Katp) * (Glc * ATP - G6P * ADP / p_HXK_Keq) / (1 + Glc / p_HXK_Kglc + G6P / p_HXK_Kg6p + T6P / p_HXK_Kit6p) / (1 + ATP / p_HXK_Katp + ADP / p_HXK_Kadp)' {}
%     'v3' 'p_PGM1_Vmax * 1 / p_PGM1_Kg6p * (G6P - G1P / p_PGM1_Keq) / (1 + G6P / p_PGM1_Kg6p + G1P / p_PGM1_Kg1p)'      {}
%     'v4' 'p_UPG_Vmax * UTP * G1P / (p_UPG_Kutp * p_UPG_Kg1p) / (p_UPG_Kiutp / p_UPG_Kutp + UTP / p_UPG_Kutp + G1P / p_UPG_Kg1p + UTP* G1P / (p_UPG_Kutp * p_UPG_Kg1p) + p_UPG_Kiutp / p_UPG_Kutp * UDP_glucose / p_UPG_Kiudg + G1P / p_UPG_Kg1p * UDP_glucose / p_UPG_Kiudg)'      {}
%     'v5' 'p_TPS1_Vmax * G6P / p_TPS1_Kg6p * UDP_glucose / p_TPS1_Kudg / (1 + G6P / p_TPS1_Kg6p) / (1 + UDP_glucose / p_TPS1_Kudg)'      {}
%     'v6' 'p_TPS2_Vmax * T6P / p_TPS2_Kt6p / (1 + T6P / p_TPS2_Kt6p)'      {}
%     'v7' 'p_NTH1_Vmax * Trehalose / p_NTH1_Ktrh / (1 + Trehalose / p_NTH1_Ktrh)'      {}
%     'v8' 'p_PGI_Vmax * 1 / p_PGI_Kg6p * (G6P - F6P / p_PGI_Keq) / (1 + G6P / p_PGI_Kg6p + F6P / p_PGI_Kf6p)'      {}
};