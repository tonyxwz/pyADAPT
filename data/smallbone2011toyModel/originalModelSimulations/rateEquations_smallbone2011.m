% reassigning fixed metabolite concentrations
ADP = p.fxdSps.ADP;
ATP = p.fxdSps.ATP;
diphosphate = p.fxdSps.diphosphate;
F6P = p.fxdSps.F6P;
Hplus = p.fxdSps.Hplus;
phosphate = p.fxdSps.phosphate;
UDP = p.fxdSps.UDP;
UTP = p.fxdSps.UTP;
water = p.fxdSps.water;
glucose_ext = p.fxdSps.glucose_ext;

% glucose transport (GLT)
v_GLT   = p.GLT.Vmax .* (glucose_ext - Glc) ./ p.GLT.Kglc ./ (1 + (glucose_ext + Glc) ./ p.GLT.Kglc + (p.GLT.Ki .* Glc .* glucose_ext) ./ (p.GLT.Kglc .* p.GLT.Kglc));

% hexokinase (HK)
v_HXK   = p.HXK.Vmax .* 1 ./ (p.HXK.Kglc .* p.HXK.Katp) .* (Glc .* ATP - G6P .* ADP ./ p.HXK.Keq) ./ (1 + Glc ./ p.HXK.Kglc + G6P ./ p.HXK.Kg6p + T6P ./ p.HXK.Kit6p) ./ (1 + ATP ./ p.HXK.Katp + ADP ./ p.HXK.Kadp);

% phosphoglucomutase (PGM1)
v_PGM1  = p.PGM1.Vmax .* 1 ./ p.PGM1.Kg6p .* (G6P - G1P ./ p.PGM1.Keq) ./ (1 + G6P ./ p.PGM1.Kg6p + G1P ./ p.PGM1.Kg1p);

% UDP glucose phosphorylase (UPG)
v_UDP   = p.UPG.Vmax .* UTP .* G1P ./ (p.UPG.Kutp .* p.UPG.Kg1p) ./ (p.UPG.Kiutp ./ p.UPG.Kutp + UTP ./ p.UPG.Kutp + G1P ./ p.UPG.Kg1p + UTP.* G1P ./ (p.UPG.Kutp .* p.UPG.Kg1p) + p.UPG.Kiutp ./ p.UPG.Kutp .* UDP_glucose ./ p.UPG.Kiudg + G1P ./ p.UPG.Kg1p .* UDP_glucose ./ p.UPG.Kiudg);

% T6P synthase (TPS1)
v_TPS1  = p.TPS1.Vmax .* G6P ./ p.TPS1.Kg6p .* UDP_glucose ./ p.TPS1.Kudg ./ (1 + G6P ./ p.TPS1.Kg6p) ./ (1 + UDP_glucose ./ p.TPS1.Kudg);

% T6P phosphatase (TPS2)
v_TPS2  = p.TPS2.Vmax .* T6P ./ p.TPS2.Kt6p ./ (1 + T6P ./ p.TPS2.Kt6p);

% Trehalase (NTH1)
v_NTH1  = p.NTH1.Vmax .* Trehalose ./ p.NTH1.Ktrh ./ (1 + Trehalose ./ p.NTH1.Ktrh);

% G6P isomerase (PGI)
v_PGI   = p.PGI.Vmax .* 1 ./ p.PGI.Kg6p .* (G6P - F6P ./ p.PGI.Keq) ./ (1 + G6P ./ p.PGI.Kg6p + F6P ./ p.PGI.Kf6p);

