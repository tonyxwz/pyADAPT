% initial compartment sizes
p.extra.cell        = 1; %L
p.extra.medium      = 1; %L

% initial global quantities
p.extra.heat        = 0;
p.extra.glc_0       = 0.09765;
p.extra.log10_glc_0 = 0;
p.extra.g1p_0       = 0.1;
p.extra.log10_g1p_0 = 0;
p.extra.g6p_0       = 2.675;
p.extra.log10_g6p_0 = 0;
p.extra.trh_0       = 0.05;
p.extra.log10_trh_0 = 0;
p.extra.t6p_0       = 0.02;
p.extra.log10_t6p_0 = 0;
p.extra.udg_0       = 0.7;
p.extra.log10_udg_0 = 0;

% kinetic paramters: G6P isomerase (PGI)
p.PGI.Vmax      = 1071;
p.PGI.Kg6p      = 1.4; %Mm
p.PGI.Kf6p      = 0.29; %Mm
p.PGI.Keq       = 0.3;
p.PGI.shock     = 1;
p.PGI.heat      = 0;

% kinetic paramters: glucose transport (GLT)
p.GLT.Vmax      = 97.24;
p.GLT.Kglc      = 1.1918; %mM
p.GLT.Ki        = 0.91;
p.GLT.shock     = 8;
p.GLT.heat      = 0;

% kinetic paramters: hexokinase (HK)
p.HXK.Vmax      = 289.6;
p.HXK.Kglc      = 0.08; %mM
p.HXK.Katp      = 0.15; %mM
p.HXK.Kg6p      = 30; %mM
p.HXK.Kadp      = 0.23; %mM
p.HXK.Keq       = 2000;
p.HXK.Kit6p     = 0.04; %mM
p.HXK.shock     = 8;
p.HXK.heat      = 0;

% kinetic paramters: phosphoglucomutase (PGM1)
p.PGM1.Vmax     = 0.3545;
p.PGM1.Kg6p     = 0.05; %mM
p.PGM1.Kg1p     = 0.023; %mM
p.PGM1.Keq      = 0.1667;
p.PGM1.shock    = 16;
p.PGM1.heat     = 0;

% kinetic paramters: T6P phosphatase (TPS2)
p.TPS2.Vmax     = 6.5;
p.TPS2.Kt6p     = 0.5; %mM
p.TPS2.shock    = 18;
p.TPS2.heat     = 0;

% kinetic paramters: T6p synthase (TPS1)
p.TPS1.Vmax     = 1.371;
p.TPS1.Kg6p     = 3.8; %mM
p.TPS1.Kudg     = 0.886; %mM
p.TPS1.shock    = 12;
p.TPS1.activity = 1;
p.TPS1.heat     = 0;

% kinetic paramters: trehalase (NTH1)
p.NTH1.Vmax     = 15.2;
p.NTH1.Ktrh     = 2.99; %mM
p.NTH1.shock    = 6;
p.NTH1.heat   	= 0;

% kinetic paramters: UDP glucose phosphorylase (UPG)
p.UPG.Vmax      = 36.82;
p.UPG.Kutp      = 0.11;
p.UPG.Kiutp     = 0.11;
p.UPG.Kg1p      = 0.32;
p.UPG.Kiudg     = 0.0035;
p.UPG.shock     = 16;
p.UPG.heat      = 0;

% fixed metabolite concentrations
p.fxdSps.ADP            = 1.282;
p.fxdSps.ATP            = 2.525;
p.fxdSps.diphosphate    = 1;
p.fxdSps.F6P            = 0.625;
p.fxdSps.Hplus          = 1;
p.fxdSps.phosphate      = 1;
p.fxdSps.UDP            = 0.2815;
p.fxdSps.UTP            = 0.6491;
p.fxdSps.water          = 1;
p.fxdSps.glucose_ext    = 100; % main metabolite to change