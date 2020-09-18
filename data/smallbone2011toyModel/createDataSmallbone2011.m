load('draft2_smallbone2011toyData.mat');
load('draft_smallbone2011_data.mat');

smallbone2011toy.t          = [0 0.5 1 2.3 10 20 40]';
smallbone2011toy.s1_mean    = x_exp.mean.Glc_exp;
smallbone2011toy.s1_std     = x_exp.std.std_Glc_exp;
smallbone2011toy.s2_mean    = x_exp.mean.G1P_exp;
smallbone2011toy.s2_std     = x_exp.std.std_G1P_exp;
smallbone2011toy.s3_mean    = x_exp.mean.G6P_exp;
smallbone2011toy.s3_std     = x_exp.std.std_G6p_exp;
smallbone2011toy.s4_mean    = x_exp.mean.trehalose_exp;
smallbone2011toy.s4_std     = x_exp.std.std_trehalose_exp;
smallbone2011toy.s5_mean    = x_exp.mean.T6P_exp;
smallbone2011toy.s5_std     = x_exp.std.std_T6P_exp;
smallbone2011toy.s6_mean    = x_exp.mean.UDPglucose_exp;
smallbone2011toy.s6_std     = x_exp.std.std_UDPglucose_exp;
    
saveName = 'smallbone2011_data.mat';
save(saveName,'smallbone2011toy');