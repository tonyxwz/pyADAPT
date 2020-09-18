function [v] = calculateReactionRates_smallbone2011(x,p)

% recall metabolite values
Glc         = x(:,1);
G1P         = x(:,2);
G6P         = x(:,3);
Trehalose   = x(:,4);
T6P         = x(:,5);
UDP_glucose = x(:,6);

% calculate reaction rates
rateEquations_smallbone2011;
v = [v_GLT, v_HXK, v_PGM1, v_UDP, v_TPS1, v_TPS2, v_NTH1, v_PGI];