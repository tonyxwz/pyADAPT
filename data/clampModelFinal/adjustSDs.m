data = load('DataRaw.mat');

names = fieldnames(data);

% Variable      Statistical Error	Systemic Error	Squared total
% Glucose       From data           2%              sqrt(f.d.^2 + 2^2)
% Glycerol      3%                  3%              sqrt(3^2 + 3^2)
% Insulin       5%                  5%              sqrt(5^2 + 5^2)
% FFA           1%                  9.5%            sqrt(1^2 + 9.5^2)
% Glucose TTR	5%                  5%              sqrt(5^2 + 5^2)
% Glycerol TTR	1.5%                6%              sqrt(1.5^2 + 6^2)

CV_glu = 0.01*sqrt(2^2);
CV_glu_ttr = 0.01*sqrt(5^2 + 5^2);
CV_gly = 0.01*sqrt(3^2 + 3^2);
CV_gly_ttr = 0.01*sqrt(1.5^2 + 6^2);
CV_ins = 0.01*sqrt(5^2 + 5^2);
CV_ffa = 0.01*sqrt(1^2 + 9.5^2);

for i = 1:length(names)
    
    glu_raw = data.(names{i}).glu_raw;
    glu_raw_std = nanstd(glu_raw, [], 2);
    glu_std = glu_raw_std(~isnan(glu_raw_std));
    
    CV_glu = 0.01*sqrt(glu_std.^2 + 2^2);
    
    data.(names{i}).glu_std = data.(names{i}).glu_conc .* CV_glu;
    data.(names{i}).glu_ttr_std = data.(names{i}).glu_ttr * CV_glu_ttr;
    data.(names{i}).gly_std = data.(names{i}).gly_conc * CV_gly;
    data.(names{i}).gly_ttr_std = data.(names{i}).gly_ttr * CV_gly_ttr;
    data.(names{i}).ins_std = data.(names{i}).ins_conc * CV_ins;
    data.(names{i}).FFA_std = data.(names{i}).FFA_conc * CV_ffa;
    
    data.(names{i}).glu2_time = [data.(names{i}).glu2_time(1); 1e-9; data.(names{i}).glu2_time(2:end)];
    data.(names{i}).glu2_flux = [0; data.(names{i}).glu2_flux];
    data.(names{i}).glu2_std  = [0; data.(names{i}).glu2_std];

    data.(names{i}).glu3_time2 = [data.(names{i}).glu3_time2(1); 1e-9; data.(names{i}).glu3_time2(2:end)];
    data.(names{i}).glu3_flux2 = [0; data.(names{i}).glu3_flux2];
    data.(names{i}).glu3_std2  = [0; data.(names{i}).glu3_std2];
    
    data.(names{i}).gly1_time = [data.(names{i}).gly1_time(1); 1e-9; data.(names{i}).gly1_time(2:end)];
    data.(names{i}).gly1_flux = [0; data.(names{i}).gly1_flux];
    data.(names{i}).gly1_std  = [0; data.(names{i}).gly1_std];
    
    data.(names{i}).gly2_time2 = [data.(names{i}).gly2_time2(1); 1e-9; data.(names{i}).gly2_time2(2:end)];
    data.(names{i}).gly2_flux2 = [0; data.(names{i}).gly2_flux2];
    data.(names{i}).gly2_std2  = [0; data.(names{i}).gly2_std2];
end

save('DataFinal2.mat','-struct','data');