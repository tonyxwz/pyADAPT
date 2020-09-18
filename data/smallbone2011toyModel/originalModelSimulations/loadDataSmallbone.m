%% Load Data
load('vHeerden_trehalose_data_micromolgdw.mat'); 
    data.nucleotides=data.nucleotides/2;
    data.metabolites=data.metabolites/2;
    data.fluxes=data.fluxes/2;
        AXP=sum(data.nucleotides(:,2:4)'); 
        time=[0:1:340];
        AXP_intp=interp1(data.time_nucleotides,AXP,time,'linear','extrap');
        dAXPdt=diff(AXP_intp);
    data.time_dAXPdt=time(2:end);
    data.AXP=AXP;
    data.dAXPdt=dAXPdt;
    clear AXP time AXP_intp dAXPdt
    
%% plot this data

% % % % metabolite profile vs time
figure(3)
% figure(103)

% G6P
subplot(2,3,1)
plot(data.time_metabolites,data.metabolites(:,2),'r+')
title('G6P')
xlabel('time [s]')
ylabel('concentration [mM]')
xlim([-100 340])

% F6P
subplot(2,3,2)
plot(data.time_metabolites,data.metabolites(:,3),'r+')
title('F6P')
% xlabel('time [s]')
% ylabel('concentration [mM]')
xlim([-100 340])

% G1P
subplot(2,3,3)
plot(data.time_metabolites,data.metabolites(:,13),'r+')
title('G1P')
% xlabel('time [s]')
% ylabel('concentration [mM]')
xlim([-100 340])

% T6P
subplot(2,3,4)
plot(data.time_metabolites,data.metabolites(:,15),'r+')
title('T6P')
% xlabel('time [s]')
% ylabel('concentration [mM]')
xlim([-100 340])

% TRE
subplot(2,3,5)
plot(data.time_metabolites, data.metabolites(:,16),'r+')
title('TRE')
xlim([-100 340])

% UDP_glc
subplot(2,3,6)
plot(data.time_metabolites, data.metabolites(:,14),'r+')
title('UDP_glc')
xlim([-100 340])

suptitle({'Metabolite concentrations vs time';'Glucose pulse from 0.1 to 110 mM';'starting conentrations follow from steadt state'})

%%

figure(4)
% figure(104)

% v_GLT
subplot(2,4,1)
plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
title('GLT')
xlabel('time [s^{-1}]')
ylabel('reaction rate [mM s^{-1}]')
xlim([-100 340])

% v_HK
subplot(2,4,2)
plot(data.time_fluxes(1:end),data.fluxes(:,1), 'r+')
title('HK')
xlim([-100 340])

% v_PGI
subplot(2,4,3)
plot(data.time_fluxes(1:end),data.fluxes(:,2)-data.fluxes(:,3), 'r+')
title('PGI')
xlim([-100 340])

% v_PFK
subplot(2,4,4)
plot(data.time_fluxes(1:end),data.fluxes(:,4), 'r+')
title('PFK')
xlim([-100 340])

% v_PGM1
subplot(2,4,5)
plot(data.time_fluxes(:), data.fluxes(:,6)-data.fluxes(:,7), 'r+')
title('PGM1')
xlim([-100 340])

% v_TPS1
subplot(2,4,6)
plot(data.time_fluxes(:), data.fluxes(:,9), 'r+')
title('TPS1')
xlim([-100 340])

% v_TPS2
subplot(2,4,7)
plot(data.time_fluxes(:), data.fluxes(:,10), 'r+')
title('TPS2')
xlim([-100 340])

% v_NTH1
subplot(2,4,8)
plot(data.time_fluxes(:), data.fluxes(:,11), 'r+')
title('NTH1')
xlim([-100 340])

suptitle({'van Heerden GS dataset. Fluxes profile at changing times';'The first measurement for v_{GLT}, v_{HK}, v_{PGI} and v_{PFK} comes after the pulse already '})


