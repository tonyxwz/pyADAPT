% clc, clear, close all
dbstop if error

%% Simulation protocol

% First simulation: ss values to then start the pulse
setParameterStructure_smallbone2011;
p.fxdSps.glucose_ext    = 0.1;
t_end = 50; % time in minutes
[tss,xss] = ode15s(@Smallbone2011, [0, t_end], Smallbone2011,[],p);
[vss] = calculateReactionRates_smallbone2011(xss,p);

% pulse simulation
p.fxdSps.glucose_ext    = 10;
t_end = 50; % time in minutes
[tgp,xgp] = ode15s(@Smallbone2011, [0, t_end], xss(end,:),[],p);
[vgp] = calculateReactionRates_smallbone2011(xgp,p);

tprev = (tss - tss(end));
t = [tprev(1:end-1); tgp];
x = [xss(1:end-1,:); xgp];
v = [vss(1:end-1,:); vgp];


%% adding white noise (normally distributed to the data)

[xm,xn] = size(x);
[vm,vn] = size(v);
nsims   = 5; % the sample size is selected low enough so that the mean value is still a bit off (by trial and error)

x_data = cell(nsims,1);
v_data = cell(nsims,1);

for i = 1:nsims
    % create randomized simulations for metabolite concentrations
    x_data{i} = x;
    for j = 1:xm
        for k = 1:xn
            x_data{i}(j,k) = x(j,k) + x(j,k)*randn*0.01;
        end
    end
    % create randomized simulations for reaction rates
    v_data{i} = v;
    for j = 1:vm
        for k = 1:vn
            v_data{i}(j,k) = v(j,k) + v(j,k)*randn*0.01;
        end
    end
end


%% selecting standard deviation
% we assume few datapoints. focused in the regions where the dataset deviates
% the most.
t_exp           = [0    0.0002 	0.1     0.5     1   10  40]; 
% id for times
% t_ids           = [86   91      131     140     145 157 165]; % 5 data poins
t_ids           = [86   140 145 150 157 161 165]; 

Glc_exp_array           = zeros(length(t_ids),nsims); 
G1P_exp_array           = Glc_exp_array;
G6P_exp_array           = Glc_exp_array;
trehalose_exp_array     = Glc_exp_array;
T6P_exp_array           = Glc_exp_array;
UDPglucose_exp_array    = Glc_exp_array;
 for i = 1:nsims
     for j = 1:length(t_ids)
        Glc_exp_array(j,i)         = x_data{i}(t_ids(j),1); 
        G1P_exp_array(j,i)         = x_data{i}(t_ids(j),2); 
        G6P_exp_array(j,i)         = x_data{i}(t_ids(j),3); 
        trehalose_exp_array(j,i)   = x_data{i}(t_ids(j),4); 
        T6P_exp_array(j,i)         = x_data{i}(t_ids(j),5); 
        UDPglucose_exp_array(j,i)  = x_data{i}(t_ids(j),6);
%         disp(Glc_exp_array);
     end
end
Glc_exp         = mean(Glc_exp_array,2);
G1P_exp         = mean(G1P_exp_array,2);
G6P_exp         = mean(G6P_exp_array,2);
trehalose_exp   = mean(trehalose_exp_array,2); 
T6P_exp         = mean(T6P_exp_array,2); 
UDPglucose_exp  = mean(UDPglucose_exp_array,2); 


std_Glc_exp         = std(Glc_exp_array,0,2);
std_G1P_exp         = std(G1P_exp_array,0,2);
std_G6P_exp         = std(G6P_exp_array,0,2);
std_trehalose_exp   = std(trehalose_exp_array,0,2);
std_T6P_exp         = std(T6P_exp_array,0,2);
std_UDPglucose_exp  = std(UDPglucose_exp_array,0,2);
% std_Glc_exp         = Glc_exp.*randn        .*0.01; %0.02;
% std_G1P_exp         = G1P_exp.*randn        .*0.01; %0.02;
% std_G6P_exp         = G6P_exp.*randn        .*0.01; %0.02;
% std_trehalose_exp   = trehalose_exp.*randn  .*0.01; %0.02;
% std_T6P_exp         = T6P_exp.*randn        .*0.01; %0.05;
% std_UDPglucose_exp  = UDPglucose_exp.*randn .*0.01; %0.02;


% in structure (easier to recall later)
x_exp.mean.Glc_exp = Glc_exp;
x_exp.mean.G1P_exp = G1P_exp;
x_exp.mean.G6P_exp = G6P_exp;
x_exp.mean.trehalose_exp = trehalose_exp;
x_exp.mean.T6P_exp = T6P_exp;
x_exp.mean.UDPglucose_exp = UDPglucose_exp;
x_exp.std.std_Glc_exp = std_Glc_exp;
x_exp.std.std_G1P_exp = std_G1P_exp;
x_exp.std.std_G6p_exp = std_G6P_exp;
x_exp.std.std_trehalose_exp = std_trehalose_exp;
x_exp.std.std_T6P_exp = std_T6P_exp;
x_exp.std.std_UDPglucose_exp = std_UDPglucose_exp;
 

v_GLT_array    = zeros(length(t_ids),nsims); 
v_HXK_array    = v_GLT_array; 
v_PGM1_array   = v_GLT_array;
v_UDP_array    = v_GLT_array;
v_TPS1_array   = v_GLT_array;
v_TPS2_array   = v_GLT_array;
v_NTH1_array   = v_GLT_array;
v_PGI_array    = v_GLT_array;
 for i = 1:nsims
     for j = 1:length(t_ids)
        v_GLT_array(j,i)           = v_data{i}(t_ids(j),1); 
        v_HXK_array(j,i)           = v_data{i}(t_ids(j),2); 
        v_PGM1_array(j,i)          = v_data{i}(t_ids(j),3); 
        v_UDP_array(j,i)           = v_data{i}(t_ids(j),4); 
        v_TPS1_array(j,i)          = v_data{i}(t_ids(j),5); 
        v_TPS2_array(j,i)          = v_data{i}(t_ids(j),6); 
        v_NTH1_array(j,i)          = v_data{i}(t_ids(j),7); 
        v_PGI_array(j,i)           = v_data{i}(t_ids(j),8); 
     end
end
v_GLT           = mean(v_GLT_array,2);
v_HXK           = mean(v_HXK_array,2);
v_PGM1          = mean(v_PGM1_array,2);
v_UDP           = mean(v_UDP_array,2);
v_TPS1          = mean(v_TPS1_array,2);
v_TPS2          = mean(v_TPS2_array,2);
v_NTH1          = mean(v_NTH1_array,2);
v_PGI           = mean(v_PGI_array,2);

std_v_GLT         = std(v_GLT_array,0,2);
std_v_HXK         = std(v_HXK_array,0,2);
std_v_PGM1        = std(v_PGM1_array,0,2);
std_v_UDP         = std(v_UDP_array,0,2);
std_v_TPS1        = std(v_TPS1_array,0,2);
std_v_TPS2        = std(v_TPS2_array,0,2);
std_v_NTH1        = std(v_NTH1_array,0,2); 
std_v_PGI         = std(v_PGI_array,0,2);
% std_v_GLT           = v_GLT.*randn  .*0.01; 
% std_v_HXK           = v_HXK.*randn  .*0.01; 
% std_v_PGM1          = v_PGM1.*randn .*0.01; 
% std_v_UDP           = v_UDP.*randn  .*0.01; 
% std_v_TPS1          = v_TPS1.*randn .*0.01; 
% std_v_TPS2          = v_TPS2.*randn .*0.01; 
% std_v_NTH1          = v_NTH1.*randn .*0.01; 
% std_v_PGI           = v_PGI.*randn  .*0.01; 

% in structure (easier to recall later)
v_exp.mean.v_GLT = v_GLT;
v_exp.mean.v_HXK = v_HXK;
v_exp.mean.v_PGM1 = v_PGM1;
v_exp.mean.v_UDP = v_UDP;
v_exp.mean.v_TPS1 = v_TPS1;
v_exp.mean.v_TPS2 = v_TPS2;
v_exp.mean.v_NTH1 = v_NTH1;
v_exp.mean.v_PGI = v_PGI;
v_exp.std.std_v_GLT = std_v_GLT;
v_exp.std.std_v_HXK = std_v_HXK;
v_exp.std.std_v_PGM1 = std_v_PGM1;
v_exp.std.std_v_UDP = std_v_UDP;
v_exp.std.std_v_TPS1 = std_v_TPS1;
v_exp.std.std_v_TPS2 = std_v_TPS2;
v_exp.std.std_v_NTH1 = std_v_NTH1;
v_exp.std.std_v_PGI = std_v_PGI;


A1 = x_exp.mean;    B1 = (struct2cell(A1));
A2 = x_exp.std;     B2 = (struct2cell(A2));
A3 = v_exp.mean;    B3 = (struct2cell(A3));
A4 = v_exp.std;     B4 = (struct2cell(A4));


%% visualization

legendaMets     = {'glucose_{in}','G1P','G6P','trehalose','T6P','UDPglucose'};
legendaFluxes   = {'v_{GLT}', 'v_{HXK}', 'v_{PGM1}', 'v_{UDP}', 'v_{TPS1}', 'v_{TPS2}', 'v_{NTH1}', 'v_{PGI}'};

% figure(1)
% plot(t,x)
% legend(legendaMets)
% xlim([-10 t_end])
% title('metabolites vs time (all together)')

% figure(2)
% for i = 1:6
%     subplot(2,3,i)
%     plot(t,x(:,i))
%     title(legendaMets{i})
%     xlim([-10 t_end])
% end
% suptitle('metabolites vs time (one per plot)')
figure(2)
for i = 1:6
    subplot(2,3,i)
    for j = 1:nsims
        plot(t,x_data{j}(:,i), 'Color', [0.7 0.7 0.7])
        hold on        
    end
    hold on
    plot(t,x(:,i), 'k', 'LineWidth', 2)
    % exp data
    hold on
    errorbar(t(t_ids),B1{i},B2{i},B2{i},'o','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','red','color','black')
    title(legendaMets{i})
    xlim([-10 t_end])    
end
suptitle('metabolites vs time (one per plot)')

% figure(3)
% plot(t,v)
% legend(legendaFluxes)
% xlim([-10 t_end])
% title('fluxes vs time (all together)')

% figure(4)
% for i = 1:8
%     subplot(3,3,i)
%     plot(t,v(:,i))
%     title(legendaFluxes{i})
%     xlim([-10 t_end])
% end
% suptitle('fluxes vs time (one per plot)')
figure(4)
for i = 1:8
    subplot(3,3,i)
    for j = 1:nsims
        plot(t,v_data{j}(:,i), 'Color', [0.7 0.7 0.7])
        hold on        
    end
    hold on
    plot(t,v(:,i), 'k', 'LineWidth', 2)
    hold on
    errorbar(t(t_ids),B3{i},B4{i},B4{i},'o','MarkerSize',4,'MarkerEdgeColor','black','MarkerFaceColor','red','color','black')
    title(legendaFluxes{i})
    xlim([-10 t_end])    
end
suptitle('fluxes vs time (one per plot)')

%% saving data and figures
savefig(2,'Smallbone2011_toy_metabolites.fig')
savefig(4,'Smallbone2011_toy_fluxes.fig')
save('smallbone2011_data','x_exp','v_exp')


