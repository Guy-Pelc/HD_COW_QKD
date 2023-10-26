% generate combined results file from all computed data

clear r;
h = helper;
BLUE = "#0072BD";
RED = "#D95319";
ORANGE = "#EDB120";
%%
x = load('res_m1D2_0.50_err_2023.02.21.1447.mat'); r.m1d2_050 = x.res;
x = load('res_m1D3_0.50_err_2023.02.21.1448.mat'); r.m1d3_050 = x.res;
x= load('res_m1D4_0.50_err_2023.02.21.1450.mat'); r.m1d4_050 = x.res;



%% Rate Per Block
figure(5); clf; 
% subplot(2,1,1); cla
hold on; grid on; set(gca,'YScale','log');
% ylim([1e3,2e5])
% title('BER = 0.5%, Visibility = 99.75%')
leg = [];
x = r.m1d2_050; res = h.rpb(x); plot(x.lossDBs,res,'-','color',BLUE,LineWidth=2); leg = [leg,"d=2"];
x = r.m1d3_050; res = h.rpb(x); plot(x.lossDBs,res,'-','color',RED,LineWidth=2); leg = [leg,"d=3"];
x = r.m1d4_050; res = h.rpb(x); plot(x.lossDBs,res,'-','color',ORANGE,LineWidth=2); leg = [leg,"d=4"];

legend(leg)
ylabel('Secure bits per block')
xlabel('Loss (dB)')
xlim([1,12])

fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 6;
set(fig.Children,     'FontSize',     9);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig.PaperPositionMode   = 'auto';
box on
%%  Rate Per Second
figure(6); clf; 

hold on; grid on; set(gca,'YScale','log');
% ylim([5e3,8e4])
xlim([1,12])
leg = [];

% test for different saturations
T = 0;
title("deadtime="+T)
% lower bound
x = r.m1d2_050; res = h.rps(x,'ps',1,'T',T); plot(x.lossDBs,res,'-','color',BLUE,LineWidth=2); leg = [leg,"d=2 lower bound"];
x = r.m1d3_050; res = h.rps(x,'ps',1,'T',T); plot(x.lossDBs,res,'-','color',RED,LineWidth=2); leg = [leg,"d=3 lower bound"];
x = r.m1d4_050; res = h.rps(x,'ps',1,'T',T); plot(x.lossDBs,res,'-','color',ORANGE,LineWidth=2); leg = [leg,"d=4 lower bound"];

% upper bound
x = r.m1d2_050; res = h.rps_upper(x,T); plot(x.lossDBs,res,'--','color',BLUE,LineWidth=2); leg = [leg,"d=2 upper bound"];
x = r.m1d3_050; res = h.rps_upper(x,T); plot(x.lossDBs,res,'--','color',RED,LineWidth=2); leg = [leg,"d=3 upper bound"];
x = r.m1d4_050; res = h.rps_upper(x,T); plot(x.lossDBs,res,'--','color',ORANGE,LineWidth=2); leg = [leg,"d=4 upper bound"];

legend(leg,'Position',[0.37 0.275 0 0.05])
ylabel('Secure bits per second')
xlabel('Loss (dB)')

% yticks([1 3  5]*1e4)
% yticklabels({'1','3','5'});

ax = gca;
ax.YAxis.Exponent = 2;

%
annotation('textbox', [0.11, 1, 0, 0], 'string', 'x10^4','FontSize',     9)
fig = gcf;
fig.Units               = 'centimeters';
fig.Position(3)         = 8;
fig.Position(4)         = 8;
set(fig.Children,     'FontSize',     9);
fig.PaperPositionMode   = 'auto';
box on