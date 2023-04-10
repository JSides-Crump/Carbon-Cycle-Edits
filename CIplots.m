
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Need to transpose initially to get Confidence intervals
d13CCI = sens.Aiso.';
CO2CI = sens.CO2atm.';
O2ACI = sens.O2_A.';
O2DPCI = sens.O2_DP.';
Dist_Preac_BurialCI = sens.Dist_Preac_Burial.';
SRP_DPCI = sens.SRP_DP.';
fanoxicdistCI = sens.fanoxicdist.';
GASTCI = sens.GAST.';
forgCI = sens.forg.';

% FoxidwCI = sens.Foxidw.';
% FcarbwCI = sens.Fcarbw.' ;
% FccdegCI = sens.Fccdeg.' ;
% FrgfCI = sens.Frgf.' ;
% FmocbCI = sens.Fmocb.' ;
% FmccbCI = sens.Fmccb.' ;
% FsfwCI = sens.Fsfw.' ;
% FlocbCI = sens.Flocb.' ;


d13C95quant = quantile(d13CCI,[0.05 0.95]);
CO295quant = quantile(CO2CI,[0.05 0.95]);
O2A95quant = quantile(O2ACI,[0.05 0.95]);
O2DP95quant = quantile(O2DPCI,[0.05 0.95]);
Dist_Preac_Burial95quant = quantile(Dist_Preac_BurialCI,[0.05 0.95]);
SRP_DP95quant = quantile(SRP_DPCI,[0.05 0.95]);
fanoxicdist95quant = quantile(fanoxicdistCI,[0.05 0.95]);
GAST95quant = quantile(GASTCI,[0.05 0.95]);
forg95quant = quantile(forgCI,[0.05 0.95]);

d13C95median = quantile(d13CCI,[0.5]);
CO295median = quantile(CO2CI,[0.5]);
O2A95median = quantile(O2ACI,[0.5]);
O2DP95median = quantile(O2DPCI,[0.5]);
Dist_Preac_Burial95median = quantile(Dist_Preac_BurialCI,[0.5]);
SRP_DP95median = quantile(SRP_DPCI,[0.5]);
fanoxicdist95median = quantile(fanoxicdistCI,[0.5]);
GAST95median = quantile(GASTCI,[0.5]);
forg95median = quantile(forgCI,[0.5]);

% Foxidwmedian = quantile(FoxidwCI,[0.5]) ;
% Fcarbwmedian = quantile(FcarbwCI,[0.5]) ;
% Fccdegmedian = quantile(FccdegCI,[0.5]) ;
% Frgfmedian = quantile(FrgfCI,[0.5]) ;
% Fmocbmedian = quantile(FmocbCI,[0.5]) ;
% Fmccbmedian = quantile(FmccbCI,[0.5]) ;
% Fsfwmedian = quantile(FsfwCI,[0.5]) ;
% Flocbmedian = quantile(FlocbCI,[0.5]) ;


% % semilogy((sens.time_myr),(atmoso295quant/3.7e19),'linewidth',1,'color',c_mean)


global starting
%%%%%% define colours
c_mean = [0.2 0.6 0.6] ;
c_std = [0.3 0.7 0.7] ;
c_range = [ 0.4 0.8 0.8] ;

%%%% output to screen
fprintf('running sens plotting script... \t')
tic

%%%% make column vector
sens.time_myr = sens.time(:,1) /1e6 ;
%% Call in data
figure('Color',[0.80 0.80 0.70])

%d13C
load('Havigd13C.mat')
d13Ctime = -1 * Havigd13C(:,1);
d13Ccalcite = Havigd13C(:,2);
d13Cdolomite = Havigd13C(:,3);
d13Cother = Havigd13C(:,4);

subplot(4,2,1)
scatter(d13Ctime, d13Ccalcite,5,'o')
hold on
scatter(d13Ctime, d13Cdolomite,5,'+')
hold on
scatter(d13Ctime, d13Cother,5,'*')
hold on
% plot((sens.time_myr),mean((sens.Aiso),2),'linewidth',3,'color',c_mean)
plot((sens.time_myr),(d13C95median),'linewidth',3,'color',c_mean)
hold on
plot((sens.time_myr),(d13C95quant),'linewidth',1,'color',c_range)

% plot((sens.time_myr),max((sens.Aiso),[],2),'linewidth',2,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.Aiso),[],2),'linewidth',2,'color',c_range)
xlim([-4e3 0]);
ylim([-25 20])
legend('calcite','dolomite','other')
%CO2
load('CO2proxy')


subplot(4,2,2)
% hold on
box on
xlabel('Time (Ma)')
ylabel('CO2')
%%%% plot data comparison
% plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
%%%% plot this model
% semilogy((sens.time_myr),mean((sens.CO2atm*1e6),2),'linewidth',3,'color',c_mean)
semilogy((sens.time_myr),(CO295median*1e6),'linewidth',3,'color',c_mean)

hold on
semilogy((sens.time_myr),(CO295quant*1e6),'linewidth',1,'color',c_range)

% semilogy((sens.time_myr),max((sens.CO2atm*1e6),[],2),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),min((sens.CO2atm*1e6),[],2),'linewidth',0.5,'color',c_range)
hold on
semilogy(1e3*HighTime, HighCO2*1e6,'linewidth',2,'color','k')
hold on
semilogy(1e3*LowTime, LowCO2*1e6,'linewidth',2,'color','k')
xlim([-4e3 0]);
ylim([10 1e7]);

title('CO2')

%Atmos O2
load('AtmosO2proxy.mat')
O2_A_min = min((sens.O2_A/3.7e19),[],2) ;
O2_A_max = max((sens.O2_A/3.7e19),[],2) ;

subplot(4,2,3)
% hold on
box on
xlabel('Time (Ma)')
ylabel('PAL O2')
%%%% plot data comparison
% plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
%%%% plot this model
semilogy((sens.time_myr),(O2A95quant/3.7e19),'linewidth',1,'color',c_range)

% semilogy((sens.time_myr),(O2_A_max),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),(O2_A_min),'linewidth',0.5,'color',c_range)
% hold on
% patch([sens.time_myr fliplr(sens.time_myr)], [O2_A_min fliplr(O2_A_max)],
% 'g') %% NEED TO REMOVE SPIN UP FROM DATASET IN ORDER TO USE
hold on
% semilogy((sens.time_myr),mean((sens.O2_A/3.7e19),2),'linewidth',3,'color',c_mean)
semilogy((sens.time_myr),(O2A95median/3.7e19),'linewidth',3,'color',c_mean)

hold on
semilogy(AtmosO2proxtime, AtmosO2proxlow,'linewidth',2,'color','k')
hold on
semilogy(AtmosO2proxtime, AtmosO2proxhigh,'linewidth',2,'color','k')
xlim([-4e3 0]);
ylim([1e-7 10])
title('O2 atmosphere PAL')

subplot(4,2,4)
box on
xlabel('Time (Ma)')
ylabel('deep O2 relative')
% semilogy((sens.time_myr),mean((sens.O2_DP/2.21e17),2),'linewidth',3,'color','b')
semilogy((sens.time_myr),(O2DP95median/2.21e17),'linewidth',3,'color',c_mean)

hold on
semilogy((sens.time_myr),(O2DP95quant/2.21e17),'linewidth',1,'color',c_range)

% semilogy((sens.time_myr),max((sens.O2_DP/2.21e17),[],2),'linewidth',0.5,'color','b')
% hold on
% semilogy((sens.time_myr),min((sens.O2_DP/2.21e17),[],2),'linewidth',0.5,'color','b')
xlim([-4e3 0]);
ylim([1e-9 10])
title('Relative O2 Deep')

%P burial
load('ReinhardPdata.mat')

XXX = [Time Reinhard_P] ;
ZZ = XXX(~any(isnan(XXX)| isinf( XXX ), 2 ),: ) ;
x = ZZ(:,1) ;
y = ZZ(:,2) ;
edges = (0:50:4e3);
[~,~,loc]=histcounts(x,edges);
meany = accumarray(loc(:),y(:))./accumarray(loc(:),1);
xmid = 0.5*(edges(1:end-1)+edges(2:end));
xmid = xmid.' ;
xmid = xmid(1:end-10,:);

%forg
% subplot(3,2,4)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('PAL O2')
% 
% % patch([sens.time_myr fliplr(sens.time_myr)], [O2_A_min fliplr(O2_A_max)],
% % 'g') %% NEED TO REMOVE SPIN UP FROM DATASET IN ORDER TO USE
% plot((sens.time_myr),mean(sens.forg,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.forg,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.forg,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% ylim([0 1])
% title('forg')
% figure
subplot(4,2,5)
box on
xlabel('Time (Ma)')
yyaxis right
scatter((-1*Time), logReinhard_P,5)
hold on
plot((-1*xmid),log10(meany),'linewidth',2,'color','k')
ylim([-3 1])

yyaxis left
set(gca, 'YScale', 'log')
ylabel('P burial')
% plot((sens.time_myr),mean(sens.Dist_Preac_Burial,2),'linewidth',3,'color',c_mean)
plot((sens.time_myr),(Dist_Preac_Burial95median),'linewidth',3,'color',c_mean)

hold on
plot((sens.time_myr),(Dist_Preac_Burial95quant),'linewidth',1,'color',c_range)

% plot((sens.time_myr),max((sens.Dist_Preac_Burial),[],2),'linewidth',2,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.Dist_Preac_Burial),[],2),'linewidth',2,'color',c_range)
ylim([1e7 1e13])

xlim([-4e3 0]);
title('Total P burial')

%P deep
subplot(4,2,6)
% hold on
box on

% patch([sens.time_myr fliplr(sens.time_myr)], [O2_A_min fliplr(O2_A_max)],
% 'g') %% NEED TO REMOVE SPIN UP FROM DATASET IN ORDER TO USE
% semilogy((sens.time_myr),mean((sens.SRP_DP/2790e12),2),'linewidth',3,'color',c_mean)
semilogy((sens.time_myr),(SRP_DP95median/2790e12),'linewidth',3,'color',c_mean)

hold on
semilogy((sens.time_myr),(SRP_DP95quant/2790e12),'linewidth',1,'color',c_range)

% semilogy((sens.time_myr),max((sens.SRP_DP/2790e12),[],2),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),min((sens.SRP_DP/2790e12),[],2),'linewidth',0.5,'color',c_range)
xlim([-4e3 0]);
ylim([1e-3 5])
title('P deep relative')



% 
% subplot(4,2,7)
% plot((sens.time_myr),mean((sens.fanoxicprox),2),'linewidth',3,'color',c_mean)
% hold on
% plot((sens.time_myr),max((sens.fanoxicprox),[],2),'linewidth',2,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.fanoxicprox),[],2),'linewidth',2,'color',c_range)
% xlim([-4e3 0]);
% ylim([0 1]);
% 
subplot(4,2,7)
% plot((sens.time_myr),mean((sens.fanoxicdist),2),'linewidth',3,'color',c_mean)
plot((sens.time_myr),(fanoxicdist95median),'linewidth',3,'color',c_mean)

hold on
% plot((sens.time_myr),max((sens.fanoxicdist),[],2),'linewidth',2,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.fanoxicdist),[],2),'linewidth',2,'color',c_range)
plot((sens.time_myr),(fanoxicdist95quant),'linewidth',1,'color',c_range)
title('fanoxicdist')
xlim([-4e3 0]);
ylim([0 1]);

subplot(4,2,8)
% plot((sens.time_myr),mean((sens.GAST-273),2),'linewidth',3,'color',c_mean)
plot((sens.time_myr),(GAST95median-273),'linewidth',3,'color',c_mean)

hold on
% plot((sens.time_myr),max((sens.fanoxicdist),[],2),'linewidth',2,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.fanoxicdist),[],2),'linewidth',2,'color',c_range)
plot((sens.time_myr),(GAST95quant-273),'linewidth',1,'color',c_range)
title('GAST')
xlim([-4e3 0]);
% ylim([-20 40]);
% %% Figure
% %%%%%%% make figure
% figure('Color',[0.80 0.80 0.70])
% 
% %%%% load geochem data
% % load('data/data.mat')
% 
% 
% %%%% d13C record
% subplot(3,4,1)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('PAL O2')
% 
% % patch([sens.time_myr fliplr(sens.time_myr)], [O2_A_min fliplr(O2_A_max)],
% % 'g') %% NEED TO REMOVE SPIN UP FROM DATASET IN ORDER TO USE
% semilogy((sens.time_myr),mean((sens.O2_A/3.7e19),2),'linewidth',1,'color',c_mean)
% hold on
% semilogy((sens.time_myr),(O2_A_max),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),(O2_A_min),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% ylim([1e-7 10])
% title('Atmos O2')
% 
% subplot(3,4,2)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('Degassing')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.D,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.D,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.D,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('Degassing')
% subplot(3,4,3)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('carbon build up')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.C,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.C,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.C,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('Carbon build up')
% subplot(3,4,4)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('anoxic C:P')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.CPanoxic,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.CPanoxic,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.CPanoxic,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('anoxic C:P')
% subplot(3,4,5)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('fbiota')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.fbiota,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.fbiota,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.fbiota,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('fbiota')
% subplot(3,4,6)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('LAND MASS')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.EXPOSED,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.EXPOSED,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.EXPOSED,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('land mass')
% 
% subplot(3,4,7)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('d13C')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.Aiso,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.Aiso,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.Aiso,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('carbon iso')
% 
% 
% 
% figure('Color',[0.80 0.80 0.70])
% 
% 
% subplot(3,3,1)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('Temp')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.GAST - 273,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.GAST -273,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.GAST -273,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('Temp C')
% 
% subplot(3,3,2)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('mocb molyr')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.Fmocb,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.Fmocb,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.Fmocb,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('mocb')
% 
% subplot(3,3,3)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('PP')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean((sens.PP_P + sens.PP_D + sens.PP_S),2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max((sens.PP_P + sens.PP_D + sens.PP_S),[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min((sens.PP_P + sens.PP_D + sens.PP_S),[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('PP')
% 
% subplot(3,3,4)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('CO2')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% semilogy((sens.time_myr),mean((sens.CO2atm*1e6),2),'linewidth',1,'color',c_mean)
% hold on
% semilogy((sens.time_myr),max((sens.CO2atm*1e6),[],2),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),min((sens.CO2atm*1e6),[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('CO2')
% 
% subplot(3,3,5)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('O2')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% semilogy((sens.time_myr),mean(sens.O2_DP,2),'linewidth',1,'color',c_mean)
% hold on
% semilogy((sens.time_myr),max(sens.O2_DP,[],2),'linewidth',0.5,'color',c_range)
% hold on
% semilogy((sens.time_myr),min(sens.O2_DP,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('Deep O2')
% 
% 
% subplot(3,3,6)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('foxic proximal')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.fanoxicprox,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.fanoxicprox,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.fanoxicprox,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('foxic proximal')
% 
% subplot(3,3,7)
% % hold on
% box on
% xlabel('Time (Ma)')
% ylabel('foxic distal')
% %%%% plot data comparison
% % plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
% %%%% plot this model
% plot((sens.time_myr),mean(sens.fanoxicdist,2),'linewidth',1,'color',c_mean)
% hold on
% plot((sens.time_myr),max(sens.fanoxicdist,[],2),'linewidth',0.5,'color',c_range)
% hold on
% plot((sens.time_myr),min(sens.fanoxicdist,[],2),'linewidth',0.5,'color',c_range)
% xlim([-4e3 0]);
% title('foxic distal')

figure
subplot(3,1,1)
plot(sens.time_myr, forg95median)
hold on
plot(sens.time_myr, forg95quant)

xlim([-4e3 0])
ylim([0 1])
title('forg')

subplot(3,1,2)
plot(sens.time_myr,Foxidwmedian)
hold on
plot(sens.time_myr,Fcarbwmedian)
plot(sens.time_myr,Fccdegmedian)
plot(sens.time_myr,Frgfmedian)
plot(sens.time_myr,Fmocbmedian)
plot(sens.time_myr,Fmccbmedian)
plot(sens.time_myr,Fsfwmedian)
plot(sens.time_myr,Flocbmedian)
title('Carbon fluxes')
legend('oxidw','carbw','ccdeg','rgf','mocb','mccb','sfw','locb')
xlim([-4e3 0])
% 
% subplot(3,1,3)
% plot(sens.time_myr,mean(sens.doxidw.'))
% hold on
% plot(sens.time_myr,mean(sens.dcarbw.'))
% plot(sens.time_myr,mean(sens.dccdeg.'))
% plot(sens.time_myr,mean(sens.drgf.'))
% plot(sens.time_myr,mean(sens.dmocb.'))
% plot(sens.time_myr,mean(sens.dmccb.'))
% plot(sens.time_myr,mean(sens.dsfw.'))
% plot(sens.time_myr,mean(sens.dlocb.'))
% legend('oxidw','carbw','ccdeg','rgf','mocb','mccb','sfw','locb')
% xlim([-4e3 0])