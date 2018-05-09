%%%%%%%
% Dual rate model for equilibrating DIC solutions
% Description: Reads experimental data and solves for equilibrium and initial 
% isotopic compositions for D47 and d18O values. Models equilibrating d18O
% values using a first-order exponential fit, models D47 values using a
% kinetic difference in equilibration rates for C12 and C13 (KIE). 
% Output: Figure 2 (single expeiment), initial and equilibrium
% compositions, equilibration rate, KIE. All output values are ± 95%,
% unless the parameter 'Conf' is changed.
% Script written by: Philip Staudigel
% Date: May 2018
%%%%%%%


%% Load Data
COTime = [1 20 60 120 210 300 1230 2790 4320 7200 10134 14514]
d18O = [-18.565 -18.036 -16.939 -15.356 -14.384 -13.406 -7.130 -2.138 0.013 0.549 1.103 1.080];
d13C = [-16.053 -15.974 -15.480 -15.300 -16.027 -15.941 -15.677 -15.632 -15.272 -15.715 -15.174 -14.956];

D47Time = [1 20 120 210 300 1230 2790 4320 14514];
D47  = [0.566 0.565 0.516 0.523 0.525 0.543 0.634 0.690 0.700];
D47err = [0.013 0.011 0.010 0.005 0.011 0.039 0.003 0.001 0.042];

COLOR = [1 0 0];
SYMBOL = 'o';
SYMBOLSIZE = 10

Conf = 0.05 ; % desired confidence interval for statistics "Alpha"

% Fits function d18O(time)=d18Oeq - (d18Oeq-d18Oi)*exp(-lambda* time)
d18Ofun = @(OC,time)OC(1)-(OC(1)-OC(2))*exp(-OC(3)*time);
OC0 = [0 -18 .007]; % initial values for constants, Eq d18O, initial d18O and lambda
% Values of [0 -18 and .007] works for all experimental results presented
% herein, but may need to be refined to fit other measurements. These
% values are a coarse approximation which is refined using the nlinfit
% function
[d18OC,ro,Jo,covo,mseo] = nlinfit(COTime,d18O,d18Ofun,OC0); % Calculate best fits for constants
d18Oci = nlparci(d18OC,ro,'jacobian',Jo,'alpha',Conf); % calulate 95% confidence for constants
d18OC % These are the refined parameters from OC0
d18O95 = d18Oci(:,2)-d18OC(:) % 95% (or other) CI for each parameter for d18O model

d13Cm=mean(d13C);



%% Fit model to calculate equilibrium and initia D47 and KIE for Mass-13 
% Cfun is simply a paired d18O and d18O13 model which outputs D47 relative
% to time using inputs CC [EQ INIT KIE]
Cfun = @(CC,time)clump(d13Cm,d18Ofun(d18OC,time),d18Ofun([unclump(d13Cm,d18OC(1),CC(1)) unclump(d13Cm,d18OC(2),CC(2)) d18OC(3)*(CC(3))],time));

% Initial parameters CC0 are equilibrium D47, Initial D47 and Kinetic
% isotope effect for D47. Initial values of [.78 .64 and .02] work for all
% experimental results, and are tuned using the nlinfit function.
CC0 = [.78 .64 .98]; % initial values for constants, Eq d18O, initial d18O and KIE

[CC,rc,Jc,covc,msec] = nlinfit(D47Time,D47,Cfun,CC0); % Calculate best fits for constants
Cci = nlparci(CC,rc,'jacobian',Jc,'alpha',Conf); % calulate 95% confidence for constants
CC
C95 = Cci(:,2)-CC(:) 

MT = ([min(D47Time):1:max(D47Time)]);
YO = d18Ofun(d18OC,MT);
YC = Cfun(CC,MT);
YD = Cfun([CC(1), CC(2), 1],MT);

% close all
% figure('Units','pixels','position',[00 300 500 375]);

subplot(2,1,1)
semilogx(COTime,d18O,'.')
hold on
line(MT,YO,'color',COLOR,'linewidth',2);

plot(COTime,d18O,SYMBOL,'markeredgecolor',[0 0 0],'markerfacecolor',COLOR,'markersize',SYMBOLSIZE);
set(gca,'tickdir','both','linewidth',1,'xcolor','black','ycolor','black','fontsize',12);

subplot(2,1,2)

semilogx(D47Time,D47,'.');
hold on
errorbar(D47Time,D47,D47err,'linewidth',1,'linestyle','none','color',[0 0 0])
semilogx(D47Time,D47,SYMBOL,'markeredgecolor',[0 0 0],'markerfacecolor',COLOR,'markersize',SYMBOLSIZE);
plot(MT,YC,'color',COLOR,'linewidth',2);
plot(MT,YD,'color',[COLOR 0.5],'linewidth',2,'linestyle','--')
set(gca,'tickdir','both','linewidth',1,'xcolor','black','ycolor','black','fontsize',12);


function [D47] = clump(d13C,d18O,d18O13)
% Values for PDB and lambda
R13PDB = 0.0118; %Chang and Li, 1990
R18PDB = 0.00208839; % Modified from Baertschi 1976
R17PDB = 0.0003931; % Assonov and Brenninkmeijer, 2003
l = 0.528; % Lambda: Barkan and Luz, 2005. R17 = (R18s/R18pdb)^l x R17pdb


r13 = (d13C./1000+1).*R13PDB;
r18 = (d18O./1000+1).*R18PDB;
r17 = R17PDB*(r18/R18PDB).^l;
r1813 = (d18O13./1000+1).*R18PDB;
r1713 = R17PDB.*(r1813./R18PDB).^l;

C12 = 1./(1+r13);
C13 = r13./(1+r13);

C16 = 1./(1+r17+r18);
C17 = r17./(1+r17+r18);
C18 = r18./(1+r17+r18);

C1613 = 1./(1+r1713+r1813);
C1713 = r1713./(1+r1713+r1813);
C1813 = r1813./(1+r1713+r1813);

R47stoc = (C13.*C17.*C17 + 2.*C13.*C18.*C16 + 2.*C12.*C17.*C18)./(C12.*C16.*C16);
R47real = (C13.*C1713.*C1713 + 2.*C13.*C1813.*C1613 + 2.*C12.*C17.*C18)./(C12.*C16.*C16);

D47 = (R47real./R47stoc-1).*1000;

end
function [d18O13] = unclump(d13C,d18O,D47)
 steps = 5; % number of iterations for solving d18O13
 d18O13 = d18O+D47; % an initial approximation for d18O13, to be refined iteratively
 for i=1:steps
     Diff = clump(d13C,d18O,d18O13)-D47;
     d18O13 = d18O13-Diff;
 end
 
end


