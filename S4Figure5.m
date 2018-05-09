%%%%%%%%%%%%%%
% Quiver Plot for equilibrating DIC D47 and d18O
% Script written by Philip Staudigel
% May 2018
%%%%%%%%%%%%%%
clear
close all

KIE = 0.98; % Value for kinetic Isotope Effect
EQd18O = -1.38;
EQD47 = 0.689;
d13C = -16;

EQd18O13 = unclump(d13C,EQd18O,EQD47);

%% Plot Values
d18OMAX = 15; % Maximum plotted d18O value
d18OMIN = -25; % Minimum plotted d18O value
D47MAX  = .85; % Maximum plotted D47 value
D47MIN  = .45; % Minimum plotted D47 value


%% Number of vectors, creates a X-by-Y grid
d18Osteps = 11; % Number of d18O values for vectors (X)
D47steps  = 11; % Number of D47 values for vectors (Y)

dt = 0.001; % Timestep for approximating derivitive. Just needs to be small (<0.1)

%% Generate data for quiver plot
for i = 1:d18Osteps
    for j= 1:D47steps
        d18O(i,j) = d18OMIN+(i-1)*(d18OMAX-d18OMIN)/(d18Osteps-1);
        D47(i,j) = D47MIN +(j-1)*(D47MAX-D47MIN)/(D47steps-1);
        d18O13temp1 = unclump(d13C,d18O(i,j),D47(i,j));
        d18Otemp = EQd18O - (EQd18O-d18O(i,j))*exp(-dt);
        d18O13temp2 = EQd18O13 - (EQd18O13-d18O13temp1)*exp(-dt*(KIE));
        D47temp = clump(d13C,d18Otemp,d18O13temp2);
        dd18O(i,j) = d18Otemp-d18O(i,j);
        dD47(i,j) = D47temp-D47(i,j);     
    end
end

%% Real data to be displayed (In this case, experimental results for 15°C)
rd18O=[-18.98 -17.46 -15.60 -11.84 -9.25 -2.50 -1.77 -0.99]; % d18O values
rd18Oerr=[0.11 0.09 0.19 0.19 0.02 0.32 0.08 0.02]; % d18O error bars
rD47=[0.548 0.535 0.511 0.493 0.501 0.657 0.672 0.676]; % D47 values
rD47err = [0.026 0.006 0.023 0.028 0.000 0.012 0.017 0.002]; % D47 error bars

%% Model showing best fit for real data
%Best fit line 
bfd18Oi = -18.55;
bfD47i  = .55;
bfd18O13i = unclump(d13C,bfd18Oi,bfD47i)

%% Additional line (if desired)
secondline = 1; % 1 if yes,  0 if no second line is desired
ffd18Oi = +12
ffD47i = 0.7
ffd18O13i = unclump(d13C,ffd18Oi,ffD47i)

k12 = .00259;
timestep = 10;


for i=1:201
    bfd18O(i) = EQd18O - (EQd18O-bfd18Oi)*exp(-k12*timestep*(i-1));
    bfD47(i) = clump(d13C,bfd18O(i),EQd18O13-(EQd18O13-bfd18O13i)*exp(-k12*(KIE)*timestep*(i-1)));
    ffd18O(i) = EQd18O - (EQd18O-ffd18Oi)*exp(-k12*timestep*(i-1));
    ffD47(i) = clump(d13C,ffd18O(i),EQd18O13-(EQd18O13-ffd18O13i)*exp(-k12*(KIE)*timestep*(i-1)));
end


%% Plot data
figure('Units','pixels','position',[00 300 500 375]);
hold on
 ylim([D47MIN D47MAX])
 xlim([d18OMIN d18OMAX])
Quiv = quiver(d18O,D47,dd18O,dD47,.3)
 
set(Quiv,'color','blue','MaxHeadSize',0.00,'linewidth',2);
set(gca,'tickdir','both','linewidth',1,'xcolor','black','ycolor','black','fontsize',12);


dots=plot(d18O,D47,'.','color','blue','markersize',1);
bestfit = plot(bfd18O,bfD47,'LineWidth',2,'color',[0 0.7 0]);
if secondline == 1
fakefit = plot(ffd18O,ffD47,'LineWidth',2,'color',[0.7 .0 0.0]);
end
EQ = plot(EQd18O,EQD47,'o', 'MarkerEdgeColor',[0 0 0], 'Markerfacecolor', [1 0 0], 'markersize', 16);

rerr = errorbar(rd18O,rD47,rD47err,rD47err,rd18Oerr,rd18Oerr,'color',[0 0 0],'LineStyle','none');

real = plot(rd18O,rD47,'o','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'markerfacecolor',[.2 .7 .2]);
% pTitle = title('Vector map for equilibrating DIC \delta^{18}O and \Delta_{47} values')
xlabel('\delta^{18}O (  ^{\fontsize{8}o}/{\fontsize{8}o o} VPDB)','fontname','calibri','fontsize',12);
ylabel('\Delta_{47}(^{ \fontsize{8}o}/{\fontsize{8}o o})','fontname','calibri','fontsize',12);
box on


%% Two functions for clumping and unclumping
function [d18O13] = unclump(d13C,d18O,D47)
 steps = 5; % number of iterations for solving d18O13
 d18O13 = d18O+D47; % an initial approximation for d18O13, to be refined iteratively
 for i=1:steps
     Diff = clump(d13C,d18O,d18O13)-D47;
     d18O13 = d18O13-Diff;
 end
 
end

function [D47] = clump(d13C,d18O,d18O13)
% Values for PDB and lambda
R13PDB = 0.0118; %Chang and Li, 1990
R18PDB = 0.00208839; % Modified from Baertschi 1976
R17PDB = 0.0003931; % Assonov and Brenninkmeijer, 2003
l = 0.528; % Lambda: Barkan and Luz, 2005. R17 = (R18s/R18pdb)^l x R17pdb


r13 = (d13C/1000+1)*R13PDB;
r18 = (d18O/1000+1)*R18PDB;
r17 = R17PDB*(r18/R18PDB)^l;
r1813 = (d18O13/1000+1)*R18PDB;
r1713 = R17PDB*(r1813/R18PDB)^l;

C12 = 1/(1+r13);
C13 = r13/(1+r13);

C16 = 1/(1+r17+r18);
C17 = r17/(1+r17+r18);
C18 = r18/(1+r17+r18);

C1613 = 1/(1+r1713+r1813);
C1713 = r1713/(1+r1713+r1813);
C1813 = r1813/(1+r1713+r1813);

R47stoc = (C13*C17*C17 + 2*C13*C18*C16 + 2*C12*C17*C18)/(C12*C16*C16);
R47real = (C13*C1713*C1713 + 2*C13*C1813*C1613 + 2*C12*C17*C18)/(C12*C16*C16);

D47 = (R47real/R47stoc-1)*1000;
end

