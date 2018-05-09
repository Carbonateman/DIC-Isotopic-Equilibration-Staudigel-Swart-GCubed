%%%%%%%
% D47 mixing model
% Graphs the mixing line for endmember d13C, d18O and D47 values
% Calculates mixing hyperbola for d18OC13 using C-13 concentration
% Script written by Philip Staudigel
% May 2018
%%%%%%%

clear all
%% Input end-member values here

d13CA = -15; % d13C of end-member A
d13CB = 0; % d13C of end-member B

d18OA = +15; % d18O of end-member A
d18OB = 0; % d18O of end-member B

D47A  = .7; % D47 of end-member A
D47B  = .75; % D47 of end-member B

%% Calculate Concentrations of 13C and d18OC13 values for end-members
r13A = (d13CA/1000+1)*0.0118; %Chang and Li, 1990
r13B = (d13CB/1000+1)*0.0118;

C13A = r13A/(1+r13A);
C13B = r13B/(1+r13B);

d18O13A= unclump(d13CA, d18OA, D47A)
d18O13B= unclump(d13CB, d18OB, D47B);

%% Fraction of component 1 for model
FA = [0:.1:1];

%% Simple mixing of d13C and d18O
d13Cm = d13CA*FA + d13CB*(1-FA);
d18Om = d18OA*FA + d18OB*(1-FA);

%% Calculate mixing hyperbola for d18OC13
r13m = (d13Cm./1000+1).*0.0118; %Chang and Li, 1990
C13m = r13m./(1+r13m);

d18O13m = (FA*d18O13A*C13A + (1-FA)*d18O13B*C13B)./C13m;

D47m = clump(d13Cm, d18Om, d18O13m);


%% Plot the mixed system ?47 relative to d18O
close all
plot(d18Om,D47m,'color',[0 0 0])  % plot the line
hold on
plot(d18Om,D47m,'o','markeredgecolor',[0 0 0],'markerfacecolor',[.7 .7 .7],'markersize',10) % plot the dots

% make it pretty
set(gca,'tickdir','both','linewidth',1,'xcolor','black','ycolor','black','fontsize',16);

% axis labels
xlabel('\delta^{18}O (  ^{\fontsize{8}o}/{\fontsize{8}o o} VPDB)','fontname','calibri','fontsize',16);
ylabel('\Delta_{47}(^{ \fontsize{8}o}/{\fontsize{8}o o})','fontname','calibri','fontsize',16);

%% Two scripts for clumping / unclumping

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