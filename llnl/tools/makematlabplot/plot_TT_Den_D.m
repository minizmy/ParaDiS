clc
close all
clear all
    
% User setting
LoadData = 1;
saveJPG= 1;

% Diameter = 150;
% TTrange=[0 3.0 0 300];     % [MPa]
% TPrange=[0 0.5 0 300];     % [MPa]
% Denrange=[0 0.5 0 60];     % 1e13 [m^{-2}]
% Diameter = 300;
% TTrange=[0 3.0 0 300];     % [MPa]
% TPrange=[0 0.4 0 250];     % [MPa]
% Denrange=[0 0.4 0 20];     % 1e13 [m^{-2}]
% Diameter = 600;
% TTrange=[0 3.0 0 200];     % [MPa]
% TPrange=[0 0.4 0 200];     % [MPa]
% Denrange=[0 0.4 0 15];     % 1e13 [m^{-2}]
Diameter = 1000;
TTrange=[0 2.0 0 150];     % [MPa]
TPrange=[0 0.4 0 150];     % [MPa]
Denrange=[0 0.4 0 6];     % 1e13 [m^{-2}]

mu = 4.800000e+10;
Enormal = mu*pi/16; 

% Load data
if (LoadData)
Data_1 = load(sprintf('Torsion_D%d_1',Diameter));
Data_2 = load(sprintf('Torsion_D%d_2',Diameter));
Data_3 = load(sprintf('Torsion_D%d_3',Diameter));
Data_4 = load(sprintf('Torsion_D%d_4',Diameter));
Data_5 = load(sprintf('Torsion_D%d_5',Diameter));
end
C = {'k','b','r','g','c',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
% Plot TT
figure('Position',[100,100,1200,800]);
L1=plot(Data_1(:,7)*100,Data_1(:,8)./1e6,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,7)*100,Data_2(:,8)./1e6,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,7)*100,Data_3(:,8)./1e6,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,7)*100,Data_4(:,8)./1e6,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,7)*100,Data_5(:,8)./1e6,'color',C{5},'LineWidth',2);hold on
L6=plot(Data_5(:,7)*100,Data_5(:,7)*Enormal./1e6,'-.','color',C{1},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5 L6]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,sprintf('Elastic')...    
    ,'Location','NorthWest')   
    xlabel('\gamma_{surf} [%]', 'Fontsize',20,'FontName','Times New Roman')
    ylabel('\sigma (=T/D^3) [MPa]', 'Fontsize',20,'FontName','Times New Roman')
% title(sprintf('FCC under torsion'),'Fontsize',20);
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(TTrange)

if (saveJPG)
        hgexport(gcf, sprintf('TT_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end

% Plot TP
figure('Position',[100,1200,1200,800]);
L1=plot(Data_1(:,6)*100,Data_1(:,8)./1e6,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,6)*100,Data_2(:,8)./1e6,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,6)*100,Data_3(:,8)./1e6,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,6)*100,Data_4(:,8)./1e6,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,6)*100,Data_5(:,8)./1e6,'color',C{5},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,'Location','SouthEast')   
    xlabel('\gamma_{surf}^{P} [%]', 'Fontsize',20,'FontName','Times New Roman')
    ylabel('\sigma (=T/D^3) [MPa]', 'Fontsize',20,'FontName','Times New Roman')
% title(sprintf('FCC under torsion'),'Fontsize',20);
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(TPrange)

if (saveJPG)
        hgexport(gcf, sprintf('TP_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end

% Plot Density
figure('Position',[1500,100,1200,800]);
L1=plot(Data_1(:,6)*100,Data_1(:,9)./1e13,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,6)*100,Data_2(:,9)./1e13,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,6)*100,Data_3(:,9)./1e13,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,6)*100,Data_4(:,9)./1e13,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,6)*100,Data_5(:,9)./1e13,'color',C{5},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,'Location','SouthEast')  
    xlabel('\gamma_{surf} [%]', 'Fontsize',20,'FontName','Times New Roman')
    ylabel('\rho[10^{13}m^{-2}]', 'Fontsize',20,'FontName','Times New Roman')

% title(sprintf('FCC under torsion'),'Fontsize',20);
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(Denrange)

if (saveJPG)
        hgexport(gcf, sprintf('Den_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end
