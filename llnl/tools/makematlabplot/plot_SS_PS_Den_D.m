clc
close all
clear all
    
% User setting
LoadData = 1;
saveJPG= 1;

% Diameter = 150;
% SSrange=[0 1.5 0 1500];     % [MPa]
% SPrange=[0 0.25 0 1500];     % [MPa]
% Denrange=[0 1.5 0 300.0];     % 1e13 [m^{-2}]

% Diameter = 300;
% SSrange=[0 1.5 0 1000];     % [MPa]
% SPrange=[0 0.4 0 1000];     % [MPa]
% Denrange=[0 1.5 0 100.0];     % 1e13 [m^{-2}]

% Diameter = 600;
% SSrange=[0 1.5 0 800];     % [MPa]
% SPrange=[0 0.4 0 800];     % [MPa]
% Denrange=[0 1.5 0 15.0];     % 1e13 [m^{-2}]

% Diameter = 1000;
% SSrange=[0 0.50 0 500];     % [MPa]
% SPrange=[0 0.2 0 500];     % [MPa]
% Denrange=[0 0.50 0 5.0];     % 1e13 [m^{-2}]

mu = 4.800000e+10;
Enormal = mu*pi/16; 

% Load data
if (LoadData)
Data_1 = load(sprintf('Tension_D%d_1',Diameter));
Data_2 = load(sprintf('Tension_D%d_2',Diameter));
Data_3 = load(sprintf('Tension_D%d_3',Diameter));
Data_4 = load(sprintf('Tension_D%d_4',Diameter));
Data_5 = load(sprintf('Tension_D%d_5',Diameter));
end
C = {'k','b','r','g','c',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
% Plot TT
figure('Position',[100,100,1200,800]);
L1=plot(Data_1(:,1)*100,Data_1(:,3)./1e6,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,1)*100,Data_2(:,3)./1e6,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,1)*100,Data_3(:,3)./1e6,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,1)*100,Data_4(:,3)./1e6,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,1)*100,Data_5(:,3)./1e6,'color',C{5},'LineWidth',2);hold on
% L6=plot(Data_5(:,1)*100,Data_5(:,7)*Enormal./1e6,'-.','color',C{1},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,'Location','NorthWest')   
%     ,sprintf('Elastic')...    
    xlabel('\epsilon [%]', 'Fontsize',20,'FontName','Times New Roman')
    ylabel('\sigma [MPa]', 'Fontsize',20,'FontName','Times New Roman')
    set(gca,'FontSize',20,'FontName','Times New Roman') ;
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(SSrange)

if (saveJPG)
        hgexport(gcf, sprintf('SS_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end

% Plot TP
figure('Position',[100,1200,1200,800]);
L1=plot(Data_1(:,2)*100,Data_1(:,3)./1e6,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,2)*100,Data_2(:,3)./1e6,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,2)*100,Data_3(:,3)./1e6,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,2)*100,Data_4(:,3)./1e6,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,2)*100,Data_5(:,3)./1e6,'color',C{5},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,'Location','SouthEast')   
xlabel('\epsilon_P [%]', 'Fontsize',20,'FontName','Times New Roman')
ylabel('\sigma [MPa]', 'Fontsize',20,'FontName','Times New Roman')
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(SPrange)

if (saveJPG)
        hgexport(gcf, sprintf('SP_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end

% Plot Density
figure('Position',[1500,100,1200,800]);
L1=plot(Data_1(:,1)*100,Data_1(:,4)./1e13,'color',C{1},'LineWidth',2);hold on
L2=plot(Data_2(:,1)*100,Data_2(:,4)./1e13,'color',C{2},'LineWidth',2);hold on
L3=plot(Data_3(:,1)*100,Data_3(:,4)./1e13,'color',C{3},'LineWidth',2);hold on
L4=plot(Data_4(:,1)*100,Data_4(:,4)./1e13,'color',C{4},'LineWidth',2);hold on
L5=plot(Data_5(:,1)*100,Data_5(:,4)./1e13,'color',C{5},'LineWidth',2);hold on
legend([L1 L2 L3 L4 L5]...
    ,sprintf('case 1')...
    ,sprintf('case 2')...
    ,sprintf('case 3')...
    ,sprintf('case 4')...
    ,sprintf('case 5')...
    ,'Location','SouthEast')  
    xlabel('\epsilon [%]', 'Fontsize',20,'FontName','Times New Roman')
    ylabel('\rho[10^{13}m^{-2}]', 'Fontsize',20,'FontName','Times New Roman')

% title(sprintf('FCC under torsion'),'Fontsize',20);
set(gca,'FontSize',20,'FontName','Times New Roman') 
grid on  
axis(Denrange)

if (saveJPG)
        hgexport(gcf, sprintf('Den_%d',Diameter),  ...
        hgexport('factorystyle'), 'Format', 'jpeg'); 
end
