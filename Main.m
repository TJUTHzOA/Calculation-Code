% Title: Calculation of time-domain THz-OA signal intensities of NaCl solutions with different concentration
% Institute: Center for Terahertz Waves of Tianjin University
% Author: Ke Zhang 
% Date: 2022/4/1

clear all;
close all;
clc;

%%  Specific heat capacity, Sound velocity, Thermal coefficient of volume expansion 
c_molL = [0,0.5,0.75,1,1.5,2,3]; % Volume molar concentration [mol/L]
rho = 0.03726*c_molL + 1.001; % Density [g/cm^3]
c_per = c_molL*58.44./rho/1e3*100; % Mass percentage [%]
vs = 12.18*c_per+1492.7; % Sound velocity [m/s]
Cp = (4.157-0.036*c_per)*1e3; % Specific heat capacity [J/kg/K]
c_molKg = c_molL./rho; % Mass molar concentration [mol/kg]
beta = (257.27+80.35*c_molKg-13.991*c_molKg.^2+1.0033*c_molKg.^3)*1e-6; % Thermal coefficient of volume expansion [/K]

%% ua(v)
data = xlsread('NaClAbsorption.xlsx','Sheet1'); %Absorption coefficient of NaCl solution
fre = data(:,1)/1e3; % [THz]
for n = 1:1:7
    eval(["Abs"+num2str(n)+"=data(:,"+num2str(n+1)+");"]); % [/cm]
end
% 9th order polynomial fitting
for n = 1:1:7
    eval(["Coef"+num2str(n)+"=polyfit(fre,Abs"+num2str(n)+",9);"]);
end
% Polynomial expression
syms v;
for n = 1:1:7
    eval(["a=Coef"+num2str(n)+"(1);"]);
    eval(["b=Coef"+num2str(n)+"(2);"]);
    eval(["c=Coef"+num2str(n)+"(3);"]);
    eval(["d=Coef"+num2str(n)+"(4);"]);
    eval(["e=Coef"+num2str(n)+"(5);"]);
    eval(["f=Coef"+num2str(n)+"(6);"]);
    eval(["g=Coef"+num2str(n)+"(7);"]);
    eval(["h=Coef"+num2str(n)+"(8);"]);
    eval(["i=Coef"+num2str(n)+"(9);"]);
    eval(["j=Coef"+num2str(n)+"(10);"]);
    eval(["Ua"+num2str(n)+"=a*v^9+b*v^8+c*v^7+d*v^6+e*v^5+f*v^4+g*v^3+h*v^2+i*v^1+j*v^0;"]);
end
%% The spectrum of wide-spectrum terahertz pulse
data = xlsread('Spectrum.xlsx','Sheet1');
x = data(:,1);
y = data(:,2);
% Double Gaussian function fitting
a1 = 0.5874;b1 = 0.4871;c1 = 0.2868;
a2 = 0.6371;b2 = 0.2727;c2 = 0.1826;
Jv_temp =  a1*exp(-((v-b1)/c1).^2) + a2*exp(-((v-b2)/c2).^2); % Fitting expression
Int = int(Jv_temp,v,fre(1),fre(end)); 
Jv = 4e-6/Int*Jv_temp; % Jv(v) [J*s]
Te = Jv/1.5e-3/1.5e-3; % Fv(v) [J/m/m*s]

%% Process
Tau = zeros(1,7); % Grüneisen parameter
uF = zeros(1,7); % Specific absorption
p0 = zeros(1,7);  % Local pressure
for n = 1:1:7
    eval(["Flag = Ua"+num2str(n)+"*100*Te;"]);
    uF(n) = int(Flag,v,fre(1),fre(end));
    Tau(n) = beta(n)*vs(n)^2/Cp(n);
    p0(n) = Tau(n)*uF(n);
end

%% Plot
figure(1);
subplot(311);plot(c_molL,Cp,'k-','Linewidth',2.5);title('Specific heat capacity (J/kg/K)');xlabel('NaCl concentration (mol/L)');
subplot(312);plot(c_molL,vs,'k-','Linewidth',2.5);title('Sound velocity (m/s)');xlabel('NaCl concentration (mol/L)');
subplot(313);plot(c_molL,beta,'k-','Linewidth',2.5);title('Thermal coefficient of volume expansion (/K)');xlabel('NaCl concentration (mol/L)');

figure(2);
subplot(211);
hold on;
color = ["'r-'","'g-'","'b-'","'k-'","'y-'","'m-'","'c-'"];
for i = 1:1:7
    eval(["plot(fre,Abs"+num2str(i)+","+color(i)+",'Linewidth',2);"]);
end
xlabel('NaCl concentration (mol/L)');box on;
title('Absorption coefficient of NaCl solution (Experiment) (/cm)');
subplot(212);
hold on;
color = ["'r-'","'g-'","'b-'","'k-'","'y-'","'m-'","'c-'"];
for i = 1:1:7
    eval(["fplot(v,Ua"+num2str(i)+","+color(i)+",'Linewidth',2);"]);
end
xlabel('NaCl concentration (mol/L)');box on;
axis([fre(1),fre(end),0,300]);
title('Absorption coefficient of NaCl solution (Fitting) (/cm)');

figure(3);
hold on;
plot(x,y,'ro','Linewidth',2.5);
fplot(v,Jv_temp,'k-','Linewidth',2.5);
box on;
legend('Experiment','Fitting');
xlabel('Frequency (THz)');
title('Spectrum of terahertz pulse');
axis([0,1.8,0,1]);

figure(4);
box on;
subplot(311);plot(c_molL,Tau,'k-','Linewidth',2.5);title('Grüneisen parameter');xlabel('NaCl concentration (mol/L)');
subplot(312);plot(c_molL,uF,'k-','Linewidth',2.5);title('Specific absorption (J/m^3)');xlabel('NaCl concentration (mol/L)');
subplot(313);plot(c_molL,p0,'k-','Linewidth',2.5);title('Local pressure (Pa)');xlabel('NaCl concentration (mol/L)');
