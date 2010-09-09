clear all; 

% toggle problem graphs
p1 = false;
p2 = false;
p3 = false;
p4 = false;

% this is really problem 2
p5 = true;




close all;
load dsddata_20050513.mat
len=length(dsd_data(1,1,:));
dD=0.2;
Ct_all=squeeze(dsd_data(:,2,:));
D_all=squeeze(dsd_data(:,3,:));
Nd_all=squeeze(dsd_data(:,6,:));
Vel_all=squeeze(dsd_data(:,9,:));


%%%%%%%%%%%%%
%%% P 1.1 %%%
%%%%%%%%%%%%%

% Surface distribution A(D)
Ad_all = pi.*D_all.^2 .* Nd_all;

% Mass distribution M(D)
Md_all = pi/6 * D_all.^3 .* Nd_all;

% Reflectivity distribution Z(D)
Zd_all = D_all.^6 .* Nd_all;

% create the time matrix. 481 data points between
% 7 and 12 o'clock
t1 = linspace(7,12,481);

% only plot if enabled at the top
if (p1 == true)
  % now plot the whole thing
  subplot(4,1,1);
  pcolor([t1], [D_all], [log10(Nd_all)]);
  shading interp;
  colorbar;
  xlabel('Time, UTC');
  ylabel('D,mm');
  title('N(D), log_{10} (m^{-3} mm^{-1})');

  subplot(4,1,2);
  pcolor([t1], [D_all], [log10(Ad_all)]);
  shading interp;
  colorbar;
  xlabel('Time, UTC');
  ylabel('D,mm');
  title('A(D), log_{10} (mm^2 m^{-3} mm^{-1})');

  subplot(4,1,3);
  pcolor([t1], [D_all], [log10(Md_all)]);
  shading interp;
  colorbar;
  xlabel('Time, UTC');
  ylabel('D,mm');
  title('M(D), log_{10} (g m^{-3} mm^{-1})');

  subplot(4,1,4);
  pcolor([t1], [D_all], [log10(Zd_all)]);
  shading interp;
  colorbar;
  xlabel('Time, UTC');
  ylabel('D,mm');
  title('Z(D), log_{10} (mm^6 m^{-3} mm^{-1})');
end


%%%%%%%%%%%%%
%%% P 1.2 %%%
%%%%%%%%%%%%%

% total number concentration
% = 0th moment
M_0 = sum(Nd_all .* dD);
N_t = M_0; 


% water content
% related to 3rd moment
M_3 = sum(D_all.^3 .* Nd_all .* dD);
W = pi/6 * 10^-3 .* M_3;


% reflectivity factor
% = 6th moment
M_6 = sum(D_all.^6 .* Nd_all .* dD);
Z = M_6;

% rainfall rate
R = 6e-4 * pi * sum(D_all.^3 .* Vel_all .* Nd_all .* dD);

if(p2 == true)
  figure;
  subplot(4,1,1);
  plot(t1, N_t);
  xlabel('Time, UTC');
  ylabel('N (m^{-3})');
  title('Total Number Concentration N');

  subplot(4,1,2);
  plot(t1, W);
  xlabel('Time, UTC');
  ylabel('W (g m^{-3})');
  title('Water Content W');

  subplot(4,1,3);
  plot(t1, Z);
  xlabel('Time, UTC');
  ylabel('Z (mm^6 m^{-3})');
  title('Reflectivity factor Z');

  subplot(4,1,4);
  plot(t1, R);
  xlabel('Time, UTC');
  ylabel('R (mm hr^{-1})');
  title('Rainfall rate R');
end



%%%%%%%%%%%%%
%%% P 1.3 %%%
%%%%%%%%%%%%%

% get the rest of the moment we need
M_2 = sum(D_all .^2 .* Nd_all .* dD);
M_4 = sum(D_all .^4 .* Nd_all .* dD);


%% N(D)
% exponential model estimator
lambda = sqrt(12 * M_2 ./ M_4);
N_0 = M_2 .* lambda.^3 ./ gamma(3);

aSlice = 51;
Nd_comp = N_0(:,aSlice) .* exp(-lambda(:,aSlice) .* D_all(:,aSlice));

% gamma model estimator
eta = M_4 .^ 2 ./ M_2 ./ M_6;
mu = ((7-11.*eta)-(eta.^2 + 14 .* eta + 1).^0.5)/2.*(eta-1);
lambda2 = (M_2 ./ M_4 .* (mu + 3).*(mu+4)).^0.5;
N_02 = ( M_2 .* lambda2.^(mu+3))./gamma(mu+3);

Nd_comp_gamma = N_02(:,aSlice) .* D_all(:,aSlice).^mu(:,aSlice) .* exp(-lambda2(:,aSlice) .* D_all(:,aSlice));


%% Surface distribution A(D)
Ad_all_comp = pi.*D_all(:,aSlice).^2 .* Nd_comp;
Ad_all_comp_gamma = pi.*D_all(:,aSlice).^2 .* Nd_comp_gamma;

% Mass distribution M(D)
Md_all_comp = pi/6 * D_all(:,aSlice).^3 .* Nd_comp;
Md_all_comp_gamma = pi/6 * D_all(:,aSlice).^3 .* Nd_comp_gamma;

% Reflectivity distribution Z(D)
Zd_all_comp = D_all(:,aSlice).^6 .* Nd_comp;
Zd_all_comp_gamma = D_all(:,aSlice).^6 .* Nd_comp_gamma;

if(p3 == true)
  figure;
  subplot(2,2,1);
  plot(D_all(:,aSlice), Nd_all(:,aSlice), '*', D_all(:,aSlice), Nd_comp, D_all(:,aSlice), Nd_comp_gamma);
  axis([0 4 0 1300]);
  legend('N(D)', 'Exponential', 'Gamma');
  title('N(D)');
  xlabel('Dropsize (mm)');
  ylabel('N(D) (m^{-3} mm^{-1})');

  subplot(2,2,2);
  plot(D_all(:,aSlice), Ad_all(:,aSlice), '*', D_all(:,aSlice), Ad_all_comp, D_all(:,aSlice), Ad_all_comp_gamma);
  %axis([0 4 0 1300]);
  legend('A(D)', 'Exponential', 'Gamma');
  title('A(D)');
  xlabel('Dropsize (mm)');
  ylabel('A(D) (mm^2 m^{-3} mm^{-1})');
  
  subplot(2,2,3);
  plot(D_all(:,aSlice), Md_all(:,aSlice), '*', D_all(:,aSlice), Md_all_comp, D_all(:,aSlice), Md_all_comp_gamma);
  %axis([0 4 0 1300]);
  legend('M(D)', 'Exponential', 'Gamma');
  title('M(D)');
  xlabel('Dropsize (mm)');
  ylabel('M(D) (g m^{-3} mm^{-1})');

  subplot(2,2,4);
  plot(D_all(:,aSlice), Zd_all(:,aSlice), '*', D_all(:,aSlice), Zd_all_comp, D_all(:,aSlice), Zd_all_comp_gamma);
  %axis([0 4 0 1300]);
  legend('Z(D)', 'Exponential', 'Gamma');
  title('Z(D)');
  xlabel('Dropsize (mm)');
  ylabel('Z(D) (mm^6 m^{-3} mm^{-1})');
end


%%%%%%%%%%%%%
%%% P 1.4 %%%
%%%%%%%%%%%%%


% We'll need M_0, 3, and 6
M_0_fit = N_0 .* lambda2.^-(mu+1) .* gamma(mu + 1);
M_3_fit = N_0 .* lambda2.^-(mu+4) .* gamma(mu + 4);
M_6_fit = N_0 .* lambda2.^-(mu+7) .* gamma(mu + 7);


% total number concentration
% = 0th moment
N_t_fit = M_0_fit; 


% water content
% related to 3rd moment
W_fit = pi/6 * 10^-3 .* M_3_fit;


% reflectivity factor
% = 6th moment
Z_fit = M_6_fit;


% rainfall rate
%R_fit = 6e-4 * pi * sum(D_all.^3 .* Vel_all .* Nd_comp_gamma .* dD);

if(p4 == true)

  figure;
  subplot(3,1,1);
  plot(t1, N_t, t1, N_t_fit);
  xlabel('Time, UTC');
  ylabel('N (m^{-3})');
  title('Total Number Concentration N');
  legend('Data','Fit');

  subplot(3,1,2);
  plot(t1, W, t1, W_fit);
  xlabel('Time, UTC');
  ylabel('W (g m^{-3})');
  title('Water Content W');
  legend('Data','Fit');

  subplot(3,1,3);
  plot(t1, Z, t1, Z_fit);
  xlabel('Time, UTC');
  ylabel('Z (mm^6 m^{-3})');
  title('Reflectivity factor Z');
  legend('Data','Fit');

  %subplot(4,1,4);
  %plot(t1, R, t1, R_fit);
  %xlabel('Time, UTC');
  %ylabel('R (mm hr^{-1})');
  %title('Rainfall rate R');
  %legend('Data','Fit');

end



%%%%%%%%%%%%%%%
%%% P 2 %%%%%%%
%%%%%%%%%%%%%%%

% we know:
% Z = a * R ^ b
%log(Z) = log(a) + b .* log(R)
%log(R) = log(c) + d .* log(Z)

x = log(R);
y = log(Z);
result = polyfit(x,y,1);


if(p5 == true)
  plot(10*log10(Z), 10*log10(R), '.');
  xlabel('Reflectivity Z (dbZ)');
  ylabel('Rainfall rate R (mm hr^{-1})');
  title('Z-R relation');
  legend('Data');

end
