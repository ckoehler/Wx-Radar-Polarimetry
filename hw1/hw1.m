clear all; 

% toggle problem graphs
p1 = false;
p2 = true;





close all;
load dsddata_20050513.mat
len=length(dsd_data(1,1,:));
dD=0.2;
Ct_all=squeeze(dsd_data(:,2,:));
D_all=squeeze(dsd_data(:,3,:));
Nd_all=squeeze(dsd_data(:,6,:));
Vel_all=squeeze(dsd_data(:,9,:));
%for n=1:len 
  %D=D_all(:,n);
  %Nd=Nd_all(:,n);
%end
D=D_all(:,1);
Nd=Nd_all(:,1);

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
  subplot(4,1,1);
  plot(t1, N_t);
  xlabel('Time, UTC');
  ylabel('# (m^{-3})');
  title('Total Number Concentration');

  subplot(4,1,2);
  plot(t1, W);
  xlabel('Time, UTC');
  ylabel('# (g m^{-3})');
  title('Water Content');

  subplot(4,1,3);
  plot(t1, Z);
  xlabel('Time, UTC');
  ylabel('(mm^6 m^{-3})');
  title('Reflectivity factor');

  subplot(4,1,4);
  plot(t1, R);
  xlabel('Time, UTC');
  ylabel('(mm hr^{-1})');
  title('Rainfall rate');
end
