clear all; 
close all;
load dsddata_20050513.mat
len=length(dsd_data(1,1,:));
dD=0.2;
Ct_all=squeeze(dsd_data(:,2,:));
D_all=squeeze(dsd_data(:,3,:));
Nd_all=squeeze(dsd_data(:,6,:));
%for n=1:len 
  %D=D_all(:,n);
  %Nd=Nd_all(:,n);
%end
D=D_all(:,1);
Nd=Nd_all(:,1);

% Surface distribution A(D)
Ad = pi.*D_all.^2 .* Nd_all;

% Mass distribution M(D)
Md = pi/6 * D_all.^3 .* Nd_all;

% Reflectivity distribution Z(D)
Zd = D_all.^6 .* Nd_all;

% create the time matrix. 481 data points between
% 7 and 12 o'clock
t1 = linspace(7,12,481);
t = zeros(41,481);
for k=1:41
  t(k,:) = t1;
end

% now plot the whole thing
pcolor([t1], [D_all], [Ct_all]);
colormap(jet);
colorbar;

