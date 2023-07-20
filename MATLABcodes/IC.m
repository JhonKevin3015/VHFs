% Initial Conditions ========================================
%Single front in center of domain, following Fox-Kemper et al. 2008a
%% ===========
Dirdata='../DATA/'
load([Dirdata 'IC.mat'])
% theta and u from 2d to 3d
for i=1:nx % nx=300
    theta3d(i,:,:)=theta2d;
    u3d(i,:,:)=u2d;
end
%===============================Plots========
figure ()
contourf([0:deltaxy:(ny-1)*deltaxy]/1000-100,zc,theta2d'), shading flat
xlabel('Y-axis (Km)')
ylabel('Depth (m)')
title('Potential temperature (C)')
colorbar
%===============================Plots========
figure ()
contourf([0:deltaxy:(ny-1)*deltaxy]/1000-100,zc,u2d'), shading flat
xlabel('Y-axis (Km)')
ylabel('Depth (m)')
title('Zonal velocity (m/s)')
colorbar
%=================================================