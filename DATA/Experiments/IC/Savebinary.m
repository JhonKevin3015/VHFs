function Savebinary(dirsave,tem3d,ETA2d,SA3d,uvel,vvel,x,y,z,delR) 
%   cuando cambia la resolucion de debe cambiar la pared
%   de la batimetria y los archivos 
%..................SAVE INFORMACION in 2D.......................................
% tem3d, SA3d, u3d ,v3d % data for save.................
for file=1:4
if file==1, dat2d=tem3d(:,:,1);  %temperatura superficial....
    fid=fopen([dirsave 'front_SST_relax.bin'],'w','b'); end
if file==2, dat2d(1:length(x),1:length(y))=-z(end); % batimetria plana.....
    dat2d(1:length(x),1:2)=0; % ponemos pared en el sur ............
    dat2d(1:length(x),end-1:end)=0; % ponemos pared en el norte ..NEW..........
    fid=fopen([dirsave 'front_bathy.bin'],'w','b'); end
if file==3, dat2d(1:length(x),1:length(y))=0; % viento , sin viento ........
    fid=fopen([dirsave 'front_zonal_wind.bin'],'w','b'); end
if file==4, dat2d(1:length(x),1:length(y))=ETA2d(:,:); % ETA inicial ........
    fid=fopen([dirsave 'front_eta_INI.bin'],'w','b'); end
fwrite(fid,dat2d,'float32');fclose(fid);
end
%..................3D........................................
for file=1:5
if file==1,dat3d=tem3d;
    fid=fopen([dirsave 'front_temperature_INI.bin'],'w','b'); end
if file==2,dat3d=SA3d;
    fid=fopen([dirsave 'front_salinity_INI.bin'],'w','b'); end
if file==3,dat3d=tem3d*0; % Restauracion de Temp 
%    dat3d(:,end-1:end,:)=1; dat3d(:,end-4:end-3,:)=0.25; % walt sponge
    fid=fopen([dirsave 'front_T_relax_mask.bin'],'w','b'); end
if file==4,dat3d=uvel; % velocidad zonal
    fid=fopen([dirsave 'front_uvel_INI.bin'],'w','b'); end
if file==5,dat3d=vvel; % vel meridional
    fid=fopen([dirsave 'front_vvel_INI.bin'],'w','b'); end
fwrite(fid,dat3d,'float32');fclose(fid);
end
%..................3D........................................
for file=1:2
if file==1,dat3d=tem3d*0; % % Restauracion de U vel
%    dat3d(:,end-1:end,:)=1; dat3d(:,end-4:end-3,:)=0.25; % walt sponge
    fid=fopen([dirsave 'front_U_relax_mask.bin'],'w','b'); end
if file==2,dat3d=tem3d*0; % % Restauracion de V vel
%    dat3d(:,end-1:end,:)=1; dat3d(:,end-4:end-3,:)=0.25; % walt sponge
    fid=fopen([dirsave 'front_V_relax_mask.bin'],'w','b'); end
fwrite(fid,dat3d,'float32');fclose(fid);
end
% save deltaZ
fid=fopen([dirsave 'INI_DelR.bin'],'w','b'); 
fwrite(fid,delR,'float32'); fclose(fid);