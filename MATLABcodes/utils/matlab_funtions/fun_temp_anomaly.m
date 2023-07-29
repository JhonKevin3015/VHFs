function [Tprima] = fun_temp_anomaly(dirdata,dirPC,nt,caso)
CAMPO='THETA';
dat=rdmds([dirdata CAMPO],nt);
dat(:,1:2,:)=NaN;
dat(:,end-1:end,:)=NaN;
%=======new-NOva===========================================================
[nx ny nz]=size(dat);
VEC(1:nx,1:ny,1:nz)=NaN;
datper(:,1)=nanmean(reshape(dat,[nx*ny,nz]),1);
%====================================
        for nnz=1:nz
            VEC(1:nx,1:ny,nnz)=datper(nnz);
        end
%====================================
Tprima=dat-VEC;
%=====================Save data============================================
%dirsave=[dirPC 'MIT-JOB-TESIS/caso' num2str(caso) '-data-p/'];
%fid=fopen([dirsave 'anomaly-' CAMPO '-' num2str(nt) '.bin'],'w','b');
%fwrite(fid,datPRIMA,'float32');fclose(fid);