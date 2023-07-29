function [by,bx,bz,b] = fun_bonyance_and_gradient(dirdata,nt,LX,LY,LZ)
%%
% boyance using linear density state equation 
% and gradient 
%=========== input====model data====================================  
theta=rdmds([dirdata 'THETA'],nt);
theta(:,1:2,:)=NaN;
theta(:,end-1:end,:)=NaN;
%===================================================================
tAlpha=-2e-4;g=9.81;
b=-g*(theta*tAlpha); % boyance using linear density state equation
[by,bx,bz]=gradient(b,LY,LX,LZ);
end



