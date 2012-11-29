%
% [x,y,z,info]=csdp(At,b,c,K,pars)
%
% Uses CSDP to solve a problem in SeDuMi format.
%
% Input:
%        At, b, c, K      SDP problem in SeDuMi format.
%        pars             CSDP parameters (optional parameter.)
%
% Output:
%
%        x, y, z          solution.
%        info             CSDP return code.
%                         info=100 indicates a failure in the MATLAB
%                         interface, such as inability to write to 
%                         a temporary file.
%
% Note: This interface makes use of temporary files with names given by the
% tempname function.  This will fail if there is no working temporary
% directory or there isn't enough space available in this directory.  
%
% Note: This code writes its own param.csdp file in the current working
% directory.  Any param.csdp file already in the directory will be deleted.
%
% Note: It is assumed that csdp is the search path made available through
% the "system" or "dos" command.  Typically, having the csdp executable in
% current working directory will work, although some paranoid system
% administrators keep . out of the path.  In that case, you'll need to
% install csdp in one of the directories that is in the search path.
% A simple test is to run csdp from a command line prompt.  
%
function [x,y,z,info]=csdp(At,b,c,K,pars)
%
% First, put a dummy pars in place if no argument was given.  Also
% set pars.printlevel if not given.
%
if (nargin<5)
  pars.printlevel=1;
else
  if (isfield(pars,'printlevel'))
    pars.printlevel=pars.printlevel;
  else
    pars.printlevel=1;
  end
end
%
% Write out the param.csdp file.
%
fid=fopen('param.csdp','w');
if (fid==-1)
  if (pars.printlevel ~= 0) 
    fprintf('Could not open param.csdp\n');
  end
  info=100;
  return 
end
%
% Now, go through the parameters.
%

if (isfield(pars,'axtol'))
  fprintf(fid,'axtol= %e\n',pars.axtol);
else
  fprintf(fid,'axtol=%e\n',1.0e-8);
end

if (isfield(pars,'atytol'))
  fprintf(fid,'atytol= %e\n',pars.atytol);
else
  fprintf(fid,'atytol=%e\n',1.0e-8);
end

if (isfield(pars,'objtol'))
  fprintf(fid,'objtol= %e\n',pars.objtol);
else
  fprintf(fid,'objtol=%e\n',1.0e-8);
end

if (isfield(pars,'pinftol'))
  fprintf(fid,'pinftol= %e\n',pars.pinftol);
else
  fprintf(fid,'pinftol=%e\n',1.0e8);
end

if (isfield(pars,'dinftol'))
  fprintf(fid,'dinftol= %e\n',pars.dinftol);
else
  fprintf(fid,'dinftol=%e\n',1.0e8);
end

if (isfield(pars,'maxiter'))
  fprintf(fid,'maxiter= %d\n',pars.maxiter);
else
  fprintf(fid,'maxiter=%d\n',100);
end
if (isfield(pars,'minstepfrac'))
  fprintf(fid,'minstepfrac= %e\n',pars.minstepfrac);
else
  fprintf(fid,'minstepfrac=%e\n',0.90);
end

if (isfield(pars,'maxstepfrac'))
  fprintf(fid,'maxstepfrac= %e\n',pars.maxstepfrac);
else
  fprintf(fid,'maxstepfrac=%e\n',0.97);
end

if (isfield(pars,'minstepp'))
  fprintf(fid,'minstepp= %e\n',pars.minstepp);
else
  fprintf(fid,'minstepp=%e\n',1.0e-8);
end

if (isfield(pars,'minstepd'))
  fprintf(fid,'minstepd= %e\n',pars.minstepd);
else
  fprintf(fid,'minstepd=%e\n',1.0e-8);
end

if (isfield(pars,'usexzgap'))
  fprintf(fid,'usexzgap= %d\n',pars.usexzgap);
else
  fprintf(fid,'usexzgap=%d\n',1);
end

if (isfield(pars,'tweakgap'))
  fprintf(fid,'tweakgap= %d\n',pars.tweakgap);
else
  fprintf(fid,'tweakgap=%d\n',0);
end

if (isfield(pars,'affine'))
  fprintf(fid,'affine= %d\n',pars.affine);
else
  fprintf(fid,'affine=%d\n',0);
end

if (isfield(pars,'printlevel'))
  fprintf(fid,'printlevel= %d\n',pars.printlevel);
else
  fprintf(fid,'printlevel=%d\n',1);
end

if (isfield(pars,'perturbobj'))
  fprintf(fid,'printlevel= %d\n',pars.perturbobj);
else
  fprintf(fid,'printlevel=%d\n',1);
end

if (isfield(pars,'fastmode'))
  fprintf(fid,'printlevel= %d\n',pars.fastmode);
else
  fprintf(fid,'printlevel=%d\n',0);
end


%
% close the parameter file.
%
fclose(fid);
%
% Write the problem out.
%
fname=tempname;
ret=mvu_writesdpa([fname '.dat-s'],At,b,c,K,pars);
if (ret==1)
  info=100;
  x=[];
  y=[];
  z=[];
  return
end
%
% Solve the problem.
%
if (ispc==1)
  info=dos(['csdp ' fname '.dat-s' ' ' fname '.sol'],'-echo');
else
  if (pars.printlevel ~=0)
    info=system(['time csdp ' fname '.dat-s' ' ' fname '.sol']);
  else
    info=system(['csdp ' fname '.dat-s' ' ' fname '.sol']);
  end
end
%
% Read back the solution.
%
[x,y,z]=mvu_readsol([fname '.sol'],K,length(b));
%
% Delete the temporary files, including param.csdp if we wrote one!
%
delete([fname '.dat-s']);
delete([fname '.sol']);
delete('param.csdp');
