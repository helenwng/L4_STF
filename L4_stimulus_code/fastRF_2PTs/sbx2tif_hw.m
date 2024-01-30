function sbx2tif_hw(fname,varargin)
%HW: modified so that makes a new tiff once it reaches max 1000 frames per
%tif to avoid max file sizse
% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is converted


z = sbxread(fname,1,1);
global info;

if(nargin>1)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end

k = 0;
done = 0;
j= -1;
maxsize=1000;
%create new folder if on first trial
if exist([fname '_trials'],'dir')~=7
    mkdir([fname '_trials']);
end
while(~done && k<=N)
    try
        q = sbxread(fname,k,1);
        q = squeeze(q(1,:,:));
        if(k==1) || (mod(k,maxsize)==0 && k>0) %k==1 in orig sbx2tif for some reason?
            j=j+1;
            imwrite(q,[fname '_trials/' fname '_' num2str(j,'%02d') '.tif'],'tif');
        else
            imwrite(q,[fname '_trials/' fname '_' num2str(j,'%02d') '.tif'],'tif','writemode','append');
        end
    catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        
        disp(['check fifnamele size! end of tif conversion at k= ' int2str(k) '/ N= ' int2str(N)]);
        j
        done = 1;
    end
    k = k+1;
end