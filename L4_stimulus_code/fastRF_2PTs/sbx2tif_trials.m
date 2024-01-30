function sbx2tif_trials(fname,vol,varargin)
%HW modify sbx2tif to create multiple TIFs per trial (like scanimage)
% sbx2tif
% Generates tif file from sbx files
% Argument is the number of trials to convert
% If no argument is passed the whole file is converted
%discpath= 'H:\2PT_data\5364\20190131\'; %where you want tiffs saved
discpath='';
z = sbxread(fname,1,1);
global info;

if(nargin>2)
   N = min(varargin{1},numel(info.frame)/2); 
else
    N = floor(numel(info.frame)/2);
end
% vol=0; %if vol scan that you want converted all together (interleaved)
%vol=1; if vol scan that you want converted with each plane seperate.
rigid=[]; nslices=1;
if info.volscan~=0 && vol==1

    if exist( [fname '_ot_' num2str(0,'%03.f') '_rigid.sbx']) ==2
        rigid='_rigid'; %set to [] if not needed
    else
        rigid=[]; 
        if exist( [fname '_ot_' num2str(0,'%03.f') '.sbx'])==0
            sbxsplit(fname); %split if volumetric
            max_idx=info.max_idx;
        end
    end
    nslices = info.otparam(3);
end

%divides up total frames in session to which plane they match up to
frames=0:info.frame(end);
frames=frames';
frames_pp=mod(frames,nslices);

plane=nan(ceil(size(frames,1)/nslices),nslices);
for i=1:nslices
   plane(1:numel(find(frames_pp==i-1)),i)=frames(frames_pp==i-1);
end

%get on/off times, but throw error if extra pulse found
if info.frame(1)==0 %delete first frame if at zero (bug)
    info.frame(1)=[];
    info.line(1)=[];
    info.event_id(1)=[];
    
end

if rem(numel(info.frame),2)~=0
    disp('ODD NUMBER OF TRIALS FOUND...FIXING');
    [x,in]=min(diff(info.frame));
    info.frame(in+1)=[];
end
on=info.frame(1:2:end); off=info.frame(2:2:end);
disp('frame on/off/diffs'); [on off off-on]
disp(['total trials =' int2str(numel(on))]);
if min(off-on)<80
    disp('SOMETHING STILL WRONG, CHECK FRAMES');
    return;
end
%save([fname '_info'],'on','off','info'); %save info about sbx

%%
tic;
N
otplanes={};
for n=1:nslices %loop through each plane
    k = 1; %index of which trial we are on
    done = 0;
    i=1;
while(~done && k<=N)
    try
        if k==1
            otplanes{n}=[]; 
        end
            
        in=find(plane(:,n)==min(plane((plane(:,n)-on(k))>=0,n)));
        in2=find(plane(:,n)==max(plane((off(k)-plane(:,n))>=0,n)));
        otplanes{n}=vertcat(otplanes{n},[in;in2]); 
        
        if vol==1 && info.volscan~=0
            fn=[fname '_ot_' num2str(n-1,'%03.f') rigid];
        else
            fn=fname;
        end
        q = squeeze(sbxread(fn,in-1,in2-in-1)); %-1 bc index starts at 0
        q = uint16(smoothdata(q,3,'gaussian',3)); 
%         q = sbxread(fname,k,1);
%         q = squeeze(q(1,:,:));
        if(k==1)
            %create new folder if on first trial
            if exist([discpath fn '_trials'],'dir')~=7
                mkdir([discpath fn '_trials']);
            end
        end
            for i=1:size(q,3)
                if(i==1)
                    imwrite(q(:,:,i),[discpath fn '_trials/' fn '_' num2str(k,'%03d') '.tif'],'tif');
                else
                    imwrite(q(:,:,i),[discpath fn '_trials/' fn '_' num2str(k,'%03d') '.tif'],'tif','writemode','append');
                end
            end

    catch e %e is an MException struct
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        
         disp(['something happened at plane=' int2str(n) ', k trial=' int2str(k)...
             ' and i slice =' int2str(i) ]);
         done = 1;
    end
    k = k+1;
    i=1; %reset i
    if rem(k,50)==0
        disp(['working on plane ' int2str(n) ', trial ' int2str(k) '/' int2str(N)]); toc;
    end
end

end

end
