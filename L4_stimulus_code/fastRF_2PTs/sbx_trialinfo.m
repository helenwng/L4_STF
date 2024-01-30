function [on, off,otplanes,filel]= sbx_trialinfo(info)
%HW 8/13/18 Gets trial start/stop information from scanbox frames
% TO DO 8/20/20- include line #?
%check if volume scan
if info.volscan~=0
    nslices = info.otparam(3); 
else
    nslices=1
end

%divides up total frames in session to which plane they match up to

frames=0:info.max_idx;
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

%check if info.frame number reset to zero
%if reset, add 2^16 
if sum(diff(info.frame)<0)>0 %sudden negative change in index means reset occurred
    info.frame(find(diff(info.frame)<0)+1:end)=65536+info.frame(find(diff(info.frame)<0)+1:end);
end
if rem(numel(info.frame),2)~=0
    disp('ODD NUMBER OF TRIALS FOUND...FIXING');
    [x,in]=min(diff(info.frame));
    info.frame(in+1)=[];
end
% on=info.frame(1:2:end)+info.line(1:2:end)./info.config.lines; 
% off=info.frame(2:2:end)+info.line(2:2:end)./info.config.lines;
on=info.frame(1:2:end); 
off=info.frame(2:2:end);

disp('frame on/off/diffs'); [on off off-on]
disp(['total trials =' int2str(numel(on))]);
if min(off-on)<80
    disp('SOMETHING STILL WRONG, CHECK FRAMES');
    return;
end

%%
N = floor(numel(info.frame)/2);
otplanes={};
for n=1:nslices %loop through each plane
    k = 1; %index of which trial we are on
    done = 0;
    tic;
    while(~done && k<=N)
        try
            if k==1
                otplanes{n}=[]; 
            end
            in=find(plane(:,n)==min(plane((plane(:,n)-on(k))>=0,n)));
            in2=find(plane(:,n)==max(plane((off(k)-plane(:,n))>=0,n)));
            otplanes{n}=vertcat(otplanes{n},[in;in2]);
        catch e %e is an MException struct
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);

             disp(['something happened at plane=' int2str(n) ', k trial=' int2str(k)...
                 ' and i slice =' int2str(i) ]);
             done = 1;
        end
        k = k+1;
        if rem(k,50)==0
            disp(['working on plane ' int2str(n) ', trial ' int2str(k) '/' int2str(N)]); toc;
        end
    end
end
filel= size(plane,1);