function playcohmotion_HW

global Mstate screenPTR screenNum loopTrial daq

%global daq  %Created in makeAngleTexture

global Stxtr %Created in makeSyncTexture

global DotFrame %created in makeRandomDots

global DotColor %created in makerandom dots HW

%% HW from flash cartesia

% doSession = daq.createSession('ni');
% doSession.addDigitalChannel('Dev1', 'Port0/Line0', 'OutputOnly');
% doSession.Rate = 20000;
% 
% doSession_stim = daq.createSession('ni'); %sends pulse for stim only?
% doSession_stim.addDigitalChannel('Dev1', 'Port0/Line1', 'OutputOnly');
% doSession_stim.Rate = 20000; %25000;
% %Wake up the daq:
% %DaqDOut(daq, 0, 0); %I do this at the beginning because it improves timing on the call to daq below
% 
% aoSession = daq.createSession('ni');
% aoSession.addAnalogOutputChannel('Dev1','ao0', 'Voltage');
% aoSession.Rate = 20000; %25000;
% doSession.Rate = 20000; % 25000;
% doSession_stim.Rate = 20000; %25000;

%% pulses for opto

% Pstruct = getParamStruct;
% if Pstruct.light_bit ~= 0 %if 0 no light, if not 0, value is the amplitude of light signal (set to -1 if using TTL not analog/mod)
%     SampleRate=20000; %25000;
%     data=makelightpulse_HW(Pstruct, SampleRate);
%     
%     % Stuff from before? HW- I can't remember if I commented this out, or
%     % Anupam did long ago.
%     % daqregister('mcc');
%     % dio = digitalio('mcc', 0);
%     % addline(dio, 0:7, 0, 'Out');
%     % 
%     % ao = analogoutput('mcc');
%     % chan = addchannel(ao, 0);
%     % 
%     % set(ao, 'SampleRate', SampleRate);
%     % set(ao, 'TriggerType', 'Manual');
% 
%     % t = 0 : 1/SampleRate : duration;
%     % d = PulseDuration / 2 : 1/freq : duration; %frequency of pulse (ie this is one if you just want 1 rect during stim, but a vector in fraction steps for trains
% 
%     %data = pulstran(t, d, 'rectpuls', PulseDuration)' *10;
%     % data = pulstran(t, d, 'tripuls', PulseDuration, 1)' * 10;
%     %data(length(data)) = 0;
% 
%     queueOutputData(aoSession, data)
%     prepare(aoSession) %NOTE: aoSession/data will "play" during prestim or on stim depending on which light style timing it is, CHECK below.
% end



%%
Pstruct = getParamStruct;


[screenRes, pixpercmX, pixpercmY, xN, yN]=StimSizeinPixel(Pstruct.x_size, Pstruct.y_size);

% Synchronization Stimulus size; 
syncWX = round(pixpercmX*Mstate.syncSize);
syncWY = round(pixpercmY*Mstate.syncSize);

%in all of the code, we treat the screen as if it is round. this means that
%a stimulus of size x deg ends up having a size in cm of 2pi/360*x deg*monitor
%distance (this is simply the formula for the length of an arc); then
%transform from cm to pixels

sizeDotsCm=Pstruct.sizeDots*2*pi/360*Mstate.screenDist;
sizeDotsPx=round(sizeDotsCm*pixpercmX);





%%%%%%%%%%%%%%%%%%


Npreframes = ceil(Pstruct.predelay*screenRes.hz);
Nstimframes = ceil(Pstruct.stim_time*screenRes.hz);
Npostframes = ceil(Pstruct.postdelay*screenRes.hz);


%%%%%%%%%%%%%%%

SyncLoc = [0 0 syncWX-1 syncWY-1]';
SyncPiece = [0 0 syncWX-1 syncWY-1]';


%%%%%%%%%%%%%%%

Screen(screenPTR, 'FillRect', [Pstruct.backgroundR, Pstruct.backgroundG, Pstruct.backgroundB]);   % PL170812

if Pstruct.contrast==0
    r=Pstruct.backgroundR;
    g=Pstruct.backgroundG;
    b=Pstruct.backgroundB;
    if Pstruct.dotColor==2 %multiple black and white dots
       for td= 1:size(DotColor,2)
            DotColor{td}=[r;g;b];
       end
    end
else
    r=Pstruct.redgun;
    g=Pstruct.greengun;
    b=Pstruct.bluegun;
end

 %Wake up the daq:
DaqDOut(daq, 0, 0); % I do this at the beginning because it improves timing on the first call to daq below
                    % Peichao: DaqDOut(DeviceIndex,port,data), The index of USB-1208FS on our Win7 is 6
                    % "port" 0 = port A, 1 = port B
Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip'); %make sure this starts down
%% %Play predelay %%%%
%2 PHOTON SPECIFIC EDIT. Half the pre-time = gray screen

%START PREDELAY 1 = GRAY SCREEN
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip'); %1 is up

for i = 2:Npreframes/2
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
end

%start the dots. PREDEALY PART 2
if Pstruct.dotColor==2
Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, DotColor{1},...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
else
    Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, [r g b],...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip'); %1 is up

%Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);  
%Screen(screenPTR, 'Flip'); %2 is down for pd sq?


% if strcmp(Pstruct.light_type,'pre') || strcmp(Pstruct.light_type, 'prepost')
%     if Pstruct.light_bit ~= 0
%         startBackground(aoSession);
%     end
% end
%if loopTrial ~= -1
    digWord = 7;  %Make 1st bit high
    DaqDOut(daq, 0, digWord);
%end
% outputSingleScan(doSession, 1);

for i = (Npreframes/2)+1:Npreframes %2:Npreframes
if Pstruct.dotColor==2
    Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, DotColor{1},...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);    
else
        Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, [r g b],...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
end
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
end


%%%%%Play whats in the buffer (the stimulus)%%%%%%%%%%
%DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low

% outputSingleScan(doSession_stim, 1);

% if strcmp(Pstruct.light_type,'stim') || strcmp(Pstruct.light_type,'middle')
%     if Pstruct.light_bit ~= 0
%     	startBackground(aoSession);
%     end
% end

if Pstruct.dotColor==2
Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, DotColor{1},...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
else
    Screen('DrawDots', screenPTR, DotFrame{1}, sizeDotsPx, [r g b],...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip'); %1 is up
%Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);  
%Screen(screenPTR, 'Flip'); %2 is down for pd sq?

%Screen('DrawTextures', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
%Screen(screenPTR, 'Flip');
% if loopTrial ~= -1
%     digWord = 3;  %toggle 2nd bit to signal stim on
%     DaqDOut(daq, 0, digWord);
% end
for i = 2:Nstimframes
    if Pstruct.dotColor==2
    Screen('DrawDots', screenPTR, DotFrame{i}, sizeDotsPx, DotColor{i},...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
    else
            Screen('DrawDots', screenPTR, DotFrame{i}, sizeDotsPx, [r g b],...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
    end
    Screen('DrawTextures', screenPTR,Stxtr(2),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
end
% if loopTrial ~= -1
%     digWord = 1;  %toggle 2nd bit to signal stim off
%     DaqDOut(daq, 0, digWord);
% end
% outputSingleScan(doSession_stim, 0);
%Screen(screenPTR, 'Flip');  %Show sync on last frame of stimulus
%digWord = bitxor(digWord,7);  %toggle only the 3rd bit on each grating cycle
%DaqDOut(daq, 0, digWord); 
%% %Play postdelay %%%%

%Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);  
%Screen(screenPTR, 'Flip'); %2 is down for pd sq?
if Pstruct.dotColor==2
Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, DotColor{Nstimframes},...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
else
    Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, [r g b],...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip'); %1 is up
for i = 2:Npostframes-1
    if Pstruct.dotColor==2
   Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, DotColor{Nstimframes},...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
    else
        Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, [r g b],...
        [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
    end
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
end
if Pstruct.dotColor==2
Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, DotColor{Nstimframes},...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
else
    Screen('DrawDots', screenPTR, DotFrame{Nstimframes}, sizeDotsPx, [r g b],...
    [Pstruct.x_pos Pstruct.y_pos],Pstruct.dotType);
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
Screen(screenPTR, 'Flip');
% if loopTrial ~= -1
%     %digWord = bitxor(digWord,7); %toggle all 3 bits (1st/2nd bits go low, 3rd bit is flipped)
%     %DaqDOut(daq, 0,digWord);  
%     DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low
% end
%     outputSingleScan(doSession, 0);
%Last frame of post-stim Peichao
digWord = bitxor(digWord,7); %toggle all 3 bits (1st/2nd bits go low, 3rd bit is flipped)
DaqDOut(daq, 0,digWord);

%if loopTrial ~= -1 %if you hit "run", send pulse
DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);  
Screen(screenPTR, 'Flip');


