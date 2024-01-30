function [output]=coarseRetinopathy_HW_looper
%HW
%8.5.18
%modified coarse retinopathy so can divide up the screen into as many bars
%as you want and have an uneven number of x and y bars

%%MUST modified drawstimulus code so taht it only searces of ";" as the end
%%of formulas to work (before was not properly chunking the strings of
%%formulas beacuse you could also you ',' as an end of formula.


%Formulas for looper to divide up screen into 10x5 (or x by y bars instead
%of same number for all sides_
%%

xW =1280; 
yW = 720; 

yN=6;
xN=10;

An=0:(yN+xN-1); %0:14; %a<=4 are the horizontal bars, a>4 are the  vertical bars
figure(3); clf; 
h1 = axes;
hold on;
for i=1:numel(An) 
a=An(i);

v=min(max(a-(yN-1),0),1); %v=0 if a <=4, v=1 if a>4
w=max(min((yN)-a,1),0); %w=1 if a <=4, w=0 if a>4
 

x_pos=v*((a-yN)*xW/xN + xW/(xN*2)) + w*xW/2;

y_pos = v*yW/2 + w*(a*yW/yN+ yW/(yN*2));

xP=v*(xW/xN)+w*xW; %0.6 scaling bc 40in monitor/2PTS stimulus code converts xpixels to degrees weird when stimulus is big?
yP=v*yW+ w*(yW/yN);

output(i,:)=[x_pos y_pos xP yP v w];
if v==1;
    plot(x_pos,y_pos,'bx');
else
    plot(x_pos,y_pos,'o');
end

end
xlim([0 xW]);
ylim([0 yW]);
set(h1, 'Ydir', 'reverse')

% %%
% %

% screenXcm = 59.69;
% screenYcm = 33.65;
% screenDist = 16;
% pcmX = xW/screenXcm; 
% pcmY = yW/screenYcm;  
% % 
% x_size = 2*180/pi*atan(xP/(2*screenDist*pcmX))
% y_size = 2*180/pi*atan(yP/(2*screenDist*pcmY))
% 



end