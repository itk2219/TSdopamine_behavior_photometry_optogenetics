function analyze_blink_1810

% this function analyze blink in Iku/Will's rig
% First, run 'detect_video_LED_1808_Iku' (this takes time)
% Next, run 'videosync_every10_1808_Iku'
% Also run 'detect_blink_1810' (this takes time)

% c{1} = 'E:\2024-Hokudai\behavior\headfix\eachAnimal\IK149_biSNL\CCrp1';
 c{1} = 'E:\2024-Hokudai\behavior\headfix\eachAnimal\IK216_CrhAi9\CCrp2';
% c{1} = 'E:\2024-Hokudai\behavior\headfix\eachAnimal\IK220_DA2m+tdT\CCrp5';
 % c{1} = 'E:\2024-Hokudai\behavior\headfix\eachAnimal\IK218_CrhAi9\CCrp1';
% c{1} = 'D:\2024-Hokudai\behavior\headfix\eachAnimal\IK134_DA2m+tdT_VSDLSTS\CCrptest';
% c{1} = 'E:\2023-Keio\behavior\headfix\eachAnimal\IK96_CRF+tdT_VS-TS-VTA-SNL\videos\CCrp3';
% c{1} = 'E:\2023-Keio\behavior\headfix\eachAnimal\IK97_CRF+tdT_VS-TS-VTA-SNL\videos\CCrp5';

cd(c{1});
% foldername = dir('*-09-14_*');
% cd(foldername.name);

% read fiber photometry
file = strcat('IK216_CCrp_250515');
file_ID = fopen(file, 'r');
CC_cart_analog = fread(file_ID, inf, 'double', 0, 'b');

A = reshape(CC_cart_analog, 2, [ ]);
 % B = reshape (A, 2, 8, [ ]);
 B = reshape (A, 2, 10, [ ]);

trial_type_sig = B (:,4,:);
trial_type_sig = reshape (trial_type_sig,[],1);

odor3 = B (:,8,:);
odor3 = reshape (odor3,[],1);
odor3_on = crossing(odor3,[],2); %threshold(mV)
odor3_on_ts = (odor3_on(1:2:end)).';
odor3_off_ts = (odor3_on(2:2:end)).';

puff = B (:,10,:); %free water or puff
puff = reshape (puff,[],1);
puff_on = crossing(puff,[],2); %threshold(mV)
puff_on_ts = (puff_on(1:2:end)).';
puff_off_ts = (puff_on(2:2:end)).';

odor2 = B (:,7,:);
odor2 = reshape (odor2,[],1);
odor2_on = crossing(odor2,[],2); %threshold(mV)
odor2_on_ts = (odor2_on(1:2:end)).';
odor2_off_ts = (odor2_on(2:2:end)).';

odor4 = B (:,9,:);%bigstim
odor4 = reshape (odor4,[],1);
odor4_on = crossing(odor4,[],2); %threshold(mV)
odor4_on_ts = (odor4_on(1:2:end)).';
odor4_off_ts = (odor4_on(2:2:end)).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the start of each trial and type %%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial_type_sig_on = crossing(trial_type_sig,1:length(trial_type_sig),2.5);
% trial_type_sig_off = (trial_type_sig_on(2:2:end)).';
trial_type_sig_on = (trial_type_sig_on(1:2:end)).';


trial_type_sig_on(:,2) = trial_type_sig_on(:,1);
trial_type_sig_on(:,2) = 0;

for i = 1:size(trial_type_sig_on,1)
    
temp_B = 1;

temp_A = trial_type_sig_on(i,1);    

for f = i+1:size(trial_type_sig_on,1)
if trial_type_sig_on(f,1) < temp_A + 1000  %identify multiple signals within 1s

temp_B = temp_B + 1;  
  
end
end

trial_type_sig_on(i,2) = temp_B;  %number of signals within 1s after trial_type_sig_on(i,1)

end

% delete duplicates of the same trial

for i = 1:size(trial_type_sig_on,1)
    
temp_B = 1;

temp_A = trial_type_sig_on(i,1);    

for f = i+1:size(trial_type_sig_on,1)
if trial_type_sig_on(f,1) < temp_A + 1000

trial_type_sig_on(f,:) = 0;
  
end
end

end

odor2_water_ts = [];odor2_omission_ts=[];odor3_ts = [];odor4_puff_ts = [];odor4_omission_ts = [];

trial_type1 = find(trial_type_sig_on(:,2)==1);
trial_type1_ts = trial_type_sig_on(trial_type1);
for i = 1:length(trial_type1)
    odor2_water = find(odor2_on_ts > trial_type1_ts(i) & odor2_on_ts < trial_type1_ts(i)+1000);
    odor2_water_ts = [odor2_water_ts; odor2_on_ts(odor2_water)];
end

trial_type2 = find(trial_type_sig_on(:,2)==2);
trial_type2_ts = trial_type_sig_on(trial_type2);
for i = 1:length(trial_type2)
    odor2_omission = find(odor2_on_ts > trial_type2_ts(i) & odor2_on_ts < trial_type2_ts(i)+1000);
    odor2_omission_ts = [odor2_omission_ts; odor2_on_ts(odor2_omission)];
end

trial_type4 = find(trial_type_sig_on(:,2)==4);
trial_type4_ts = trial_type_sig_on(trial_type4);
for i = 1:length(trial_type4)
    odor4_puff = find(odor4_on_ts > trial_type4_ts(i) & odor4_on_ts < trial_type4_ts(i)+1000);
    odor4_puff_ts = [odor4_puff_ts; odor4_on_ts(odor4_puff)];
end

trial_type11 = find(trial_type_sig_on(:,2)==11);
trial_type11_ts = trial_type_sig_on(trial_type11);
for i = 1:length(trial_type11)
    odor4_omission = find(odor4_on_ts > trial_type11_ts(i) & odor4_on_ts < trial_type11_ts(i)+1000);
    odor4_omission_ts = [odor4_omission_ts; odor4_on_ts(odor4_omission)];
end

%% analyze eye area

% calculate trigger frame
% Trial type 1: odor2 water on
% Trial type 2: odor2 water off
% Trial type 3: odor3 water off
% Trial type 4: odor4 puff on
% Trial type 11: odor4 puff off

load ('video_time') % 'video_t'
frame_interval = median(diff(video_t)) %should be ~33.3 ms
figure
hist(diff(video_t))
title('frame interval')
xlabel('ms')

puff_frame = [];
for i = 1:length(puff_on_ts)
    if puff_on_ts(i)<video_t(end)
        puff_frame1 = find(video_t < puff_on_ts(i),1,'last');
        puff_frame = [puff_frame; puff_frame1];
    end
end

odor_frame = [];
odor_ts = [odor2_on_ts;odor3_on_ts;odor4_on_ts];
odor_ts = sort(odor_ts);
for i = 1:length(odor_ts)
    if odor_ts(i)<video_t(end)
        odor_frame1 = find(video_t < odor_ts(i),1,'last');
        odor_frame = [odor_frame; odor_frame1];
    end
end

odor2_frame = [];
for i = 1:length(odor2_on_ts)
%for i = 1:length(odor2_water_ts)
    if odor2_on_ts(i)<video_t(end)
    %if odor2_water_ts(i)<video_t(end)
        odor2_frame1 = find(video_t < odor2_on_ts(i),1,'last');
        %odor2_frame1 = find(video_t < odor2_water_ts(i),1,'last');
        odor2_frame = [odor2_frame; odor2_frame1];
    end
end

odor3_frame = [];
for i = 1:length(odor3_on_ts)
    if odor3_on_ts(i)<video_t(end)
        odor3_frame1 = find(video_t < odor3_on_ts(i),1,'last');
        odor3_frame = [odor3_frame; odor3_frame1];
    end
end

odor4_frame = [];
for i = 1:length(odor4_on_ts)
    if odor4_on_ts(i)<video_t(end)
        odor4_frame1 = find(video_t < odor4_on_ts(i),1,'last');
        odor4_frame = [odor4_frame; odor4_frame1];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get average eye area using 2-1s before trials

load ('eye') %'eye_Area','eye_Long','eye_Short'

eye_Area = smooth(eye_Area,10);
    ts = odor_frame;
    %ind = find( ts-168>0,1,'first');
    ind = find( ts-60>0,1,'first');
    %ind2 = find( ts-84< length(eye_Area),1,'last');
    ind2 = find( ts-30< length(eye_Area),1,'last');
    ts = ts(ind:ind2);
    % plotind = bsxfun(@plus, repmat([-168:-84],length(ts),1),ts);
    plotind = bsxfun(@plus, repmat([-60:-30],length(ts),1),ts);
    baseline = eye_Area(plotind); 
    F_each = mean(baseline,2);
    F_mean = mean(F_each);

    max_eyeArea = max(F_each);  % can use this instead of average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% average plot

% trigger = {puff_frame};
% trigger = {odor2_frame; odor3_frame; odor4_frame; puff_frame};%test
  trigger = {odor2_frame; odor3_frame; odor4_frame};%training
% triggerP = {odor2_frame+168; odor3_frame+168; odor4_puff_frame};

trigger_n = [length(odor2_frame),length(odor3_frame),length(odor4_frame)]

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
% plotColors = {'b','r','g', 'k'};%test
 plotColors = {'b','r','g'};%training
plotdata = eye_Area;
%plotWin = [-168:336];  %84 frame/second
plotWin = [-60:180];  %30 frame/second
%responseWin = [314:334];  %choose response time, 168 is trigger
responseWin = [120:150];  %choose response time, 60 is trigger, 1to3sec, CS
% responseWin = [150:210];  %choose response time, 60 is trigger, 3to5sec, US
% legend_name = {'80%water','80%puff','nothing', 'freeP'};%test
legend_name = {'100%water','100%puff','nothing'};%training
% legend_name = {'bigStim','bigWater'};


DeltaF = [];
Tmark = [];
Trial_number = [];
response_all = [];
Response = [];
ste_Response = [];
for i = 1:length(trigger)
%     ts = round(trigger{i});
ts = trigger{i};
    if ~exist('triggerB');
        triggerB = trigger;
    end
%     tsB = round(triggerB{i});
tsB = triggerB{i};

    ind = find( tsB+ plotWin(1)>0,1,'first');
    ind2 = find( ts+ plotWin(end)< length(plotdata),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),ts);
    rawTrace = plotdata(plotind);
    
    tsB = tsB(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),tsB);
    rawTraceB = plotdata(plotind); 
    
    %F = mean(rawTraceB(:,85:150),2);
    F = mean(rawTraceB(:,31:50),2);
%     deltaF = bsxfun(@minus,rawTrace,F)/2.5;
deltaF = bsxfun(@minus,rawTrace,F);
deltaF = deltaF/F_mean;
    %deltaF_F = bsxfun(@rdivide,rawTrace,F);
     m_plot = mean(deltaF);
%      m_plot = mean(deltaF(end-20:end,:));
%      m_plot = deltaF(6,:);
    s_plot = std(deltaF)/sqrt(length(ts));
errorbar_patch(plotWin,m_plot,s_plot,plotColors{i});

response = deltaF(:,responseWin);
response = mean(response');
response_all = [response_all,response];
Response = [Response mean(response)]
ste_Response = [ste_Response std(response)/sqrt(length(response))];
 if ~exist('triggerP')
        triggerP = triggerB;
    end
tsP = round(triggerP{i});
tsP = tsP(ind:ind2);
DeltaF = [DeltaF;deltaF];
Tmark = [Tmark;tsP-ts];
Trial_number = [Trial_number size(deltaF,1)];
end


legend(legend_name)
xlabel('time - event (s)')
ylabel('fraction change of eye area')
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')
h=gca;
%h.XTick = -168:84:336;
% h.XTick = -60:30:120;
% h.XTickLabel = {-2:1:4};
h.XTick = -60:30:180;
h.XTickLabel = {-2:1:6};


figure
bar(Response)
hold on
errorbar(Response, ste_Response,'b.')
hold on
title('eye area')
ylabel('responses (0-3s)')
h=gca;
h.XTickLabel = legend_name;

%% raster plot

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(DeltaF,1); binSize = 6;%7
binedF = squeeze(mean(reshape(DeltaF(:,1:(plotWin(end)-plotWin(1))),trialNum, binSize,[]),2));
%imagesc(binedF,[-1 1]);
% imagesc(binedF,[-0.5 0.5]);
 imagesc(binedF,[-1 1]);
colormap yellowblue
xlabel('time - event (s)');
h=gca;
%h.XTick = 0.5:84/binSize:5*84/binSize;
h.XTick = 0.5:30/binSize:5*30/binSize;
h.XTickLabel = {-2:6};
hold on;
% 2) plot the triggers
for k = 1:trialNum
% x = [(Tmark(k)-plotWin(1))/binSize (Tmark(k)-plotWin(1))/binSize]; % odor timing
%x = [(Tmark(k)-plotWin(1)+168)/binSize+0.5 (Tmark(k)-plotWin(1)+168)/binSize+0.5]; %water timing
% x = [(Tmark(k)-plotWin(1)+60)/binSize+0.5 (Tmark(k)-plotWin(1)+60)/binSize+0.5]; %water timing
% y = [k-0.5 k+0.5];
% plot(x,y,'r')
end
% x2 = [(-plotWin(1)+1000)/binSize (-plotWin(1)+1000)/binSize];%water
x2 = [-plotWin(1)/binSize+0.5 -plotWin(1)/binSize+0.5]; % odor
plot(x2,[0 trialNum+0.5],'c')

% 3) divide triggers
for j = 1:length(Trial_number)-1
     plot([0 (plotWin(end)-plotWin(1))/binSize],[sum(Trial_number(1:j))+0.5 sum(Trial_number(1:j))+0.5],'m','Linewidth',1)   

 save_name = strcat('blink_aligned');
 save (save_name,'DeltaF','Trial_number')

end
end

