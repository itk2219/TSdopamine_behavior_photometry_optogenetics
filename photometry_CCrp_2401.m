function photometry_water_tone_puff_1811  

% Extract behavior parameters from water, tone, puff responses in Iku/Will rig
 cd('E:\2024-Hokudai\behavior\headfix\eachAnimal\IK216_CrhAi9\CCrp2');
% cd('E:\2024-Hokudai\behavior\headfix\eachAnimal\IK207_biSNL\CCrp3');
 % cd('E:\2024-Hokudai\behavior\headfix\eachAnimal\IK219_CrhAi9\CCrp1');
% cd('E:\2024-Hokudai\behavior\headfix\eachAnimal\IK220_DA2m+tdT\CCrp5');
%filename = dir('IK3482A_WM3_A6alog'); 
filename = dir('IK216_CCrp_250515');  

for ii = 1: length(filename)
    clearvars -except filename ii
    filename(ii).name
    file_ID = fopen(filename(ii).name,'r');
% read fiber photometry

CC_analog = fread(file_ID, inf, 'double', 0, 'b');

A = reshape(CC_analog, 2, [ ]);
% B = reshape (A, 2, 8, [ ]);
  B = reshape (A, 2, 10, [ ]);%CCrp
% B = reshape (A, 2, 12, [ ]);%CCrptest2

GCaMP = B (:,5,:);
GCaMP = reshape (GCaMP,[],1);

tdTom = B (:,6,:);
tdTom = reshape (tdTom,[],1);

lick = B (:,3,:);
lick = reshape (lick,[],1);

%water = B (:,4,:);
water = B (:,2,:);
water = reshape (water,[],1);
water_on = crossing(water,[],2); %threshold(mV)
water_on_ts = (water_on(1:2:end)).';
water_off_ts = (water_on(2:2:end)).';

%trial_type_sig = B (:,8,:);
trial_type_sig = B (:,4,:);
trial_type_sig = reshape (trial_type_sig,[],1);

odor2 = B (:,7,:);
odor2 = reshape (odor2,[],1);
odor2_on = crossing(odor2,[],2); %threshold(mV)
odor2_on_ts = (odor2_on(1:2:end)).';
odor2_off_ts = (odor2_on(2:2:end)).';

odor3 = B (:,8,:);
odor3 = reshape (odor3,[],1);
odor3_on = crossing(odor3,[],2); %threshold(mV)
odor3_on_ts = (odor3_on(1:2:end)).';
odor3_off_ts = (odor3_on(2:2:end)).';

odor4 = B (:,9,:);
odor4 = reshape (odor4,[],1);
odor4_on = crossing(odor4,[],2); %threshold(mV)
odor4_on_ts = (odor4_on(1:2:end)).';
odor4_off_ts = (odor4_on(2:2:end)).';

% puff = B (:,6,:);
% puff = reshape (puff,[],1);
% puff_on = crossing(puff,[],2); %threshold(mV)
% puff_on_ts = (puff_on(1:2:end)).';
% puff_off_ts = (puff_on(2:2:end)).';

% plot the gcamp and tdtom signals for the whole day %%%%%%%%%%%%%%%

figure
plot(GCaMP, 'g');
hold on;
plot(tdTom+10, 'r');
% plot(lick+15, 'k');
plot(water,'b');
plot(trial_type_sig,'m');
%plot(puff,'k')

title('photometry (g,r), puff(k), water(b)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate lick rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%min_lick_value = -1;
min_lick_value = 3; 
max_lick_value = 5; 
% max_lick_value = 1.2;
lick_inhibition_width = 50;

figure
plot(lick,'k')
hold on

filtered_lick = lowpass_signal(lick,1000,40);
plot(filtered_lick,'r')

[extrema_a,extrema_b,extrema_c,extrema_d] = extrema(filtered_lick);
filtered_lick = zeros(size(filtered_lick));
plot(extrema_d,extrema_c,'o')

for f = 1:length(extrema_d)    %extrema_d is ts for min(filtered_lick)
if extrema_c(f) > min_lick_value && extrema_c(f) < max_lick_value
    filtered_lick(extrema_d(f),1) = 1;
end
end

for f = 1:(length(filtered_lick)-100)
for i = 1:lick_inhibition_width
    if filtered_lick(f) == 1
    if filtered_lick(f+i) == 1
        filtered_lick(f+i) = 0;
    end
    end
end
end

plot(max_lick_value*filtered_lick,'b')
title('lick raw (k), filtered1(r), filtered2(b)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean photometry signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %remove 60Hz noise
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
               'DesignMethod','butter','SampleRate',1000);
filtered_GCaMP = filtfilt(d,GCaMP);
GCaMP = filtered_GCaMP;
filtered_tdTom = filtfilt(d,tdTom);
tdTom = filtered_tdTom;

%remove sudden rise noise

% diff_GCaMP = diff(GCaMP);
% std_GCaMP = std(diff_GCaMP);
% ind_noise = find(diff_GCaMP>5*std_GCaMP);
% exclude_green = size(ind_noise)
% 
% for i = 1:length(ind_noise)
% GCaMP(ind_noise(i))=GCaMP(ind_noise(i)-1);
%     for k=1:1000
%         if GCaMP(ind_noise(i)+k)-GCaMP(ind_noise(i)+k-1)>5*std_GCaMP
%             GCaMP(ind_noise(i)+k) = GCaMP(ind_noise(i)+k-1);
%         end
%     end
% end
% 
% diff_tdTom = diff(tdTom);
% std_tdTom = std(diff_tdTom);
% ind_noise = find(diff_tdTom>5*std_tdTom);
% exclude_red = size(ind_noise)
% 
% for i = 1:length(ind_noise)
% tdTom(ind_noise(i))=tdTom(ind_noise(i)-1);
%     for k=1:1000
%         if tdTom(ind_noise(i)+k)-tdTom(ind_noise(i)+k-1)>5*std_tdTom
%             tdTom(ind_noise(i)+k) = tdTom(ind_noise(i)+k-1);
%         end
%     end
% end

%smoothing
normG = smooth(GCaMP,50);
GCaMP = normG;
normR = smooth(tdTom,50);
tdTom = normR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the start of each trial and type %%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial_type_sig_on = crossing(trial_type_sig,1:length(trial_type_sig),2.5);
trial_type_sig_off = (trial_type_sig_on(2:2:end)).';
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

ind = find(trial_type_sig_on(:,1)>0);
trial_type_sig_on = trial_type_sig_on(ind,:);

TTL_end = [];
for i = 1:size(trial_type_sig_on,1)-1
    ind = find(trial_type_sig_off>trial_type_sig_on(i,1) & trial_type_sig_off<trial_type_sig_on(i+1,1));
    if length(ind)==0
        TTL_end = [TTL_end,0];
    else
        ind_max = max(ind);
    TTL_end = [TTL_end,trial_type_sig_off(ind_max)];
    end
end
TTL_end = [TTL_end,trial_type_sig_off(end)];
TTL_end = TTL_end';

for iii = 1: max(trial_type_sig_on(:,2))
    trial_type = find(trial_type_sig_on(:,2)==iii);
%     trial_type_ts{iii} = trial_type_sig_on(trial_type);  % use first TTL
    trial_type_ts{iii} = TTL_end(trial_type);  % use last TTL
end
    
% check trial_type_ts
figure
plot(trial_type_sig,'m');
hold on
for iii = 1: max(trial_type_sig_on(:,2))
    if size(trial_type_ts{iii},1)>0
    plot(trial_type_ts{iii},ones,'o')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get average F using 2-1s before trials

all_trial = trial_type_sig_on(:,1);
    ts = round(all_trial);
    ts = sort(ts);
    ind = find( ts-2000>0,1,'first');
    ind2 = find( ts-1000< length(GCaMP),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat([-2000:-1000],length(ts),1),ts);
    baseline = GCaMP(plotind); 
    F_each = mean(baseline,2);
    F_mean = mean(F_each);
    
    baseline_lick = filtered_lick(plotind); 
    F_each_lick = mean(baseline_lick,2);
    F_mean_lick = mean(F_each_lick);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make matrix of GCaMP data

% trial_number = [length(trial_type_ts{1}),length(trial_type_ts{2}),length(trial_type_ts{3}),...
    % length(trial_type_ts{4}),length(trial_type_ts{5})]
% trial_number = [length(trial_type_ts{1}),length(trial_type_ts{2}),length(trial_type_ts{3}),...
    % length(trial_type_ts{4}),length(trial_type_ts{5}),length(trial_type_ts{7}),length(trial_type_ts{9})]

%trigger = {odor2_on_ts;odor3_on_ts;odor4_on_ts};

    % trigger = {trial_type_ts{1};trial_type_ts{2};trial_type_ts{3};trial_type_ts{4};...
    %     trial_type_ts{5};trial_type_ts{6}}; %test
   % 
   % trigger = {trial_type_ts{1};trial_type_ts{7};trial_type_ts{8};trial_type_ts{3};...
   %      trial_type_ts{5}}; %test2

   trigger = {trial_type_ts{1};trial_type_ts{3};...
        trial_type_ts{5}}; %training

% trigger = {trial_type5_ts};
    
%figure
% plotColors = color_select(4);
% plotColors = {'r','m','y','g','c','k'};
% plotColors = colormap(parula(7));
% plotColors = {plotColors(7,:),plotColors(6,:),plotColors(5,:),plotColors(4,:),plotColors(3,:),plotColors(2,:),plotColors(1,:)};
% plotColors = {[0,0,1],[0,0,0.5],[0,1,0],[0,0.5,0],[0.5,0.5,0], [0.5,0.5,0.5]}; %test
% plotColors = {[0,0,1],[0,0,0.5],[0,0.5,0],[0,1,0],[0.5,0.5,0]}; %test2
  plotColors = {[0,0,1],[1,0,0],[0,1,0]}; %training
plotdata = GCaMP;
plotWin = [-2000:10000];
responseWin = [3001:5000];  %choose response time, 2000 is trigger, 1to3sec
%responseWin = [2001:7000];  %choose response time, 2000 is trigger
% legend_name ={'80rewON','80rewOff','80pufON','80pufOFF','0rew','freePuff'};%test
% legend_name = {'100rewON','80rewOon','20pufON','100pufON','0rew'};%test2
legend_name = {'100rewON','100puffON','nothing'};
M_plot = [];
DeltaF = [];
Trial_number = [];

Response = [];
ste_Response = [];
for i = 1:length(trigger)
    ts = round(trigger{i});
    if ~exist('triggerB');
        triggerB = trigger;
    end
    tsB = round(triggerB{i});

    ind = find( tsB+ plotWin(1)>0,1,'first');
    ind2 = find( ts+ plotWin(end)< length(plotdata),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),ts);
    rawTrace = plotdata(plotind);
    
    tsB = tsB(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),tsB);
    rawTraceB = plotdata(plotind);  
        
    F = mean(rawTraceB(:,1:1000),2);        %using this time window in plotwin for baseline
     deltaF = bsxfun(@minus,rawTrace,F);
%     deltaF = bsxfun(@minus,rawTrace,F)/2.5;
%     deltaF = bsxfun(@rdivide,rawTrace,F);
deltaF = 100*deltaF/F_mean;

%remove outlier  %do not use to save data
% mean_intensity = mean(deltaF(:,:)');
% mean_mean_intensity = mean(mean_intensity);
% std_intensity = std(mean_intensity);
% remove_outlier = find(mean_intensity<mean_mean_intensity+3*std_intensity & ...
%     mean_intensity>mean_mean_intensity-3*std_intensity);
% % remove_outlier = find(mean_intensity>-0.06);  
% if size(deltaF,1)>length(remove_outlier)
% outlier_green = size(deltaF,1)-length(remove_outlier)
% deltaF = deltaF(remove_outlier,:);
% end

      m_plot = mean(deltaF);
%      m_plot = mean(deltaF(end-20:end,:));
%      m_plot = deltaF(6,:);
    s_plot = std(deltaF)/sqrt(length(ts));
%     errorbar_patch(plotWin,m_plot,s_plot,plotColors(i,:));
errorbar_patch(plotWin,m_plot,s_plot,plotColors{i});
M_plot = [M_plot; m_plot];

DeltaF = [DeltaF;deltaF];
Trial_number = [Trial_number size(deltaF,1)];
response = deltaF(:,responseWin);
response = mean(response');
Response = [Response mean(response)];
ste_Response = [ste_Response std(response)/sqrt(length(response))];
end
legend(legend_name)
xlabel('time - stimulus on (ms)')
ylabel('dF/F (%)')
title('DA sensor')
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')


%figure
bar(Response)
hold on
errorbar(Response, ste_Response,'b.')
title('DA sensor')
xlabel('stimulus')
ylabel('response (0-1s)')
h=gca;
h.XTickLabel = legend_name;
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',16)
set(gcf,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raster plot green

scrsz = get(groot,'ScreenSize');
%figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(DeltaF,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(DeltaF(:,1:length_x),trialNum, binSize,[]),2));
% imagesc(binedF,[-5 5]);
imagesc(binedF,[-20 20]);
colormap yellowblue
xlabel('time - water (s)');
ylabel('trials')
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title('DA sensor')
hold on;

% 2) plot the triggers
% % mark triggerP
% for k = 1:trialNum;
% % x = [(Tmark(k)-plotWin(1))/binSize (Tmark(k)-plotWin(1))/binSize]; % odor timing
% x = [(Tmark(k)-plotWin(1)+water_delay)/binSize (Tmark(k)-plotWin(1)+water_delay)/binSize]; %water timing
% y = [k-0.5 k+0.5];
% plot(x,y,'c')
% end

% mark trigger
% x2 = [(-plotWin(1)+1000)/binSize (-plotWin(1)+1000)/binSize];%water
x2 = [-plotWin(1)/binSize -plotWin(1)/binSize]; % odor
plot(x2,[0 trialNum+0.5],'r','Linewidth',1.5)

% divide trigger
for j = 1:length(Trial_number)-1  
     plot([0 length_x/binSize],[sum(Trial_number(1:j))+0.5 sum(Trial_number(1:j))+0.5],'m','Linewidth',1.5)   
end
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make matrix of tdTom data

% figure

plotdata = filtered_lick;

DeltaF_lick = [];
Trial_number = [];

Response = [];
ste_Response = [];
for i = 1:length(trigger)
    ts = round(trigger{i});
    if ~exist('triggerB');
        triggerB = trigger;
    end
    tsB = round(triggerB{i});

    ind = find( tsB+ plotWin(1)>0,1,'first');
    ind2 = find( ts+ plotWin(end)< length(plotdata),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),ts);
    rawTrace = plotdata(plotind);
    
    tsB = tsB(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),tsB);
    rawTraceB = plotdata(plotind); 
    
    F = mean(rawTraceB(:,1:1000),2);        %using this time window in plotwin for baseline
     % deltaF = bsxfun(@minus,rawTrace,F);
deltaF = bsxfun(@minus,rawTrace,F)/2.5;
%     deltaF = bsxfun(@rdivide,rawTrace,F);
deltaF = 100*deltaF/F_mean_lick;


%remove outlier  % do not use to save data
% mean_intensity = mean(deltaF(:,:)');
% mean_mean_intensity = mean(mean_intensity);
% std_intensity = std(mean_intensity);
% remove_outlier = find(mean_intensity<mean_mean_intensity+3*std_intensity & ...
%     mean_intensity>mean_mean_intensity-3*std_intensity);
% % remove_outlier = find(mean_intensity>-0.06);
% if size(deltaF,1)>length(remove_outlier)
% outlier_red = size(deltaF,1)-length(remove_outlier)
% deltaF = deltaF(remove_outlier,:);
% end


     m_plot = mean(deltaF);
%      m_plot = mean(deltaF(end-20:end,:));
%      m_plot = deltaF(6,:);
    s_plot = std(deltaF)/sqrt(length(ts));
%     errorbar_patch(plotWin,m_plot,s_plot,plotColors(i,:));
 errorbar_patch(plotWin,m_plot,s_plot,plotColors{i});

DeltaF_lick = [DeltaF_lick;deltaF];
Trial_number = [Trial_number size(deltaF,1)];
response = deltaF(:,responseWin);
response = mean(response');
Response = [Response mean(response)];
ste_Response = [ste_Response std(response)/sqrt(length(response))];
end
legend(legend_name)
xlabel('time - stimulus on (ms)')
ylabel('dF/F (%)')
title('lick')
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')

figure
bar(Response)
hold on
errorbar(Response, ste_Response,'b.')
xlabel('stimulus')
ylabel('response (0-1s)')
h=gca;
h.XTickLabel = legend_name;
title('RCaMP')
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',16)
set(gcf,'color','w')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raster plot red

scrsz = get(groot,'ScreenSize');
% figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(DeltaF_lick,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(DeltaF_lick(:,1:length_x),trialNum, binSize,[]),2));
imagesc(binedF,[-5 5]);
colormap yellowblue
xlabel('time - water (s)');
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
h.YTickLabel = {};
title('RCaMP')
hold on;

% 2) plot the triggers
% % mark triggerP
% for k = 1:trialNum;
% % x = [(Tmark(k)-plotWin(1))/binSize (Tmark(k)-plotWin(1))/binSize]; % odor timing
% x = [(Tmark(k)-plotWin(1)+water_delay)/binSize (Tmark(k)-plotWin(1)+water_delay)/binSize]; %water timing
% y = [k-0.5 k+0.5];
% plot(x,y,'c')
% end

% mark trigger
% x2 = [(-plotWin(1)+1000)/binSize (-plotWin(1)+1000)/binSize];%water
x2 = [-plotWin(1)/binSize -plotWin(1)/binSize]; % odor
plot(x2,[0 trialNum+0.5],'r')

% divide trigger
for j = 1:length(Trial_number)-1  
     plot([0 length_x/binSize],[sum(Trial_number(1:j))+0.5 sum(Trial_number(1:j))+0.5],'m','Linewidth',1)   
end
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make matrix of licking data

figure

plotdata = filtered_lick; 

Lick_trace2 = [];
Lick_trace = [];
Response_lick = [];
ste_Response_lick = [];
Response_lick_each = [];
for i = 1:length(trigger)
    ts = round(trigger{i});

    ind = find( ts+ plotWin(1)>0,1,'first');
    ind2 = find( ts+ plotWin(end)< length(plotdata),1,'last');
    ts = ts(ind:ind2);
    plotind = bsxfun(@plus, repmat(plotWin,length(ts),1),ts);
    rawTrace = plotdata(plotind);
    rawTrace2 = smooth(rawTrace,1000);
 

     m_plot = mean(rawTrace)*1000;
     m_plot = smooth(m_plot,1000);%500
      m_plot = m_plot';
    s_plot = std(rawTrace)*100/sqrt(length(ts));
    s_plot = smooth(s_plot,1000);%500
    s_plot = s_plot';

errorbar_patch(plotWin,m_plot,s_plot,plotColors{i});

Lick_trace = [Lick_trace;rawTrace];
Lick_trace2 = [Lick_trace2;rawTrace2];


response_lick = rawTrace(:,4000:5000); %2000 is trigger, 2to3
%response_lick = rawTrace(:,2000:5000); %2000 is trigger
response_lick = mean(response_lick');
Response_lick = [Response_lick mean(response_lick)]
Response_lick2 = 1000*Response_lick
Response_lick_each = [Response_lick_each response_lick];
ste_Response_lick = [ste_Response_lick std(response_lick)/sqrt(length(response_lick))];
end
legend(legend_name)
xlabel('time - water in (ms)')
title('lick')

figure
%bar(1000*Response_lick)
bar(1000*Response_lick)
hold on
%errorbar(1000*Response_lick, 1000*ste_Response_lick,'b.')
errorbar(1000*Response_lick, 1000*ste_Response_lick,'b.')
hold on

title('lick')
ylabel('lick/s')
h=gca;
h.XTickLabel = legend_name;
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',16)
set(gcf,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raster plot (lick)

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(Lick_trace,1); binSize = 100;
length_x = plotWin(end)-plotWin(1);
binedF = squeeze(mean(reshape(Lick_trace(:,1:length_x),trialNum, binSize,[]),2));
% imagesc(binedF,[-1 1]);
imagesc(binedF,[-0.03 0.03]);
colormap yellowblue
xlabel('time - water (s)');
h=gca;
h.XTick = 0:10:(length_x/binSize);
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
title('lick')
hold on;

% 2) plot the triggers
% % mark triggerP
% for k = 1:trialNum;
% % x = [(Tmark(k)-plotWin(1))/binSize (Tmark(k)-plotWin(1))/binSize]; % odor timing
% x = [(Tmark(k)-plotWin(1)+water_delay)/binSize (Tmark(k)-plotWin(1)+water_delay)/binSize]; %water timing
% y = [k-0.5 k+0.5];
% plot(x,y,'c')
% end

% mark trigger
% x2 = [(-plotWin(1)+1000)/binSize (-plotWin(1)+1000)/binSize];%water
x2 = [-plotWin(1)/binSize -plotWin(1)/binSize]; % odor
plot(x2,[0 trialNum+0.5],'g')
x2 = [(-plotWin(1)+3000)/binSize (-plotWin(1)+3000)/binSize]; % water
plot(x2,[0 trialNum+0.5],'b')

% divide trigger
for j = 1:length(Trial_number)-1  
     plot([0 length_x/binSize],[sum(Trial_number(1:j))+0.5 sum(Trial_number(1:j))+0.5],'m','Linewidth',1)   
end
box off
set(gca,'tickdir','out')
set(gca,'TickLength',2*(get(gca,'TickLength')))
set(gca,'FontSize',20)
set(gcf,'color','w')

save_name = strcat(filename(ii).name,'_lick_aligned');
save (save_name,'DeltaF','DeltaF_lick','Trial_number','Lick_trace','legend_name','Response_lick','Response_lick_each','rawTrace','Lick_trace2')

end
end
