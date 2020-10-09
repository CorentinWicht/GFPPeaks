function slider_plot(f,EEGTEMP,MaxPos,Chanlocs,Param,TEMPGFP)
% Function for generation of topoplots using a slidebar

% Corentin Wicht
% 05.06.2020

%% MAIN CODE

% FIGURE 
% f = figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(2,3,[3 6]); 
title(sprintf('Topography at %d TF',MaxPos));            
topoplotIndie(EEGTEMP(MaxPos,:),Chanlocs); % build topoplot
bgcolor = f.Color; % background color

% SLIDER
PosSlide = [0.64 0.1 0.30 0.03];
S.aSlider  = uicontrol('Parent',f,'Style','slider','unit','normalized','Position',PosSlide,...
              'value',MaxPos, 'min',0, 'max',size(EEGTEMP,1),'SliderStep',[1/size(EEGTEMP,1) 10/size(EEGTEMP,1)],...
              'callback', {@SliderCB,EEGTEMP,Chanlocs,Param,TEMPGFP});
          
% LABELS
bl1 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.63 0.095 0.01 0.03],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.94 0.095 0.02 0.03],...
                'String',num2str(size(EEGTEMP,1)),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.75 0.05 0.1 0.03],...
                'String','Select a Time Frame (TF)','BackgroundColor',bgcolor);
end

% Callback for all sliders defined above
function SliderCB(aSlider, EventData,EEGTEMP,Chanlocs,Param,TEMPGFP)
    Value = round(get(aSlider, 'Value'));  % Get the slidebar value
    assignin('base','SlideBarValue',Value);
    
    % 1. Electrodes level data
    subplot(2,3,[1 2]); 
    plot(EEGTEMP,'Color','k'); 
    title('All electrodes');
    xlabel('Time (ms)'), ylabel('uV'); axis tight; Yl=ylim(gca);  
    set(gca,'XTick',Param.Ticks,'XTickLabel',num2cell(Param.XTStr)); 
    PosinMS = round(Value/Param.SamplingRate*1000);
    LabelValue = [num2str(PosinMS + Param.Epoch(1)) 'ms ' '[' num2str(Value) 'TF' ']']; 
    xline(Value,'-',LabelValue,'Color','b','LabelVerticalAlignment','bottom');
    PosRectinTF = [Param.TEMPCompinTF(1) Yl(1) Param.TEMPCompinTF(2)-Param.TEMPCompinTF(1) Yl(2)*0.1]; 
    rectangle('Position',PosRectinTF,'FaceColor',[0 0 1 0.2],'EdgeColor',[0 0 1 0.2]);
    
    % 2. GFP
    subplot(2,3,[4 5]); plot(TEMPGFP, 'Color','k'); 
    title('GFP'); % axis tight;
    axis tight; Yl=ylim(gca); % retrieve auto y-limits
    xlabel('Time (ms)'), ylabel('GFP');
    set(gca,'XTick',Param.Ticks,'XTickLabel',num2cell(Param.XTStr));
    xline(Value,'-',LabelValue,'Color','b','LabelVerticalAlignment','bottom');
    PosRectinTF = [Param.TEMPCompinTF(1) Yl(1) Param.TEMPCompinTF(2)-Param.TEMPCompinTF(1) Yl(2)*0.1]; 
    rectangle('Position',PosRectinTF,'FaceColor',[0 0 1 0.2],'EdgeColor',[0 0 1 0.2]);
    LegendLabel = {'GFP';Param.CompN};
    Lgd = legend(LegendLabel); title(Lgd,'Component','Color','k');

    % 3. TOPOPLOTS
    subplot(2,3,[3 6]); 
    title(sprintf('Topography at %d TF',Value));            
    topoplotIndie(EEGTEMP(Value,:),Chanlocs);
end
