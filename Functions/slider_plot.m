function slider_plot(f,EEGTEMP,MaxPos,Chanlocs)
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
              'value',MaxPos, 'min',0, 'max',size(EEGTEMP,1),'callback', {@SliderCB,EEGTEMP,Chanlocs});
          
% LABELS
bl1 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.63 0.095 0.01 0.03],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.94 0.095 0.02 0.03],...
                'String',num2str(size(EEGTEMP,1)),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','unit','normalized','Position',[0.75 0.05 0.1 0.03],...
                'String','Select a Time Frame (TF)','BackgroundColor',bgcolor);
end

% Callback for all sliders defined above
function SliderCB(aSlider, EventData,EEGTEMP,Chanlocs)
    Value = round(get(aSlider, 'Value'));  % Get the slidebar value
    assignin('base','SlideBarValue',Value);
    subplot(2,3,[3 6]); 
    title(sprintf('Topography at %d TF',Value));            
    topoplotIndie(EEGTEMP(Value,:),Chanlocs);
end
