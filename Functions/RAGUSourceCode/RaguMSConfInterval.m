function RaguMSConfInterval(out)

% Copyright 2009-2011 Thomas Koenig
% distributed under the terms of the GNU AFFERO General Public License

[FileName,PathName] = uiputfile('*.txt','Enter file name');

[fid,err] = fopen(fullfile(PathName,FileName),'wt');

if fid == -1
    error(err);
end

FeatLabels = {'Onset', 'Offset','Duration','Area under curve','Center of gravity','Mean GFP'};

[bootstrapsamples,labels] = ComputeConfidenceIntervals(out);

bootstrapsamples(1:5,:,:,:) = bootstrapsamples(1:5,:,:,:) * out.DeltaX;
bootstrapsamples([1,2,5],:,:,:) = bootstrapsamples([1,2,5],:,:,:) + out.TimeOnset;


bootstrapmean = mean(bootstrapsamples,2);
BootstrapLowerCI = prctile(bootstrapsamples,  out.Threshold / 2,2);
BootstrapUpperCI = prctile(bootstrapsamples,1-out.Threshold / 2,2);

for f = 1:6
    fprintf(fid,'Class');
    for cnd = 1:numel(labels);
        fprintf(fid,'\t%s: %s(-CI)\t%s: %s(+CI)',FeatLabels{f},labels{cnd},FeatLabels{f},labels{cnd});
    end
end


for c = 1:size(bootstrapsamples,4);
    fprintf(fid,'\n%i',c);
    for f = 1:6
        for cnd = 1:numel(labels);
            fprintf(fid,'\t%f\t%f',BootstrapLowerCI(f,1,cnd,c),BootstrapUpperCI(f,1,cnd,c));
        end
    end
end

save dummy.mat

fclose(fid);

msgbox('Done!')



function [BootFeatAll,BootLabelAll] = ComputeConfidenceIntervals(out)

if out.ContBetween == true
    error('Cannot compute microstate confidence intervals for continuous predictors');
end
Group = out.IndFeature;
GroupLabels = out.GroupLabels;

GroupIndices = unique(Group);
GroupIndices(isnan(GroupIndices)) = [];
SelDesign = out.Design;
SelDesign(isnan(SelDesign(:,1)),:) = [];

DoF1    = numel(unique(SelDesign(:,1)))> 1;
DoF2 = out.TwoFactors;

if (numel(unique(Group(~isnan(Group)))) > 1)
    DoGroup = 1;
else
    DoGroup = 0;
end

nBins = 1;
GMIndexS = nBins;
GMIndexE = nBins;

GIndexS = nBins+1; 
if DoGroup
    nBins = nBins + numel(GroupIndices);
end
GIndexE = nBins;

F1IndexS = nBins + 1;
if DoF1
    nBins = nBins + numel(unique(SelDesign(:,1)));
end
F1IndexE = nBins;

F1GIndexS = nBins + 1;
if DoGroup && DoF1
    nBins = nBins + numel(GroupIndices) * numel(unique(SelDesign(:,1)));
end
F1GIndexE = nBins;

F2IndexS = nBins + 1;
if DoF2
    nBins = nBins + numel(unique(SelDesign(:,2)));
end
F2IndexE = nBins;

F2GIndexS = nBins + 1;
if DoGroup && DoF2
    nBins = nBins + numel(GroupIndices) * numel(unique(SelDesign(:,2)));
end
F2GIndexE = nBins;

F1F2IndexS = nBins + 1;
if DoF2 && DoF1
    nBins = nBins + numel(unique(SelDesign(:,1))) * numel(unique(SelDesign(:,2)));
end
F1F2IndexE = nBins;

F1F2GIndexS = nBins + 1;
if DoF2 && DoF1 && DoGroup
    nBins = nBins + numel(unique(SelDesign(:,1))) * numel(unique(SelDesign(:,2))) * numel(GroupIndices);
end
F1F2GIndexE = nBins;

BootFeatAll  = nan(6,out.Iterations,nBins,size(out.MSMaps,1));
BootLabelAll = cell(nBins,1);


AllC = ones(size(out.Design,1),1);
AllC(isnan(out.Design(:,1)),1) = nan;
AllSubjectsIn = find(~isnan(Group));
nAllSubjectsIn = numel(AllSubjectsIn);

tic

h = waitbar(0,'Computing microstate stats, please wait...');

set(h,'WindowStyle','modal');

for b = 1:out.Iterations
    AllGRandIndex = ceil(rand(nAllSubjectsIn,1) * nAllSubjectsIn);
    GroupedRandIndex = AllSubjectsIn;
    for g = 1:numel(GroupIndices)
        SubjectsIn = find(Group == GroupIndices(g));
        nSubjectsIn = numel(SubjectsIn);
        GroupRandIndex = ceil(rand(nSubjectsIn,1)*nSubjectsIn);
        GroupedRandIndex(SubjectsIn) = GroupedRandIndex(SubjectsIn(GroupRandIndex));
    end
    
    
    % The overall stuff
    [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,AllGRandIndex,~isnan(Group),AllC,out.DLabels1,out.DLabels2,{'All subjects'},out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

    BootFeatAll(:,b,GMIndexS:GMIndexE,:) = feat;
    BootLabelAll(GMIndexS:GMIndexE) = lbl;
    
    % The Group stuff
    if DoGroup
        [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,GroupedRandIndex,Group,AllC,out.DLabels1,out.DLabels2,GroupLabels,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);
        
        BootFeatAll(:,b,GIndexS:GIndexE,:) = feat;
        BootLabelAll(GIndexS:GIndexE) = lbl;

    end
    
    % The F1Stuff
    if DoF1
        [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,AllGRandIndex,~isnan(Group),out.Design(:,1),out.DLabels1,out.DLabels2,{'All subjects'},out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

        BootFeatAll(:,b,F1IndexS:F1IndexE,:) = feat;
        BootLabelAll(F1IndexS:F1IndexE) = lbl;
    end

    % The Group * F1 stuff
    if DoGroup && DoF1
        [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,GroupedRandIndex,Group,out.Design(:,1),out.DLabels1,out.DLabels2,GroupLabels,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

        BootFeatAll(:,b,F1GIndexS:F1GIndexE,:) = feat;
        BootLabelAll(F1GIndexS:F1GIndexE) = lbl;
    end

    if DoF2
        [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,AllGRandIndex,~isnan(Group),out.Design(:,2),out.DLabels1,out.DLabels2,{'All subjects'},out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

        BootFeatAll(:,b,F2IndexS:F2IndexE,:) = feat;
        BootLabelAll(F2IndexS:F2IndexE) = lbl;
    end
    
    % The Group * F2 stuff
    if DoGroup && DoF2
        [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,GroupedRandIndex,Group,out.Design(:,2),out.DLabels1,out.DLabels2,GroupLabels,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);

        BootFeatAll(:,b,F2GIndexS:F2GIndexE,:) = feat;
        BootLabelAll(F2GIndexS:F2GIndexE) = lbl;
    end
    
    % The F1 * F2 stuff
    if DoF2 && DoF1
       [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,AllGRandIndex,~isnan(Group),out.Design,out.DLabels1,out.DLabels2,{'All subjects'},out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);
       
        BootFeatAll(:,b,F1F2IndexS:F1F2IndexE,:) = feat;
        BootLabelAll(F1F2IndexS:F1F2IndexE) = lbl;
    end
    % The F1 * F2 * Group stuff
    if DoF2 && DoF1 && DoGroup
       [feat,lbl] = getBootStrapSample(out.V,out.MSMaps,GroupedRandIndex,Group,out.Design,out.DLabels1,out.DLabels2,GroupLabels,out.bSmoothLabels,out.nWindowSize,out.LabelPenalty,out.StartFrame,out.EndFrame);
        BootFeatAll(:,b,F1F2GIndexS:F1F2GIndexE,:) = feat;
        BootLabelAll(F1F2GIndexS:F1F2GIndexE) = lbl;
    end
    waitbar(b/out.Iterations,h);
    set(h,'Name',sprintf('Remaining time: %01.0f:%02.0f min',floor(toc()*(out.Iterations/b-1)/60),rem(toc()*(out.Iterations/b-1),60)));
end
close(h)

function [feat,lbl] = getBootStrapSample(V,Maps,idx,Group,Condition,L1,L2,GL,bSmooth,ws,lp,sf,ef)

    [gm,lbl] = RaguGrandMeans(V(idx,:,:,:),Group,Condition,L1,L2,GL); % TK 14.3.2013
    gfp = std(gm,1,3);
    [MSClass,MSFit] = RaguFitMicrostates(gm,Maps,bSmooth,ws,lp,sf,ef);

    feat = nan(6,size(MSClass,1) * size(MSClass,2),size(Maps,1));
    
    rs = [size(MSClass,1) * size(MSClass,2),size(Maps,1)];
    [on,off,dur,auc,cog,msgfp] = RaguMSOnOffsets(MSClass,size(Maps,1),MSFit,gfp, 0);
    feat(1,:,:) = reshape(on,rs);
    feat(2,:,:) = reshape(off,rs);
    feat(3,:,:) = reshape(dur,rs);
    feat(4,:,:) = reshape(auc,rs);
    feat(5,:,:) = reshape(cog,rs);
    feat(6,:,:) = reshape(msgfp,rs);
   
     

    

