function rd = CollapseAcrossXAxis(rd, Border, Label)


V = zeros(size(rd.V,1),size(rd.V,2),size(rd.V,3),size(Border,1));

for i = 1:size(Border,1) % TK 7.2.2017
    
    StartFrame = round((Border(i,1)-rd.TimeOnset) / rd.DeltaX) +1;
    EndFrame   = round((Border(i,2)-rd.TimeOnset) / rd.DeltaX) +1;

    V(:,:,:,i) = mean(rd.V(:,:,:,StartFrame:EndFrame),4);
end

rd.V = V;
rd.MeanInterval = false;
rd.TimeOnset = 1;
rd.DeltaX    = 1;
rd.StartFrame = 1;
rd.EndFrame = size(V,4);
rd.txtX = Label;
