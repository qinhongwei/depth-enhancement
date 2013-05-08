function Result=YangIteration(Image,InitialDepth,Height,Width,sigma_w,sigma_c,w,DepthInteval,IterativeTime)
     %% Initialization
    L=10000;
    k=1;
    D(:,:,1)=double(InitialDepth);
    CandidateD=0:DepthInteval:255;
    CostVolume=zeros(Height,Width,length(CandidateD));
    CostCW=zeros(Height,Width,length(CandidateD));

    %% Iterative Module
    while 1
        tic
        for i=1:length(CandidateD)
            CostVolume(:,:,i)=min(L,(CandidateD(i)-D(:,:,k)).^2);                   %Cost Volume C(i)
            CostCW(:,:,i) = bifilter2(Image,CostVolume(:,:,i),w,sigma_w,sigma_c);   %A bilateral ?ltering is performed throughout each slice of the cost volume to produce the new cost volume 
            % Compare with the reference, the color space is different  
        end
        [BestCost,BestDepthLocation]=min(CostCW,[],3);                          %Selecting the depth hypothesis with the minimal cost

        % Sub-pixel estimation
        CostUpper=zeros(Height,Width);
        CostLower=zeros(Height,Width);
        for i=1:length(CandidateD) 
            CostUpper=CostUpper+CostCW(:,:,i).*((BestDepthLocation+1)==i);
            CostLower=CostLower+CostCW(:,:,i).*((BestDepthLocation-1)==i);
        end
        k=k+1;
        D(:,:,k)=CandidateD(BestDepthLocation)-DepthInteval*(CostUpper-CostLower)./(2*(CostUpper+CostLower-2*BestCost));  
        % end of sub-pixel estimation   

        if IterativeTime==k
            BFTime=toc;
            fprintf('YI:The running time of bilateral iteration is %.5f s\n',BFTime)
            break;
        end
    end
    Result=D(:,:,IterativeTime);

end