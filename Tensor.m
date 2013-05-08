function T2=Tensor(Image,lamda,Height,Width)
    tic;
    %GrayImage = double(rgb2gray(uint8(Image)));
%     GrayImage=Image;
%      ProcessedImage = double(255*rgb2hsv(Image));
    ProcessedImage=double(Image);
    % Derivative masks 
%     dy = [-1 0 1; -1 0 1; -1 0 1];
%     dx = dy';
%     % Ix and Iy are the horizontal and vertical edges of image
%     Ix = conv2(GrayImage,dx,'same');
%     Iy = conv2(GrayImage,dy,'same');
    [Ix,Iy]=gradient(ProcessedImage);
    Ix=sqrt(Ix(:,:,1).^2+Ix(:,:,2).^2+Ix(:,:,3).^2);
    Iy=sqrt(Iy(:,:,1).^2+Iy(:,:,2).^2+Iy(:,:,3).^2);
    IxIyTime=toc;
    fprintf('LSLSTensor:The running time of computing Ix,Iy is %.5f s\n',IxIyTime)
    tic;
    % Generate Tensor
    V=cell(Height,Width);
    T1=cell(Height,Width);
    T2=cell(Height,Width);
    for i=1:Width     
        for j=1:Height     
            V(j,i)={[Ix(j,i);Iy(j,i);lamda*ProcessedImage(j,i,1);lamda*ProcessedImage(j,i,2);lamda*ProcessedImage(j,i,3)]};%
            T1(j,i)={V{j,i}*V{j,i}'+diag(exp(-20)*ones(5,1))}; % % Here the correction term has effects on the upsampling result
        end
    end
%     for i=1:Width     
%         for j=1:Height     
%             V(j,i)={[lamda*ProcessedImage(j,i,1);lamda*ProcessedImage(j,i,2);lamda*ProcessedImage(j,i,3)]};
%             T1(j,i)={V{j,i}*V{j,i}'+diag(exp(-20)*ones(3,1))}; % % Here the correction term has effects on the upsampling result
%         end
%     end
    % Gaussian Filter
    w=2;
    sigma_w=0.8;
    [X,Y] = meshgrid(-w:w,-w:w);
    G = exp(-(X.^2+Y.^2)/(2*sigma_w^2));  
    % G=fspecial('gaussian',[5,5],0.8);
    
    for i=1:Height    
        for j=1:Width
%            F=zeros(5,5);
            F=zeros(1,1);
            iMin = max(i-w,1);
            iMax = min(i+w,Height);
            jMin = max(j-w,1);
            jMax = min(j+w,Width);
            TTemp=T1(iMin:iMax,jMin:jMax);
            GTemp=G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
            w_i=iMax-iMin+1;
            w_j=jMax-jMin+1;
            for m=1:w_i
                for n=1:w_j
                    F=F+TTemp{m,n}.*GTemp(m,n);
                end
            end
            T2(i,j)={F/sum(GTemp(:))};
        end
    end
    TensorTime=toc;
    fprintf('LSLSTensor:The running time of generating all Tensors is %.5f s\n',TensorTime)
end