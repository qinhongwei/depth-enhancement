%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function Name: AnisotropicDiffusion
%Aim: Guided depth upsample based on anisotropic diffusion 
%Output: 
%   Result      -   The output depth map
%Input: 
%   color       -   Color image
%   depth 		-   Sparse depth map 
%   sigma_w     -   Coefficient of gaussian kernel for spatial
%Code Author:
%   Liu Junyi, Zhejiang University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = AnisotropicDiffusion(color,depth,sigma_w)
    if( size(color,3) ~= 3 ),
		error( 'color data must be of 3 channel' );
    end
    height = size(color,1);
    width = size(color,2);
    pixelNumber = height * width;
    
    tic;
    depth = double(depth);
    Z = sparse(reshape(depth,pixelNumber,1));
   
    color = double(color);
    S = ADMatrix(color,depth,sigma_w);
    
    fprintf('AD:The running time of getting A and b is %.5f s\n',toc);
    Result = S\Z;
    BackslashTime=toc;
    fprintf('AD:The running time of solving Ax=b by Backslash is %.5f s\n',BackslashTime)
    
    
    
    result = full(reshape(double(Result),height,width));
    fprintf('AD:Done!\n')