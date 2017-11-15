function cols=RDMcolormap
% this function provides a convenient colormap for visualizing
% dissimilarity matrices. it goes from blue to yellow and has grey for
% intermediate values.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

nCols = 256;
%% blue-cyan-gray-red-yellow with increasing V (BCGRYincV)
anchorCols=[0 0 1
    0 1 1
    .5 .5 .5
    1 0 0
    1 1 0];

anchorCols_hsv=rgb2hsv(anchorCols);
incVweight=1;
anchorCols_hsv(:,3)=(1-incVweight)*anchorCols_hsv(:,3)+incVweight*linspace(0.5,1,size(anchorCols,1))';

brightness(anchorCols);
anchorCols=hsv2rgb(anchorCols_hsv);

cols=colorScale(anchorCols,nCols);
cols1=cols;
    function b=brightness(RGBrows)
        % given triples of RGB values, returns the overall brightness vector. This
        % is performed on each row of the input (RGBrows) independently.
        %__________________________________________________________________________
        % Copyright (C) 2010 Medical Research Council
        
        RGBweights=[.241 .691 .068]';
        
        b=sqrt(RGBrows*RGBweights);
    end

    function cols=colorScale(anchorCols,nCols,monitor)
        % linearly interpolates between a set of given 'anchor' colours to give
        % nCols and displays them if monitor is set
        %__________________________________________________________________________
        % Copyright (C) 2012 Medical Research Council
        
        
        %% preparations
        if ~exist('monitor','var'), monitor=false; end
        
        
        %% define color scale
        nAnchors=size(anchorCols,1);
        cols = interp1((1:nAnchors)',anchorCols,linspace(1,nAnchors,nCols));
        
        
        %% visualise
        if monitor
            figure(123); clf;
            imagesc(reshape(cols,[nCols 1 3]));
        end
    end
end