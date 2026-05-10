function LBP= uniformLBP(inImg, ThLBP, varargin) % filtR, isRotInv, isChanWiseRot
%% efficientLBP
% The function implements LBP (Local Binary Pattern analysis).
%
%% Syntax
%  LBP= efficientLBP(inImg);
%
%% Description
% The LBP tests the relation between pixel and it's neighbors, encoding this relation into
%   a binary word. This allows detection of patterns/features.
% The function is inpired by materials published by Matti Pietik?inen in
%   http://www.cse.oulu.fi/CMV/Research/LBP . This implementation hovewer is not totally
%   allighned with the mthods proposed by Professor Pietik?inen (see Issues & Comments).
%
%% Input arguments (defaults exist):
% inImg- input image, a 2D matrix (3D color images will be converted to 2D intensity
%     value images)
% filtR- a 2D matrix representing a round/radial filter. It can be generated using
%   generateRadialFilterLBP function.
% isRotInv- a logical flag. When enabled generated rotation invariant LBP accuired via
%     fining an angle at whihc the LBP og a given pixelis minimal. Icreases run time, and
%     results in a relatively sparce hsitogram (as many combinations disappear).
% isChanWiseRot- a logical flag, when enabled (default value) allowes channel wise
%     rotation. When disabled/false rotation carried out based on roation of first color
%     channel. Supported only when "isEfficent" is enabled. When  "isEfficent" is
%     disabled "isChanWiseRot" is true.
%
%% Output arguments
%   LBP-    LBP image UINT8/UINT16/UINT32/UINT64/DOUBLE of same dimentions
%     [Height x Width] as inImg.
%
%% Issues & Comments
% - Currenlty, all neigbours are treated alike. Basically, we can use wighted/shaped
%     filter.
% - The rotation invariant LBP histogram includes less then bins then regular LBP BY
%     DEFINITION the zero trailing binary words are excluded for example, so it can be
%     reduced to a mush more component representation. Actually for 8 niegbours it's 37
%     bins, instead of 256. An efficnet way to calculate those bins value is needed.
%
%% Example
% img=imread('peppers.png');
% filtR=generateRadialFilterLBP(8, 1);
% tic;
%  % note this filter dimentions aren't legete...
% effLBP= efficientLBP(img, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false);
% effTime=toc;
%
% % verify pixel wise implementation returns same results
% tic;
% % same parameters as before
% pwLBP=pixelwiseLBP(img, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false); 
% inEffTime=toc;
% fprintf('\nRun time ratio %.2f. Same result eqaulity chesk: %o.\n', inEffTime/effTime,...
%    isequal(effLBP, pwLBP));
%
% figure;
% subplot(1, 3, 1)
% imshow(img);
% title('Original image');
%
% subplot(1, 3, 2)
% imshow( effLBP );
% title('Efficeint LBP image');
%
% subplot(1, 3, 3)
% imshow( pwLBP );
% title('Pixel-wise LBP image');
%
%% See also
% pixelwiseLBP  % a straigh forward iplmenetation of LBP, should achive same results
% generateRadialFilterLBP   % custom function generating circulat filters
%
%% Revision history
% First version: Nikolay S. 2012-05-01.
% Last update:   Nikolay S. 2014-01-09.
%
% *List of Changes:*
% 2014-01-16- support new radial filetr generation function
%   'generateRadialFilterLBP'. The new filter is 3D shapes, and it is alighned with 
%   "Gray Scale and Rotation Invariant Texture Classification with Local Binary Patterns" 
%   from http://www.ee.oulu.fi/mvg/files/pdf/pdf_6.pdf.
%   Changed filter direction (to CCW), starting point (3 o'clock instead of 12), support 
%   pixels interpolation.
% 2014-01-09- split pixel-wise implementation to a stand-alone function. Use
%   round/circular filter generated via generateRadialFilterLBP
% 2014-01-06 isChanWiseRot flag added to allow dictation of uniform rotation of all color
%   channels. Added witbar for the disabled 'isEfficent case, so the user will see it's
%   working...
% 2013-12-30 isRotInv flag added to allow minimal avalible LBP to result in LBP that is
%   rotation invaraint. When enabled, minimal possibel LBP will be calculated, via
%   rotating the neigborhood.
%   Inputs style cgange- support regulat values input, 'names' values pairs input, and
%   structure (where structure filed is the variable name, and it's conents is the value)
%   input.
% 2012-08-28 Neighbours were scanned column wise (regular Matlab way), while they should
%   be scanned clock-wise/counter clock-wise direction. A Helix/Snail indexing function
%   was written and added, to deal with this issue.
% 2012-08-27 Chris Forne comment mentioned some erros found in the code. As I haven't
%   made any use fo the code, for the last few month, I haven't noticed the mentioned
%   issues, so many thanks goes to Chris for his sharp eye. Bugs fixes, and some
%   modification intorduced.
% 2012-05-01 After writing down the primitove version, a filtering based miplementation
%   was proposed, improving run time by factor of 80-150..


%% Deafult params
filtR=generateRadialFilterLBP(8, 1);

%% Get user inputs overriding default values
funcParamsNames={'filtR'};
assignUserInputs(funcParamsNames, varargin{:});

if ischar(inImg) && exist(inImg, 'file')==2 % In case of file name input- read graphical file
    inImg=imread(inImg);
end

inImgType=class(inImg);
calcClass='uint8';

isCalcClassInput=strcmpi(inImgType, calcClass);
if ~isCalcClassInput
    inImg=cast(inImg, calcClass);
end
imgSize=size(inImg);

nNeigh=size(filtR, 3);

% Build LUT for uniform
% 0, 1, 2, 3, 4, 6, 7, 8, 12, 14, 15, 16, 24, 28, 30, 31, 32, 48, 56, 60, 62
% , 63, 64, 96, 112, 120, 124, 126, 127, 128, 129, 131, 135, 143, 159, 191, 
% 192, 193, 195, 199, 207, 223, 224, 225, 227, 231, 239, 240, 241, 243, 247,
% 248, 249, 251, 252, 253, 254, 255.
  %  0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100
LUTsqr = [0 1 1 2 1 9 2 3 1 9 9 9 2 9 3 4 %0000 0 15
       1 9 9 9 9 9 9 9 2 9 9 9 3 9 4 5 %0001 16 31
       1 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 %0010 32 47
       2 9 9 9 9 9 9 9 3 9 9 9 4 9 5 6 %0011 48 63
       1 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 %0100 64 79
       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 %0101 80 95
       2 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 %0110 96 111
       3 9 9 9 9 9 9 9 4 9 9 9 5 9 6 7 %0111 112 127
       1 2 9 3 9 9 9 4 9 9 9 9 9 9 9 5 %1000 128 143
       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 6 %1001 144 159
       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 %1010 160 175
       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 7 %1011 176 191
       2 3 9 4 9 9 9 5 9 9 9 9 9 9 9 6 %1100 192 207
       9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 7 %1101 208 223
       3 4 9 5 9 9 9 6 9 9 9 9 9 9 9 7 %1110 224 239
       4 5 9 6 9 9 9 7 5 6 9 7 6 7 7 8 %1111 240 255
       ];
LUT = uint8(LUTsqr(:));
weigthVec=reshape(2.^( (1:nNeigh) -1), 1, 1, nNeigh);
weigthMat=uint8(repmat( weigthVec, imgSize([1, 2]) ));
binaryWord=zeros(imgSize(1), imgSize(2), nNeigh, calcClass);

% Initiate neighbours relation filter and LBP's matrix
for iFiltElem=1:nNeigh
    % Rotate filter- to compare center to next neigbour
    filtNeight=filtR(:, :, iFiltElem);
    % calculate relevant LBP elements via filtering
    binaryWord(:, :, iFiltElem)=cast( ...
        abs(filter2( filtNeight, inImg(:, :), 'same' )) >= ThLBP,...
        calcClass );
    % Without rounding sometimes inaqulity happens in some pixels
    % compared to pixelwiseLBP
end % for iFiltElem=1:nNeigh
possibleLBP=uint8(sum(binaryWord.*weigthMat, 3));
LBP=intlut(possibleLBP, LUT);