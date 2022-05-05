function MovieToGIF(inputMovie)

inputMovie = outputMovie;
numFrames = numel(inputMovie);
frameSize = size(inputMovie(1).cdata);
mov = zeros([frameSize, numFrames]);
for idx=1:numFrames
    mov(:,:,:,idx) = inputMovie(idx).cdata;
end

