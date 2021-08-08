% % 
v = VideoReader('Experiment-12 040919.avi');
% % 
video = zeros(v.Height, v.Width, v.Duration * v.FrameRate);
% % 
i = 1;
while hasFrame(v)
% %     
     rgb = readFrame(v);
     
     video(:,:,i) = double(rgb(:,:,1))./double(rgb(:,:,2));
     
% %     
     i = i + 1;
% %     
end

video (video < 15) = 0;


Fs = v.FrameRate;
T = 1/Fs; 
L = v.Duration * v.FrameRate;  
t = (0:L-1)*T; 

freq = NaN (208, 212);
fs = NaN(208, 212, 52);
P1s = NaN(208, 212, 52);

for row = 1:208
    for col = 1:212
        S = permute(video(row, col, :), [3 2 1]);
        
        Y = fft(S);

        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        f = Fs*(0:(L/2))/L;
        P1(f<0.3) = NaN;

        [~, i] = max(P1);
        if ~isempty(i)
            freq(row, col) = f(i);
        end
        fs(row, col, 1:size(f,2)) = f;
        P1s(row, col, 1:size(P1,1)) = P1;

    end
    if (mod(row, 10) == 0)
        row
    end
end

figure
imshow(freq)

imwrite (freq, 'Experiment-12 2 Final.tiff');

