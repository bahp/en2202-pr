% Author : Bernard Hernandez and Fernando Iglesias
% Date : 11/09/2012

% ------------------------------------------------------------------------
%                    LAB2 - FEATURE EXTRACTION
% ------------------------------------------------------------------------
clear; clc; close all;

% -------------------------------------------------------
%                   EDIT BY USER
% -------------------------------------------------------
% Options if saving images.
save_im = 1;
lab_num = '2';
image_num = 1;
images_dir = 'reports/figures';

% Relative path to sounds directory.
sounds_dir = 'data/sounds';
% -------------------------------------------------------

% List the .wav files in sounds_path.
wav_files = dir(sprintf('%s/*.wav',sounds_dir));

% Memory pre-allocation
N = length(wav_files);
s = struct('name',cell(1),'y',[],'fs',0,'nbits',0,'t',0);
sounds = repmat(s,1,N);

% Load the sound files.
for i = 1:N
    % name.
    sounds(i).name = wav_files(i).name;
    
    % y(samples), fs(sample frequency), nbits.
    [sounds(i).y, sounds(i).fs, sounds(i).nbits] = ...
        wavread(sprintf('%s/%s',sounds_dir,wav_files(i).name));
    
    % create x-axis to plot (t=1/f).
    nsamples = length(sounds(i).y);
    sounds(i).t = linspace(0, nsamples/sounds(i).fs, nsamples); % 0:1/fs:nsamples/fs
end

% Listen the .wav files.
sound(sounds(1).y,sounds(1).fs)

%%%%%%%%%%%%%%%%%%%%%%
%     INCLUDE 1      %
%%%%%%%%%%%%%%%%%%%%%%
% Plots in time domain.
for i=1:N
    
    % If we don't want to show it skip.
    if (~strcmp(sounds(i).name,'female.wav')) && ...
       (~strcmp(sounds(i).name,'music.wav'))
       continue;        
    end
    
    h = figure;
    
    % Plots (normal and augmented).
    subplot('211');
    plot(sounds(i).t, sounds(i).y);
    title(sprintf('Sound wave %s over time\n', sounds(i).name));
    xlabel('Time [s]');
    ylabel('Signal amplitude');
    axis([sounds(i).t(1) sounds(i).t(end) ... 
        min(sounds(i).y) max(sounds(i).y)])
    
    subplot('212')
    plot(sounds(i).t, sounds(i).y);
    title(sprintf('Sound wave %s over time\n', sounds(i).name));
    xlabel('Time [s]');
    ylabel('Signal amplitude');
     
    % Annotations.
    if strcmp(sounds(i).name,'female.wav') 
        axis([0.43 0.58 -0.22 0.25])
        annotation('ellipse',[0.3 0.63 0.05 0.15]) 
        annotation('textarrow', [.19 .19], [.35 .27], 'String' , 'silence');
        annotation('textarrow', [.4 .4], [.35 .27], 'String' , 'unvoiced');
        annotation('textarrow', [.61 .61], [.35 .27], 'String' , 'voiced');
    end
    
    if strcmp(sounds(i).name,'music.wav')
        axis([0.45 0.86 -0.1 0.1])
        annotation('ellipse',[0.21 0.63 0.05 0.17])
    end

    image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);
           
end

% Spectrograms and cepstrograms.
addpath GetSpeechFeatures
win_length = 30/1000; % window length of 30ms.
ncep = 13;            % number of cepstral coefficients.

s = struct('name',cell(1),'mfccs',[],'spectgram',[],'f',[],'t',[]);
grams = repmat(s,1,N);
for i = 1:N   
   grams(i).name = wav_files(i).name;
   % compute spectrogram and cepstrogram(mfccs).
   [mfccs, grams(i).spectgram, grams(i).f, grams(i).t ] = ...
       GetSpeechFeatures(sounds(i).y,sounds(i).fs, win_length, ncep);
    
   % cepstrogram normalization.
   grams(i).mfccs = ...
       (mfccs - repmat(mean(mfccs,2), [1 length(grams(i).t)])) ./ ...
       repmat(std(mfccs,0,2), [1 length(grams(i).t)]);
end;

%%%%%%%%%%%%%%%%%%%%%%
%     INCLUDE 2      %
%%%%%%%%%%%%%%%%%%%%%%
% Plot spectrograms.
for i = 1:N
    
    % If we don't want to show it skip.
    if (strcmp(sounds(i).name,'male.wav'))
       continue;        
    end
    
    h = figure();
    imagesc(grams(i).t,grams(i).f,10*log10(grams(i).spectgram));
    colorbar; axis xy;
    title(sprintf('Spectrogram of the sound file %s\n', sounds(i).name));
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
        
    % Annotations.
    if strcmp(sounds(i).name,'female.wav') 
        annotation('textbox',[0.170 0.125 0.05 0.75],'String',('unvoiced')');
        annotation('textbox',[0.54 0.125 0.03 0.75],'String',('unvoiced')');
        annotation('textarrow', [.35 .25], [.7 .138], 'String' , 'voiced');
        annotation('textarrow', [.35 .35], [.7 .138], 'String' , 'voiced');
    end

    if strcmp(sounds(i).name,'music.wav')
        annotation('textarrow', [.35 .17], [.6 .19], 'String' , 'Harmonic');
        annotation('textarrow', [.35 .19], [.5 .16], 'String' , 'Harmonic');
        annotation('textarrow', [.35 .21], [.4 .138], 'String' , 'Harmonic');
    end
        
    image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);
end

%%%%%%%%%%%%%%%%%%%%%%
%     INCLUDE 3      %
%%%%%%%%%%%%%%%%%%%%%%
% Plot spectrogram vs cepstrogram.
for i = 1:N
   
    % If we dont want to show it skip.
    if (strcmp(sounds(i).name,'male.wav'))
       continue;        
    end

    h = figure();
      
    subplot('211')
    imagesc(grams(i).t,grams(i).f,10*log10(grams(i).spectgram));
    colorbar; axis xy;
    title(sprintf('Spectrogram of the sound file %s\n', sounds(i).name));
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');

    subplot('212')
    imagesc(grams(i).t,[],grams(i).mfccs);
    colorbar; axis xy;
    title(sprintf('Cepstrogram of the sound file %s\n', sounds(i).name));
    xlabel('Time [s]');
    ylabel('Cepstral coefficients');
        
    image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);
end

% Plot spect female vs spect male.
h = figure();

pos = 1;
for i = 1:N
    if strcmp(grams(i).name,'male.wav') || strcmp(grams(i).name,'female.wav')
        
        % Spectrogram.
        subplot(strcat('21',int2str(pos)))
        imagesc(grams(i).t,grams(i).f,10*log10(grams(i).spectgram));
        colorbar; axis xy;
        title(sprintf('Spectrogram of the sound file %s\n', sounds(i).name));
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');
        
        pos = pos + 1;  
    end
end
   
image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);

% Plot cepstrogram female vs male.
h = figure();

pos = 1;
for i = 1:N
    if strcmp(grams(i).name,'male.wav') || strcmp(grams(i).name,'female.wav')
        subplot(strcat('21',int2str(pos)))
        imagesc(grams(i).t,[],grams(i).mfccs);
        colorbar; axis xy;
        title(sprintf('Cepstrogram of the sound file %s\n', sounds(i).name));
        xlabel('Time [s]')
        ylabel('Cepstral coefficients')
        
        pos = pos + 1;
    end
end

image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);

%%%%%%%%%%%%%%%%%%%%%%
%     INCLUDE 4      %
%%%%%%%%%%%%%%%%%%%%%%
% Correlation matrices spectral an cepstral coefficients.
for i = 1:N
    spec_corr = corr(10*log10(grams(i).spectgram'));
    ceps_corr = corr(grams(i).mfccs');
    
    h = figure();
    
    subplot(211)
    imagesc(abs(spec_corr))
    colormap gray
    colorbar
    title(sprintf(['Absolute value of the correlation matrix for the '... 
      'spectral coefficient series of %s\n'], sounds(i).name));
    xlabel('Frequency [Hz]')
    ylabel('Frequency [Hz]')
    
    subplot(212)
    imagesc(abs(ceps_corr))
    colormap gray
    colorbar
    title(sprintf(['Absolute value of the correlation matrix for the ' ... 
      'cepstral coefficient series of %s\n'], sounds(i).name));
    xlabel('Cepstral coefficients')
    ylabel('Cepstral coefficients')
    
    image_num = saveReportImage(h,images_dir,lab_num,image_num,'jpg',save_im);
    
end
