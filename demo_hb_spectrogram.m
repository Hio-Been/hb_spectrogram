clear EEG; 
load('demodata.mat');

%% (1) PSD
chan_list = [1 2];
win_size = 2^10; % t_window = win_size/srate;
t_resolution = .1; % sec

EEG = hb_spectrogram( EEG, 'ePSD', chan_list, win_size, t_resolution);

chan = 1;
subplot(2,1,chan);
imagesc( EEG.PSD.t, EEG.PSD.f, ...
    nanmean(abs(EEG.PSD.data(:,:,chan,:)), 4)');
ylim([0 60]);
axis xy;

chan = 2;
subplot(2,1,chan);
imagesc( EEG.PSD.t, EEG.PSD.f, ...
    nanmean(abs(EEG.PSD.data(:,:,chan,:)), 4)');
ylim([0 60]);
axis xy; colormap jet;



%% (3) COH
chan_list = [1 2];
win_size = 2^10; % t_window = win_size/srate;
t_resolution = .1; % sec

EEG = hb_spectrogram( EEG, 'eCOH', chan_list, win_size, t_resolution);

imagesc( EEG.COH.t, EEG.COH.f, ...
    nanmean(EEG.COH.data(:,:,1,:), 4)');
ylim([0 60]);
axis xy; colormap jet;




%% (2) PLF (very slow)
chan_list = [1 2];
win_size = 2^10; % t_window = win_size/srate;
t_resolution = .1; % sec

EEG = hb_spectrogram( EEG, 'ePLF', chan_list, win_size, t_resolution);

imagesc( EEG.PLF.t, EEG.PLF.f, ...
    nanmean(EEG.PLF.data(:,:,1,:), 4)');
ylim([0 60]);
axis xy; colormap jet;


