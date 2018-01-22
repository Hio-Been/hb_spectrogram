function EEG = hb_spectrogram( EEG, option, chan_list, win_size, t_resolution )
%% Calculating spectrogram from EEG dataset
%
% paramIn) EEG=hb_spectrogram( EEG, calc_option, chan_list, win_size, t_resolution )
% example) EEG=hb_spectrogram( EEG,    'ePSD',    [1 2 3],     2^9,      0.1 )
%
% -- Possible calc_option --
% ePSD : event-related Power (amplitude) Spectral Density
% ePLF : event-related Phase Locking Factor
% eCOH : event-related spectral Coherence
% cPSD : continuous PSD (nTrials = 1)
% cPLF : continuous PLF (nTrials = 1)
% cCOH : continuous COH (nTrials = 1)
%
% if size( EEG.data, 3) == 1; then it returns continuous spectrogram.
%
% >>>> hiobeen.han@kaist.ac.kr. 2018-01-22.
% 

%% Defulat condition
if nargin < 5
    t_resolution = .1;
end; if nargin < 4
    win_size = 2^9;
end; if nargin < 3
    chan_list = 1:size(EEG.data,1);
end; if nargin < 2
    option = 'ePSD';
end;

% Make time vector 't' with second-resolution
if ~isfield(EEG, 'times')
    disp(['Time vector (EEG.times) does not exist!'])
    EEG.times = linspace(0, size(EEG.data,2)/EEG.srate, size(EEG.data,2) );
end
if mean(diff(EEG.times)) > 0.001
    disp(['Time vector (EEG.times) has milisecond resolution!'])
    t = EEG.times / 1000;
else
    t = EEG.times;
end
nfft = win_size;
t_fft = [t(1)+(((nfft*.5)+1)/EEG.srate), t(end)-(((nfft*.5)+1)/EEG.srate)];
t_vec = linspace( t_fft(1), t_fft(end), (diff(t_fft)/t_resolution) +1);

nTrials = size( EEG.data, 3 );
disp([ '-- Initiating hb_spectrogram.m --' ])
disp([ '. Calculating Spectrogram ... ' ])
disp([ '.. nTrials = ' num2str(nTrials ) ])
disp([ '... winSize = ' num2str( nfft/EEG.srate ) ' sec']);
disp([ '.... tResolution = ' num2str( t_resolution ) ' sec'])

nTrials = 30;


switch(option)
    
    case('ePSD') %% Event-related PSD
        
        EEG.PSD = [];
        EEG.PSD.data = []; EEG.PSD.nfft = nfft;
        EEG.PSD.freq_cut = 100; % HIGHST FREQUENCY
        disp(['... FreqBand = 0 to ' num2str(EEG.PSD.freq_cut) ' Hz']);

        for chanIdx = chan_list
            disp(['Calc ePSD.. Channel:' num2str(chanIdx) ] )
            for trialIdx = 1:nTrials
                epoch = EEG.data(chanIdx, :, trialIdx);
                psd = [];
                psd_t = [];
                for tIdx = 1:length(t_vec)
                    % (1) Indexing
                    idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                        max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
%                     idx(1)
                    actual_t = t(idx);
                    [x,f]= positiveFFT( hanning(length(epoch(idx)))' .* epoch(idx), EEG.srate, 0);
                    psd( size(psd,1)+1 , :) =  x;
                    psd_t( length(psd_t)+1 ) = mean(actual_t);
                end
                EEG.PSD.data( :, :, chanIdx, trialIdx ) = psd(:,1:max(find(f < EEG.PSD.freq_cut))+1);
            end
            EEG.PSD.t = psd_t; EEG.PSD.f = f(1:max(find(f < EEG.PSD.freq_cut))+1);
        end
        
    case('ePLF') %% Event-related PLF
        
        EEG.PLF = [];
        % PLF OPTION
        EEG.PLF.FreqVector = 2:1:100;
        EEG.PLF.BandWidth = [.999];
        EEG.PLF.filterOrder = 10;
        disp([...
            'n=' num2str(length( EEG.PLF.FreqVector) )...
            ',freq=' num2str(EEG.PLF.FreqVector(1)) '-' num2str(EEG.PLF.FreqVector(end)) 'Hz' ...
            ',width=' num2str( EEG.PLF.BandWidth ) 'Hz'...
            ',order=' num2str(EEG.PLF.filterOrder) ]);
        
        EEG.PLF.chanComb = combnk([chan_list],2);
        EEG.PLF.data = [];
        % Design Filters
        EEG.PLF.Filters = {};
        for jj=1:length(EEG.PLF.FreqVector)
            [EEG.PLF.Filters{jj,1}, EEG.PLF.Filters{jj,2}]= hb_getBandpassHd( ...
                ([-1 1]*EEG.PLF.BandWidth) +EEG.PLF.FreqVector(jj), EEG.PLF.filterOrder, EEG.srate );
        end
        % Calculation
        for combIdx = 1:size(EEG.PLF.chanComb,1)
            chan_a = EEG.PLF.chanComb(combIdx,1);
            chan_b = EEG.PLF.chanComb(combIdx,2);
            
            disp(['Calc ePLF.. Channel Combination:' num2str(combIdx) ] )
            for trialIdx = 1:nTrials
                
                epoch = EEG.data([chan_a, chan_b], :, trialIdx);
                plvs = [];
                for freqIdx = 1:length(EEG.PLF.FreqVector)
                    plv_t = [];
                    epoch_1 = hb_filtwithHds( epoch(1,:), EEG.PLF.Filters{freqIdx,1},  EEG.PLF.Filters{freqIdx,2});
                    epoch_2 = hb_filtwithHds( epoch(2,:), EEG.PLF.Filters{freqIdx,1},  EEG.PLF.Filters{freqIdx,2});
                    for tIdx = 1:length(t_vec)
                        % (1) Indexing
                        idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                            max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
                        
                        actual_t = t(idx);
                        
                        plv = hb_getPLV( epoch_1(idx), epoch_2(idx) );
                        plvs(length(plv_t)+1 , freqIdx) =  plv ;
                        plv_t(length(plv_t)+1) = mean(actual_t);
                    end
                end
                EEG.PLF.data( :, :, combIdx, trialIdx) = plvs;
            end
        end
        EEG.PLF.t = plv_t; EEG.PLF.f = EEG.PLF.FreqVector;
        
    case('eCOH') %% Event-related COH
        
        EEG.COH = [];
        EEG.COH.nfft = nfft;
        EEG.COH.chanComb = combnk([chan_list],2);
        EEG.COH.data = [];
        EEG.COH.freqCut=100;
        for combIdx = 1:size(EEG.COH.chanComb,1)
            chan_a = EEG.COH.chanComb(combIdx,1);
            chan_b = EEG.COH.chanComb(combIdx,2);
            disp(['Calc eCOH.. Channel Combination:' num2str(combIdx) ] )

            for trialIdx = 1:nTrials
                epoch = EEG.data([chan_a, chan_b], :, trialIdx);
                coh_t = [];
                epoch_1 = epoch(1,:);
                epoch_2 = epoch(2,:);
                cohs=[];
                for tIdx = 1:length(t_vec)
                    % (1) Indexing
                    idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                        max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
                    actual_t = t(idx);
                    [Cxy,F]=mscohere( epoch_1(idx), epoch_2(idx), ...
                        hanning(EEG.COH.nfft / 4), EEG.COH.nfft/8, EEG.COH.nfft, EEG.srate );
                    F2 = hb_findIdx([0 EEG.COH.freqCut], F);
                    cohs(length(coh_t)+1, 1:length(F2)) =  Cxy(F2) ;
                    coh_t( length(coh_t)+1 ) = mean(actual_t);
                end
                EEG.COH.data( :, :, combIdx, trialIdx) = cohs;
            end
        end
        EEG.COH.t = coh_t; EEG.COH.f = F(F2);

    case('cPSD') %% Continuous PSD

        EEG.PSD = [];
        EEG.PSD.data = []; EEG.PSD.nfft = nfft;
        EEG.PSD.freq_cut = 100;
        for chanIdx = chan_list
            epoch = EEG.data(chanIdx, :, 1);
%             if size(EEG.data,3)>1; error('Error:: nTrial is bigger than 1! use "cPSD" option!!'); end
            psd = [];
            psd_t = [];
            disp(['Calc cPSD.. Channel:' num2str(chanIdx) ] )

            for tIdx = 1:length(t_vec)
                % (1) Indexing
                idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                    max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
                
                actual_t = t(idx);
                [x,f]= positiveFFT( hanning(length(epoch(idx)))' .* epoch(idx), EEG.srate, 0);
                
                psd( size(psd,1)+1 , :) =  x;
                psd_t( length(psd_t)+1 ) = mean(actual_t);
            end
            
            EEG.PSD.data( :, :, chanIdx ) = psd(:,1:max(find(f < EEG.PSD.freq_cut))+1);
            EEG.PSD.t = psd_t; EEG.PSD.f = f(1:max(find(f < EEG.PSD.freq_cut))+1);
        end
        
    case('cPLF') %% Continuous PLF 
        
        EEG.PLF = [];
        EEG.PLF.chanComb = combnk(chan_list,2);
        EEG.PLF.FreqVector = 2:1:100;
        EEG.PLF.BandWidth = [.999];
        EEG.PLF.filterOrder = 10;
        EEG.PLF.data = [];
        EEG.PLF.Filters = {};
        for jj=1:length(EEG.PLF.FreqVector)
            [EEG.PLF.Filters{jj,1}, EEG.PLF.Filters{jj,2}]= hb_getBandpassHd( ...
                ([-1 1]*EEG.PLF.BandWidth) +EEG.PLF.FreqVector(jj), EEG.PLF.filterOrder, EEG.srate );
        end
        
        for combIdx = 1:size(EEG.PLF.chanComb,1)
            chan_a = EEG.PLF.chanComb(combIdx,1);
            chan_b = EEG.PLF.chanComb(combIdx,2);
            disp(['Calc cPLF.. Channel Combination:' num2str(combIdx) ] )
            epoch = EEG.data([chan_a, chan_b], :, 1);
%             if size(EEG.data,3)>1; error('Error:: nTrial is bigger than 1! use "cPLF" option!!'); end
            plvs = [];
            for freqIdx = 1:length(EEG.PLF.FreqVector)
                plv_t = [];
                epoch_1 = hb_filtwithHds( epoch(1,:), EEG.PLF.Filters{freqIdx,1},  EEG.PLF.Filters{freqIdx,2});
                epoch_2 = hb_filtwithHds( epoch(2,:), EEG.PLF.Filters{freqIdx,1},  EEG.PLF.Filters{freqIdx,2});
                for tIdx = 1:length(t_vec)
                    % (1) Indexing
                    idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                        max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
                    
                    actual_t = t(idx);
                    plv = hb_getPLV( epoch_1(idx), epoch_2(idx) );
                    plvs(length(plv_t)+1 , freqIdx) =  plv ;
                    plv_t( length(plv_t)+1 ) = mean(actual_t);
                end
            end
            EEG.PLF.data( :, :, combIdx) = plvs;
        end
        EEG.PLF.t = plv_t; EEG.PLF.f = EEG.PLF.FreqVector;
        
    case('cCOH') %% Continuous COH
        EEG.COH = [];
        EEG.COH.nfft = nfft;
        EEG.COH.chanComb = combnk([chan_list],2);
        EEG.COH.data = [];
        EEG.COH.freqCut=100;
        for combIdx = 1:size(EEG.COH.chanComb,1)
            chan_a = EEG.COH.chanComb(combIdx,1);
            chan_b = EEG.COH.chanComb(combIdx,2);
            disp(['Calc cCOH.. Channel Combination:' num2str(combIdx) ] )
            
            epoch = EEG.data([chan_a, chan_b], :, 1);
%             if size(EEG.data,3)>1; error('Error:: nTrial is bigger than 1! use "cCOH" option!!'); end
            coh_t = [];
            epoch_1 = epoch(1,:);
            epoch_2 = epoch(2,:);
            cohs=[];
            for tIdx = 1:length(t_vec)
                % (1) Indexing
                idx = max(find( t<(t_vec(tIdx)))) - nfft*.5 :...
                    max(find( t<(t_vec(tIdx))))   + nfft*.5 -1 ;
                
                actual_t = t(idx);
                [Cxy,F]=mscohere( epoch_1(idx), epoch_2(idx), ...
                    hanning(EEG.COH.nfft / 4), EEG.COH.nfft/8, EEG.COH.nfft, EEG.srate );
                F2 = hb_findIdx([0 EEG.COH.freqCut], F);
                cohs(length(coh_t)+1, 1:length(F2)) =  Cxy(F2) ;
                coh_t( length(coh_t)+1 ) = mean(actual_t);
            end
            EEG.COH.data( :, :, combIdx) = cohs;
            %                 end
        end
        EEG.COH.t = coh_t; EEG.COH.f = F(F2);
end

%% Subfunctions
function [X,freq]=positiveFFT(x,Fs,plotOption)
%  [X,freq]=positiveFFT(x,Fs,plotOption)

N=(length(x));  %get the number of points
k=0:N-1;        %create a vector from 0 to N-1
T=N/Fs;         %get the frequency interval
freq=k/T;       %create the frequency range
cutOff = ceil(N/2);
freq = freq(1:cutOff);


nTrials = size(x, 1);
if nTrials == 1
    X=fft(x)/N*2; % normalize the data
    %     disp(size(X))
    %only want the first half of the FFT, since it is redundant
    X = X(1:cutOff); % Single trial FFT
else
    X = [];
    for trialIdx = 1:nTrials
        X_temp = fft(x(trialIdx,:))/N*2;
        size(X_temp)
        X_temp = X_temp(1:cutOff); % Single trial FFT
        X = [X; X_temp];
    end
end

if plotOption
    for trialIdx = 1:nTrials
        plot(freq, abs((X(trialIdx,:))));
        hold on;
    end;
    xlim([1,161])
    xlabel('Freq (Hz)')
    ylabel('Normalized Amplitude')
    
    hold off;
end

return


function [yy] = hb_filtwithHds(xx, Hd1,Hd2)
% function [yy] = hb_filtwithHds(xx, Hd1,Hd2)

% Class convert
input_class= class(xx);
try
    if ~(min(input_class == 'double'))
        xx = double(xx);
    end
end

% Hd1 Filtering
filt_1 = filtfilt(Hd1.SOSMatrix, Hd1.ScaleValues, xx);
% Hd2 Filtering
filt_2 = filtfilt(Hd2.SOSMatrix, Hd2.ScaleValues, filt_1);
yy = filt_2;

return
function [highpassHd, lowpassHd] = hb_getBandpassHd(freqBand, filterOrder, Fs)
% function [highpassHd, lowpassHd] = hb_getBandpassHd(freqBand, filterOrder, Fs)


if nargin < 2
    filterOrder    = 8;     % Order
    Fs   = 2000;  % Sampling Frequency
end


h_high = fdesign.highpass('n,f3db', filterOrder, freqBand(1), Fs);
h_low = fdesign.lowpass('n,f3db', filterOrder, freqBand(2), Fs);

highpassHd = design(h_high, 'butter','SystemObject', true);
lowpassHd = design(h_low, 'butter','SystemObject', true);


return



function [PLV, angDiff, PLV_clean] = hb_getPLV(filtsig_x, filtsig_y, threshold)

if nargin < 3
    threshold = 3;
end

%  hb_getPLV function
% filtsig_x, _y, should be the shape of [ trialIdx x (time*srate) ],
% and each trial signal is hypothesized band-pass filtered already.

if min(size(filtsig_x) == size(filtsig_y)) == 0
    error('Signal X and Y should be the same size')
end

[nTrials, matDirection] = min(size(filtsig_x));
matDirection = 1;

if nTrials == 1
    xh = angle(hilbert(filtsig_x(:)));
    yh = angle(hilbert(filtsig_y(:)));
    %     xh = angle(hilbert(filtsig_x(:)) .* conj(hilbert(filtsig_x(:))));
    %     yh = angle(hilbert(filtsig_y(:)) .* conj(hilbert(filtsig_y(:))));
    angDiff = yh-xh;
    angDiff_clean = angDiff;
    angDiff_clean(find(angDiff< -1*std( angDiff )*threshold))=nan;
    angDiff_clean(find(angDiff>  1*std( angDiff )*threshold))=nan;
    PLV = abs(nanmean(exp(1i* angDiff)));
    PLV_clean = abs(nanmean(exp(1i* angDiff_clean)));
    % Batch trials
else
    PLV = zeros([nTrials,1]);
    angDiff = zeros([nTrials,size(filtsig_x,2)]);
    for trialIdx = 1:nTrials
        xh = angle(hilbert(filtsig_x(trialIdx,:)));
        yh = angle(hilbert(filtsig_y(trialIdx,:)));
        angDiff(trialIdx,:) = yh-xh;
        PLV(trialIdx,1) = abs(mean(exp(1i* angDiff(trialIdx,:))));
    end
end

return

function [idx,short] = hb_findIdx( range, fullData )
% function idx = hb_findIdx( range, fullData )

short=false;
idx = max(find( fullData < range(1)))+1:max(find(fullData<range(2)));

if length(idx)==0
%     idx
    idx = 1:max(find(fullData<range(2)));
end

if range(2) > fullData(end)
    short=1;
end

return

