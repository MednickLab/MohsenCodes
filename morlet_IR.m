function  [coef EP_phase EP_energy]= morlet_IR(sig, samp_rate, frecs, wave_num, timescale, is_plot)

%Usage
%   coef = spec_morlet(sig, samp_rate, frecs, wave_num, timescale, is_plot)
%Inupts
%   sig             1-d signal
%   samp_rate       sampling rate
%   frecs           vector with central frequencies
%   wave_num        wave number parameter
%   timescale       number of std for the initial span of the wavelet
%   (before cutting to the 1% criterion)
%   is_plot         1,plot;0,not plot
%Outputs
%   coefficients    time_frequency distribution.
%Algorithm
%   The signal was convoluted by complex Morlet's wavelets w(t,f0) having
%a Gaussian shape both in the time domain(SD_t) and in the frequency
%domain(SD_f) around its central frequency f0:w(t,f0) =
%(SD_t*sqrt(pi))^(-1/2) * exp( -t.^2/(2*SD_t^2) ) .* exp(i*2*pi*f0*t),
%with SD_f = 1/(2*pi*SD_t). f0/SD_f = wave_num. At 20 Hz,this leads to a wavelet 
% duration (2*SD_t) of 111.4 ms
%and to a spectral bandwidth (2*SD_f) of 5.8 Hz and at 100 Hz, to a duration
%of 22.2 ms and a bandwidth of 28.6 Hz.The time resolution of this method,
%therefore, increase with frequency,whereas the frequency resolution
%decreases.

% we can select dyadic scales, so in this case prior to run the fuction one
% needs to evaluate the following lines frec_power_step = 4; , max_freq_true=maximum frequency
% %     max_frec = ceil(frec_power_step*log2(max_freq_true));
% %     frecs = 2.^((0:1:max_frec)/frec_power_step);


sig = sig(:)';
len_sig = length(sig);
samp_period = 1/samp_rate;

%initialize output coefficients matrix
row_coef = length(frecs);
col_coef = len_sig;
coef = zeros(row_coef,col_coef);

%compute coefficients
for k = 1:row_coef
    SD_f = frecs(k)/wave_num;   %wave_num=4; , timescale=5
    SD_t = 1/(2*pi*SD_f);
    x = -timescale*SD_t:samp_period:timescale*SD_t;
    Morlets = ( SD_t*sqrt(pi) )^(-1/2) * exp( -x.^2/(2*SD_t^2) ) .* exp(1i*2*pi*frecs(k)*x );
%     Morlets = ( SD_t*sqrt(pi) )^(-1/2) * exp( -x.^2/(2*SD_t^2) ) .* (exp(1i*2*pi*frecs(k)*x ) - exp(((2*pi*frecs(k)*SD_t)^2)/2 ));
    Morlets = Morlets(find(abs(Morlets)>=max(abs(Morlets))/100,1):find(abs(Morlets)>=max(abs(Morlets))/100,1,'last'));
    %figure; plot(abs(Morlets))
    coef_freq = conv(sig,Morlets);
    coef(k,:) = coef_freq(round(length(Morlets)/2):col_coef+round(length(Morlets)/2)-1 );
    EP_phase(k,:)=angle(coef(k,:));
    EP_energy(k,:)=abs(coef(k,:)).^2;
end

%plot
if is_plot
    colormap('jet');
    imagesc([0:samp_period:(col_coef-1)*samp_period], frecs,abs(coef).^2);
    title('energy');
    xlabel('time (second)');
    ylabel('frequency (Hz)');
    axis('xy');
    colorbar('vert');
    shading interp;
    zoom on;
end