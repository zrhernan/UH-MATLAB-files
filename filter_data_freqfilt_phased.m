function [ chfi ,h ,tranfcoefs ] = filter_data_freqfilt_phased( chi, fsi, n_fi, Wni, len )
    fnq= 1/(2*(1/fsi));  %nyquist frequency 
    Wn = [Wni(1)/fnq Wni(2)/fnq];
    [b,a] = butter(n_fi,Wn);
    tranfcoefs=[b;a];
    [h,t] = impz(b,a,len,fsi);

    chfi=cell(size(chi));
    for i=1:size(chi,2)
        %chfi{1,i}=filtfilt(b,a,chi{1,i});
        chfi{1,i}=fftfilt(h,chi{1,i});
    end
    
end

