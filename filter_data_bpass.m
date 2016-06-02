function [ chfi ] = filter_data_bpass( chi, fsi, n_fi, Wni )

    fnq= 1/(2*(1/fsi));  %nyquist frequency 
    Wn = [Wni(1)/fnq Wni(2)/fnq];
    [b,a] = butter(n_fi,Wn);

    chfi=cell(size(chi));
    for i=1:size(chi,2)
        chfi{1,i}=filtfilt(b,a,chi{1,i});
    end

end

