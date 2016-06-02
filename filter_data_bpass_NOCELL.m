function [ chfi ] = filter_data_bpass_NOCELL( chi, fsi, n_fi, Wni )

    fnq= 1/(2*(1/fsi));  %nyquist frequency 
    Wn = [Wni(1)/fnq Wni(2)/fnq];
    [b,a] = butter(n_fi,Wn);

    chfi=zeros(size(chi));
    for i=1:size(chi,1)
        chfi(i,:)=filtfilt(b,a,chi(i,:));
    end

end

