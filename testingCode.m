n = 100000;
for i = 1:n
    tic; bin2dec([dec2bin(6,8) dec2bin(194,8)]); 
    T(i)=toc;
end


%%

m = 6000;
highbyte = dec2bin(m,16);
highbyte = highbyte(1:8);
bin2dec(highbyte)