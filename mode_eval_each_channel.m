function modes=mode_eval_each_channel(zz,mfb)
ft=fft(zz);
for k=1:length(mfb)
    ewt{k}=real(ifft(conj(mfb{k})'.*ft));
    modes(k,:)=ewt{k};
end
end