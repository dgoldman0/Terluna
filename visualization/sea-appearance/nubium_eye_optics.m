% Display preparation only: official CIE99 optics applied to direct Earth XYZ.
% XYZ channels share an achromatic optical kernel; neural colour is not solved.
args=argv(); cfg=jsondecode(fileread(args{1}));
if ~isempty(cfg.image_package_list), pkg('local_list',cfg.image_package_list); end
pkg load image;
addpath(cfg.hdrvdp); addpath(fullfile(cfg.hdrvdp,'utils'));
addpath(cfg.compat);
d=load(cfg.input);
rho=create_cycdeg_image(size(d.direct(:,:,1))*2,d.ppd,false);
mtf=hdrvdp_otf_cie99(rho,cfg.age);
assert(all(isfinite(mtf(:))) && max(mtf(:))<=1+1e-12 && min(mtf(:))>=0);
probe=hdrvdp_otf_cie99([0 logspace(-3,2,1000)],cfg.age);
assert(probe(1)==1 && all(diff(probe)<=1e-12));
optical=zeros(size(d.direct));
for ch=1:3
  optical(:,:,ch)=fast_conv_fft(d.direct(:,:,ch),mtf,0);
end
negative_fraction=-sum(min(optical(:),0))/sum(optical(:));
assert(negative_fraction<1e-8);
% Bound and report tiny FFT undershoot before zeroing it.
optical=max(optical,0);
fprintf('FFT negative-energy fraction zeroed: %.12g\n',negative_fraction);
foreground=d.foreground; ppd=d.ppd;
save('-mat7-binary',cfg.output,'optical','foreground','ppd','negative_fraction');
fprintf('CIE99 age %g: finite, nonnegative, monotone MTF; DC = 1; optical XYZ saved\n',cfg.age);
