% External HDR-VDP-3 runner. No changes to its model or calibration parameters.
args = argv();
cfg = jsondecode(fileread(args{1}));
if ~isempty(cfg.image_package_list), pkg('local_list',cfg.image_package_list); end
pkg load image;
if ~isempty(cfg.statistics_package_list), pkg('local_list',cfg.statistics_package_list); end
pkg load statistics;
addpath(cfg.hdrvdp);
addpath(fullfile(cfg.hdrvdp,'utils'));
addpath(fullfile(fileparts(mfilename('fullpath')),'octave_compat'));
results = {};
for k=1:length(cfg.cases)
  c=cfg.cases(k);
  d=load(c.file);
  if isfield(c,'identical') && ~isempty(c.identical) && c.identical
    d.comparison=d.reference;
  end
  tic;
  opts={'use_gpu',false,'mtf','cie','age',c.age,'surround',d.background};
  r=hdrvdp3('detection',d.comparison,d.reference,'luminance',d.ppd,opts);
  roi=d.roi>0;
  s=struct('name',c.name,'task','detection','model','HDR-VDP-3.0.7', ...
           'mtf','CIE99','age',c.age,'pixels_per_degree',d.ppd, ...
           'P_det',r.P_det,'C_max',r.C_max,'interior_C_max',max(r.C_map(roi)), ...
           'interior_P_max',max(r.P_map(roi)),'interior_P_median',median(r.P_map(roi)), ...
           'seconds',toc);
  results{end+1}=s;
  P_map=r.P_map; C_map=r.C_map;
  % Optical image is diagnostic radiance after ocular scatter, not a display image.
  pars=hdrvdp_parse_options('detection',opts);
  rho=create_cycdeg_image(size(d.reference)*2,d.ppd,false);
  mtf=hdrvdp_mtf(rho,pars);
  optical=fast_conv_fft(d.reference,mtf,d.background);
  save('-mat7-binary',fullfile(cfg.out,[c.name '-prediction.mat']),'P_map','C_map','optical');
  fid=fopen(fullfile(cfg.out,'predictions.json'),'w'); fprintf(fid,'%s\n',jsonencode(results)); fclose(fid);
  fprintf('%s\n',jsonencode(s)); fflush(stdout);
end
