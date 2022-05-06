% write .nii file for decoding accuracies after combining them from the
% parallelized searchlight procedure 

function combine_write_decoding_results(beta_dir,results, fname)
% fill resultsvol 4D and write 4D nifi
backgroundvalue = 0;
% get canonical hdr from first preprocesed functional file
template_file = dir(fullfile(beta_dir,'wbeta*.nii'));
template_file = fullfile(template_file(1).folder,template_file(1).name);
hdr= spm_vol(template_file); % choose canonical hdr from first classification image
hdr = rmfield(hdr,'pinfo');
%hdr = rmfield(hdr, 'dt');

resultsvol_hdr = hdr;
resultsvol_hdr.fname = fname; 
resultsvol_hdr.descrip = sprintf('Decoding accuracies pairwise');
resultsvol = backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume with background value (default: 0)
resultsvol(results.mask_index) = results.accuracy_pairwise.output;
spm_write_vol(resultsvol_hdr,resultsvol);
end 