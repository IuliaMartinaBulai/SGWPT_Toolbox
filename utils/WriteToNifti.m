%% Write volume into nifti
%     Created by: Anjali Tarun
% 
%     Input:
%     xsignal -- signal you want to save 
%     indices_wb -- mask that dictates brain volume coordinates
%     filename -- filename of the output nifti file
%     filename_in -- filename of a reference brain volume

function [] = WriteToNifti(filename_in, xsignal,indices_wb, filename)
    fHeader = spm_vol(filename_in);
    V = zeros(fHeader.dim);
    V(indices_wb) = xsignal;
    hdr=cbiReadNiftiHeader(fHeader.fname);
    cbiWriteNifti(filename,V,hdr,'float32');
end