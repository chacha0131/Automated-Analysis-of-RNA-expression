
// Identify Images
in = getDirectory("Parent directory"); 
out = in + "/preproc/rna_Sdk1_denoised/"
imgs = in + "/preproc/rna_Sdk1/"
files = getFileList(imgs);
nChar=7;
setBatchMode(true);

for (i=0; i<files.length; i++) {
	
	// Open file.
	run("Bio-Formats Importer", "open=[" + imgs + files[i] + "] color_mode=Default " + 
		"rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		
	// Get sample information. 
	id = getTitle();
	id = substring(id, 0, 0+nChar);

	// Run denoising function.
	run("Despeckle");

	// Save image. 
	title = id + "_Sdk1_denoised";
	saveAs("Tiff", out + title + ".tif");
	close();
	
}