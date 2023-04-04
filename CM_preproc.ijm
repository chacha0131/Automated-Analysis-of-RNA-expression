//// Preprocessing WF nuclei images for RNAscope_quant ////
// Expects nuclei image selection
// Then loads and opens all channels of the same sample
// Allows selection of ROI; When pressing enter, crops and saves cropped images

// Identify Images
in = getDirectory("Parent directory"); 
imgs = in + "/raw/raw/"
out = in + "/preproc/"
files = getFileList(imgs);
setBatchMode(true);

for (i=0; i<17; i++) {

	// Open file.
	run("Bio-Formats Importer", "open=[" + imgs + files[i] + "] color_mode=Default " + 
		"rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	
	// Label channels and subset to only channels we need (i.e. remove middle channel).
	if (i==0) {
		Dialog.create("Channel selection");
		Dialog.addMessage("Provide channel names and identify which channels to keep");
		Stack.getDimensions(width,height,channels,slices,frames);
		Dialog.addNumber("nChar for ID:", 13);
		
		// Dialog box.
		for (j=0; j<channels; j++) {
			Dialog.addCheckbox("Keep Channel " + j + "?", false);
			Dialog.addString("Channel " + j + " Probe", "");
		}
		Dialog.addNumber("Enter number of DAPI channel.", 4);
		Dialog.show();

		// Dialog box output. 
		nChar = Dialog.getNumber();
		for (j=0; j<channels; j++) {
			if (j==0) {
				chan = Dialog.getString();
				keep = Dialog.getCheckbox();
			}
			if (j>0) {
				chan = Array.concat(chan, Dialog.getString());
				keep = Array.concat(keep, Dialog.getCheckbox());
			}
		}
		nucl = Dialog.getNumber();

		// Create folder structure in output folder.
		File.makeDirectory(out + "nucl_raw");
		if (i == 0) {
			for (j=0; j<channels; j++) {
				if (keep[j] == 1) {
					if (j != nucl) {
						File.makeDirectory(out + "rna_" + chan[j]);
					}
				}
			}
		}
	}

	// Get sample information. 
	id = getTitle();
	id = substring(id, 0, 0+nChar);
	
	// Split channels and save.
	for (j=0; j<channels; j++) {
		if (keep[j] == 1) {
			if (j == nucl) {
				title = id + "_nuclei";
				run("Duplicate...", "title=" + title + " duplicate channels=" + j+1);
				saveAs("Tiff", out + "nucl_raw/" + title + ".tif");
				close();
			} else {
				title = id + "_" + chan[j];
				run("Duplicate...", "title=" + title + " duplicate channels=" + j+1);
				saveAs("Tiff", out + "rna_" + chan[j] + "/" + title + ".tif");
				close();
			}
		}
	}
	close();
}
