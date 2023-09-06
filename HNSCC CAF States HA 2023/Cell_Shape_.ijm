inputfolder = getDirectory("Select input folder");
outputfolder = getDirectory("Select output folder");

// Ask the user to input multiple values
Dialog.create("Input Values");
Dialog.addNumber("Nucleus channel:", 1);
Dialog.addNumber("Min nucleus size:", 300);
Dialog.addNumber("Max nucleus size:", 121634816);
Dialog.addChoice("Select nucleus thresholding method:", newArray("Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"), "Otsu");
Dialog.addCheckbox("Exclude nuclei on edge", false);
Dialog.addNumber("Actin channel:", 2);
Dialog.addNumber("Min cell size:", 1000);
Dialog.addNumber("Max cell size:", 121634816);
Dialog.addChoice("Select actin thresholding method:", newArray("Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"), "IJ_IsoData");
Dialog.addNumber("Actin Gaussian blur radius:", 2.5);
Dialog.addCheckbox("Exclude cells on edge", false);
Dialog.show();

// Get the input values from the user
nuc_channel = Dialog.getNumber();
min_nuc = Dialog.getNumber();
max_nuc = Dialog.getNumber();
nuc_thresh = Dialog.getChoice();
nuc_ex = Dialog.getCheckbox();
actin_channel = Dialog.getNumber();
min_cell = Dialog.getNumber();
max_cell = Dialog.getNumber();
cell_thresh = Dialog.getChoice();
cell_blur = Dialog.getNumber();
cell_ex = Dialog.getCheckbox();

nuc_ex1 = "";
cell_ex1 = "";
if (nuc_ex) {
    nuc_ex1 = " exclude_objects_on_edges";
}
if (cell_ex) {
    cell_ex1 = " exclude_objects_on_edges";
}

// Find out how many files are in the folder
list = getFileList(inputfolder);
num2process = list.length;

// Loop through all files in the folder
for (l = 0; l < num2process; l++) 
{
	options = "open=[" + inputfolder + list[l] + "]";
    run("Bio-Formats Windowless Importer", options);
    
    // After opening, the image is in focus, so get its properties
    filename = getInfo("image.filename");
    current_title = getTitle(); 
    current_image = getImageID();
    selectImage(current_image);   // this puts the image in focus
    getVoxelSize(width, height, depth, unit);
    
    // With the image selected, process it
    run("Split Channels");
    
    // The channels will have the following names:
    channel1 = "C" + nuc_channel + "-" + current_title; 
    channel2 = "C" + actin_channel + "-" + current_title;
    
    // Create directory to store results
    new_output_folder = outputfolder + filename;
    File.makeDirectory(new_output_folder);

    selectImage(channel1);
    run("Make Binary", "method=" + nuc_thresh + " calculate create");
    selectImage("MASK_" + channel1);
    run("Fill Holes", "stack");
    run("Despeckle", "stack");
    run("Close-", "stack");
    run("Open", "stack");
    run("Despeckle", "stack");
    run("3D Objects Counter", "threshold=128 slice=1 min.=" + min_nuc + " max.=" + max_nuc + " objects statistics" + nuc_ex1);
    selectImage("Objects map of MASK_" + channel1);
    run("Roundness 2D 3D", "xsize=1.0000 ysize=1.0000 zsize=1.0000 lower=1 upper=255 entire calculation=Volumetric");
    // Save the results and name the file appropriately
    objects1 = "Objects map of MASK_" + channel1;
    stats1 = "Statistics for MASK_" + channel1;
    selectImage(objects1);
    saveAs("Tiff", new_output_folder + "/" + objects1 + ".tiff");
    close(objects1);
    close(channel1);
    close("MASK_" + channel1);
    selectWindow(stats1);
    saveAs("Results", new_output_folder + "/" + stats1 + ".csv");
    close(stats1 + ".csv");
    selectImage("roundness");
    saveAs("Tiff", new_output_folder + "/" + objects1 + " roundness.tiff");
    close(objects1 + " roundness.tiff");
    
    selectImage(channel2);
    Stack.setXUnit(unit);
    run("Properties...", "pixel_width=" + width + " pixel_height=" + height + " voxel_depth"= + depth);
    run("Duplicate...", "title=dup duplicate");
    selectImage("dup");
    Stack.setXUnit(unit);
    run("Properties...", "pixel_width=" + width + " pixel_height=" + height + " voxel_depth"= + depth);
    run("Gaussian Blur 3D...", "x=" + cell_blur + " y=" + cell_blur + " z=" + cell_blur);
    selectImage("dup");
    run("Make Binary", "method=" + cell_thresh + " calculate create");
    selectImage("MASK_dup");
    run("Fill Holes", "stack");
    run("Remove Outliers...", "radius=5 threshold=50 which=Dark stack");
    run("Open", "stack");
    run("Close-", "stack");
    run("Despeckle", "stack");
    run("3D Objects Counter", "threshold=128 slice=1 min.=" + min_cell + " max.=" + max_cell + " objects" + cell_ex1);
    
    selectImage(channel2);
    run("Despeckle", "stack");
    run("Seeded Region Growing ...", "image=[" + channel2 + "] seeds=[Objects map of MASK_dup] stack=[3D volume]");
    srg = channel2 + "1-SRG";
    selectImage(srg);
    run("Measure");
    modeIndex = getResult("Mode");
    close("Results");
    selectImage(srg);
    run("Duplicate...", "title=dup1 duplicate");
    selectImage("dup1");
    setThreshold(1, modeIndex-1, "raw");
    run("Convert to Mask", "background=Dark black create");
    selectImage(srg);
    run("Subtract...", "value=" + modeIndex + " stack");
    run("Multiply...", "value=255 stack");
    imageCalculator("Add create stack", srg, "MASK_dup1");
    selectImage("Result of " + srg);
    run("Convert to Mask", "background=Default create");
    run("Fill Holes", "stack");
    run("Open", "stack");
    run("Close-", "stack");
    selectImage("MASK_Result of " + srg);
    Stack.setXUnit(unit);
    run("Properties...", "pixel_width=" + width + " pixel_height=" + height + " voxel_depth"= + depth);
    run("3D Objects Counter", "threshold=128 slice=1 min.=" + min_cell + " max.=" + max_cell + " statistics objects" + cell_ex1);
    selectImage("Objects map of MASK_Result of " + srg);
    run("Roundness 2D 3D", "xsize=1.0000 ysize=1.0000 zsize=1.0000 lower=1 upper=255 entire calculation=Volumetric");
    
    // Save the results and name the file appropriately
    objects2 = "Objects map of MASK_Result of " + srg;
    stats2 = "Statistics for MASK_Result of " + srg;
    selectImage(objects2);
    saveAs("Tiff", new_output_folder + "/" + objects2 + ".tiff");
    selectWindow(stats2);
    saveAs("Results", new_output_folder + "/" + stats2 + ".csv");
    selectImage("roundness");
    saveAs("Tiff", new_output_folder + "/" + objects2 + " roundness.tiff");
    run("Close All");
    close(stats2 + ".csv");
}