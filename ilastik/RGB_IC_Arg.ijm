//Run MultiStackReg on all paired stacks in the command line argument well_locations directory
//The subdirectory in dir must have R_Unreg, G_Unreg and P_Unreg subdirectories with tiff images
//This macro will create a P_Reg subdirectory
//and store image sequences that were aligned there

function safeSliceDelete(bad_image) {
	   	if (nSlices >= bad_image) {
      		setSlice(bad_image);
      		run("Delete Slice");
      		setSlice(1);
      	}
}

function registerStacks(dir) {	
      	 //create subdirectory for the registered images
      	Output_dir = dir + "/Analysis"+File.separator;
       	if (!File.exists(Output_dir))
        	File.makeDirectory(Output_dir);
  			if (!File.exists(Output_dir))
      			exit("Unable to create directory" + Output_dir);
      	print("Opening unregistered R stack in " +  dir);
      	run("Image Sequence...", "open=" + dir + "/R_Unreg starting=2 sort");
      	//remove slice where the image is blank due to plate removal
      	blank_slice_1 = 66;
      	blank_slice_2 = 73;
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	//enhance R_Unreg images
      	run("Duplicate...", "title=copy duplicate");
		selectWindow("copy");
		run("Gaussian Blur...", "sigma=1 stack");
		setAutoThreshold("Moments dark no-reset");
		call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
		setOption("BlackBackground", false);
		run("Convert to Mask", "method=Moments background=Dark calculate");
		run("Minimum...", "radius=1 stack");
		run("Gaussian Blur...", "sigma=2 stack");
		run("Convert to Mask", "method=Default background=Light calculate");
		imageCalculator("AND stack", "R_Unreg","copy");
		selectWindow("R_Unreg");
		run("Enhance Contrast", "saturated=0.35");
		rename("R_Unreg_contrast");
		print("Registering R stack in " +  dir);
		run("MultiStackReg", "stack_1=R_Unreg_contrast "
		+ "action_1=Align "
      	+ "file_1=[" + Output_dir + "Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore file_2=[] "
      	+ "transformation=Translation save");
		selectWindow("R_Unreg_contrast");
		rename("R_Reg");
		run("8-bit");
		print("Opening unregistered G stack in " +  dir);
		run("Image Sequence...", "open=" + dir + "/G_Unreg starting=2 sort");
		//remove slice where the image is blank due to plate removal
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	print("Registering G stack from " +  dir);
      	run("MultiStackReg", "stack_1=G_Unreg "
      	+ "action_1=[Load Transformation File]"
      	+ "file_1=[" + Output_dir + "Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore "
      	+ "file_2=[] "
      	+ "transformation=Translation");
		selectWindow("G_Unreg");
		rename("G_Reg");
		run("8-bit");
		print("Opening unregistered P stack in well_location " +  dir);
		run("Image Sequence...", "open=" + dir + "/P_Unreg starting=2 sort");
		//remove slice where the image is blank due to plate removal
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	print("Registering P stack from " +  dir);
      	run("MultiStackReg", "stack_1=P_Unreg "
      	+ "action_1=[Load Transformation File]"
      	+ "file_1=[" + Output_dir + "Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore "
      	+ "file_2=[] "
      	+ "transformation=Translation");
		selectWindow("P_Unreg");
		rename("P_Reg");
    	//merge registered R contrast, G and P stacks
      	print("Merging channels of registered R contrast, G and P stacks in " +  dir);
      	run("Merge Channels...", "c1=R_Reg c2=G_Reg c3=P_Reg create");
       	selectWindow("Composite");
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		print("Saving merged RGB for " +  dir);
		saveAs("Tiff",  dir + "/Analysis/Composite_reg.tif");
  		while (nImages>0) {
  			selectImage(nImages);
  			close();
  			} 
		print("Done registering stacks for " +  dir);
}

  setBatchMode(true);
  dir = getArgument();
  print("Registering images in well_location directories of " + dir);
  registerStacks(dir);
