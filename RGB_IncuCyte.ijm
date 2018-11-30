//Run MultiStackReg on all paired stacks in all well_locations directories
//Each subdirectory in dir must have R_Unreg, G_Unreg and P_Unreg subdirectories with tiff images
//This macro will create R_Reg, G_Reg and P_Reg subdirectories
//and store image sequences that were aligned to the R sequence

function safeSliceDelete(bad_image) {
	   	if (nSlices >= bad_image) {
      		setSlice(bad_image);
      		run("Delete Slice");
      		setSlice(1);
      	}
}

function registerStacks(dir_list) {	
  for (i=0; i<dir_list.length; i++) {
      if(matches(dir_list[i], "[A-Z][1-9]_[1-9]/")) {
      	 //create subdirectories for the registered images
        red_dir = dir + "/" + dir_list[i] + "R_Reg"+File.separator;
        if (!File.exists(red_dir))
        	File.makeDirectory(red_dir);
        	if (!File.exists(red_dir))
      			exit("Unable to create directory" + red_dir);
       	GFP_dir = dir + dir_list[i] + "G_Reg"+File.separator;
       	if (!File.exists(GFP_dir))
        	File.makeDirectory(GFP_dir);
  			if (!File.exists(GFP_dir))
      			exit("Unable to create directory" + GFP_dir);
      	Phase_dir = dir + dir_list[i] + "P_Reg"+File.separator;
       	if (!File.exists(Phase_dir))
        	File.makeDirectory(Phase_dir);
  			if (!File.exists(Phase_dir))
      			exit("Unable to create directory" + Phase_dir);
      	print("Opening unregistered R stack in well_location " +  dir + dir_list[i]);
      	run("Image Sequence...", "open=" + dir + "/" + dir_list[i] + "/R_Unreg starting=2 sort");
      	//remove slice where the image is blank due to plate removal
      	blank_slice_1 = 66;
      	blank_slice_2 = 73;
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	//enhance R_Unreg images
      	run("Duplicate...", "title=copy duplicate");
		selectWindow("copy");
		run("Gaussian Blur...", "sigma=1 stack");
		//setAutoThreshold("Moments dark");
		//run("Convert to Mask", "method=Moments background=Dark calculate");
		setAutoThreshold("Moments dark no-reset");
		//run("Threshold...");
		call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
		setOption("BlackBackground", false);
		run("Convert to Mask", "method=Moments background=Dark calculate");
		run("Minimum...", "radius=1 stack");
		run("Gaussian Blur...", "sigma=2 stack");
		run("Convert to Mask", "method=Default background=Light calculate");
		//print("using mask:" +  dir_list[i]);
		imageCalculator("AND stack", "R_Unreg","copy");
		selectWindow("R_Unreg");
		run("Enhance Contrast", "saturated=0.35");
		rename("R_Unreg_contrast");
		print("Registering R stack in well_location " +  dir + dir_list[i]);
      	run("MultiStackReg", "stack_1=R_Unreg_contrast "
      	+ "action_1=Align "
      	+ "file_1=[Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore file_2=[] "
      	+ "transformation=Translation save");
		selectWindow("R_Unreg_contrast");
		rename("R_Reg");
		run("8-bit");
		//print("Saving registered R stack in well_location " +  dir + dir_list[i]);
		run("Image Sequence... ", "format=TIFF save=" +  dir + dir_list[i] + "/R_Reg");
      	print("Opening unregistered G stack in well_location " +  dir + dir_list[i]);
      	run("Image Sequence...", "open=" + dir + dir_list[i] + "/G_Unreg starting=2 sort");
      	//remove slice where the image is blank due to plate removal
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	print("Registering G stack from " +  dir + dir_list[i]);
      	run("MultiStackReg", "stack_1=G_Unreg "
      	+ "action_1=[Load Transformation File]"
      	+ "file_1=[Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore "
      	+ "file_2=[] "
      	+ "transformation=Translation");
		selectWindow("G_Unreg");
		rename("G_Reg");
		run("8-bit");
		print("Saving registered G stack in well_location " +  dir + dir_list[i]);
		run("Image Sequence... ", "format=TIFF save=" +  dir + "/" + dir_list[i] + "/G_Reg");
		print("Opening unregistered P stack in well_location " +  dir + dir_list[i]);
		run("Image Sequence...", "open=" + dir + "/" + dir_list[i] + "/P_Unreg starting=2 sort");
		//remove slice where the image is blank due to plate removal
      	safeSliceDelete(blank_slice_1);
		safeSliceDelete(blank_slice_2);
      	print("Registering P stack from " +  dir + dir_list[i]);
      	run("MultiStackReg", "stack_1=P_Unreg "
      	+ "action_1=[Load Transformation File]"
      	+ "file_1=[Transformation.txt] "
      	+ "stack_2=None "
      	+ "action_2=Ignore "
      	+ "file_2=[] "
      	+ "transformation=Translation");
		selectWindow("P_Unreg");
		rename("P_Reg");
      	//merge registered R contrast, G and P stacks
      	print("Merging channels of registered R contrast, G and P stacks in well_location " +  dir + "/" + dir_list[i]);
      	run("Merge Channels...", "c1=R_Reg c2=G_Reg c3=P_Reg create");
       	selectWindow("Composite");
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		print("Saving merged RGB for " +  dir + dir_list[i]);
		saveAs("Tiff",  dir + "/" + dir_list[i] + "/P_Reg/Composite_reg.tif");
  		while (nImages>0) {
  			selectImage(nImages);
  			close();
  			} 
		print("Done registering stacks for " +  dir + "/" + dir_list[i]);
      }
  }
}
  
  dir = "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_035_01_1/"; //uncomment for eppec
  //dir = "/Users/dane/Documents/CellTracking/IncuCyte/PL42/"; //uncomment for laptop
  print("Registering images in well_location directories of " + dir);
  dirs = getFileList(dir);
  registerStacks(dirs);
