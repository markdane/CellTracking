//Run MultiStackReg on all paired stacks in all well_locations directories
//Called from the top level of a plate directory that contains well_location subdirectories
//Each subdirectory must have R_Unreg, G_Unreg and P_Unreg subdirectories with tiff images
//This macro will create R_Reg, G_Reg and P_Reg subdirectories
//and store image sequences that were aligned to the R sequence

function registerStacks(list) {	
  for (i=0; i<list.length; i++) {
      if(matches(list[i], "[A-Z][1-9]_[1-9]/")) {
      	 //create subdirectories for the registered images
        red_dir = dir + "/" + list[i] + "R_Reg"+File.separator;
        if (!File.exists(red_dir))
        	File.makeDirectory(red_dir);
        	if (!File.exists(red_dir))
      			exit("Unable to create directory" + red_dir);
       	GFP_dir = dir + list[i] + "G_Reg"+File.separator;
       	if (!File.exists(GFP_dir))
        	File.makeDirectory(GFP_dir);
  			if (!File.exists(GFP_dir))
      			exit("Unable to create directory" + GFP_dir);
      	Phase_dir = dir + list[i] + "P_Reg"+File.separator;
       	if (!File.exists(Phase_dir))
        	File.makeDirectory(Phase_dir);
  			if (!File.exists(Phase_dir))
      			exit("Unable to create directory" + Phase_dir);
      	run("Image Sequence...", "open=" + dir + "/" + list[i] + "/R_Unreg starting=1 sort");
      	run("Image Sequence...", "open=" + dir + "/" + list[i] + "/G_Unreg starting=1 sort");
      	//print("Registering stacks in well_location " + red_dir);
      	print("copying to registered directory " + red_dir);
      	//run("MultiStackReg", "stack_1=R_Unreg action_1=Align file_1=[] stack_2=G_Unreg action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body]");
		selectWindow("R_Unreg");
		run("Image Sequence... ", "format=TIFF save=" +  dir + "/" + list[i] + "/R_Reg");
		close();
		print("copying to registered directory " + GFP_dir);
		selectWindow("G_Unreg");
		run("Image Sequence... ", "format=TIFF save=" +  dir + "/" + list[i] + "/G_Reg");
		//print("Done registering stacks in " + list[i]);
		print("Done copying stacks in " + dir  + list[i] + "G_Reg");
		close();
      }
  }
}
  
  //dir = getDirectory("Current"); //doesn't work in eppec compute node
  dir = "/eppec/storage/groups/heiserlab/image_scratch/LI_I_L_034_01_1/"; //uncomment for eppec
  //dir = "/Users/dane/Documents/CellTracking/IncuCyte/LI_I_L_034_01_1/"; //uncomment for eppec
  print("Registering images in well_location directories of " + dir);
  files = getFileList(dir);
  registerStacks(files);