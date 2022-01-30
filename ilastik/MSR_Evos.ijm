//Run MultiStackReg on all paired stacks in all well_locations directories
//Called from the top level of a plate directory that contains well_location subdirectories
//Each subdirectory must have TxRed_Unreg and GFP_Unreg subdirectories with tiff images
//This macro will create TxRed_Reg and GFP_Reg subdirectories
//and store image sequences that were aligned to the TxRed sequence

function registerStacks(list) {	
  for (i=0; i<list.length; i++) {
      if(matches(list[i], "[A-Z][1-9]_[1-9]/")) {
      	 //create subdirectories for the registered images
        red_dir = dir + "/" + list[i] + "TxRed_Reg"+File.separator;
        if (!File.exists(red_dir))
        	File.makeDirectory(red_dir);
        	if (!File.exists(red_dir))
      			exit("Unable to create directory" + red_dir);
       	GFP_dir = dir + list[i] + "GFP_Reg"+File.separator;
       	if (!File.exists(GFP_dir))
        	File.makeDirectory(GFP_dir);
  			if (!File.exists(GFP_dir))
      			exit("Unable to create directory" + GFP_dir);
      	run("Image Sequence...", "open=" + dir + "/" + list[i] + "/TxRed_Unreg starting=2 sort");
      	run("Image Sequence...", "open=" + dir + "/" + list[i] + "/GFP_Unreg starting=2 sort");
      	print("Registering stacks in well_location " + red_dir);
      	run("MultiStackReg", "stack_1=TxRed_Unreg action_1=Align file_1=[] stack_2=GFP_Unreg action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body]");
		selectWindow("TxRed_Unreg");
		run("Image Sequence... ", "format=TIFF use save=" +  dir + "/" + list[i] + "/TxRed_Reg");
		close();
		selectWindow("GFP_Unreg");
		run("Image Sequence... ", "format=TIFF use save=" +  dir + "/" + list[i] + "/GFP_Reg");
		print("Done registering stacks in " + list[i]);
		close();
      }
  }
}
  
  //dir = getDirectory("Current");
  dir = "/eppec/storage/groups/heiserlab/image_scratch/HE_E_L_049_01_1/";
  print("Registering images in well_location directories of " + dir);
  files = getFileList(dir);
  registerStacks(files);
