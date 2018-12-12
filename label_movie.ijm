//

if (is("Batch Mode")) {
	args = getArgument();
	if (args=="") exit ("no command line arguments");
	args = split(args, ",");
	movie_name = args[0];
	stack_path = args[1] + "/P_Unreg";
	ligand1 = args[2];
	drug1 = args[3];
} else {
	stack_path = "/Users/dane/Documents/CellTracking/IncuCyte/LI_I_L_035_01_1/A1_1/P_Unreg";
	movie_name = "OSM_10";
	ligand1 = "OSM 10 ng/ml";
	drug1 = "";
}

if (!File.exists(stack_path))
  	exit("Unable to access directory" + stack_path);
print("Opening stack in " +  stack_path);
run("Image Sequence...", "open=" + stack_path + " file=(tif) sort");
setForegroundColor(255, 255, 255);
setBackgroundColor(10, 10, 10);
run("Label...", "format=Text starting=0 interval=1 x=1200 y=1000 font=18 text=[" + ligand1 + "]");
run("Label...", "format=Text starting=0 interval=1 x=1200 y=1018 font=18 text=[" + drug1 + "]");
run("AVI... ", "compression=JPEG frame=10 save=" + stack_path + "/" + movie_name + ".avi");

