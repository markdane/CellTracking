//

if (is("Batch Mode")) {
	args = getArgument();
	if (args=="") exit ("no command line arguments");
	args = split(args, ",");
	movie_name = args[0];
	well_location_path = args[1];
	ligand1 = args[2];
	ligand2 = args[3];
	drug1 = args[4];
} else {
	well_location_path = "/Users/dane/Documents/CellTracking/IncuCyte/LI_I_L_035_01_1/A1_1";
	movie_name = "OSM_10";
	ligand1 = "OSM 10 ng/ml";
	drug1 = "";
}

//create subdirectory for the registered images
Output_dir = well_location_path + "/Analysis"+File.separator;
if (!File.exists(Output_dir)) File.makeDirectory(Output_dir);
if (!File.exists(Output_dir)) exit("Unable to read or create directory" + Output_dir);
      			
print("Opening stack in " +  well_location_path + "/P_Unreg");
run("Image Sequence...", "open=" + well_location_path + "/P_Unreg" + " file=(tif) sort");
setForegroundColor(255, 255, 255);
setBackgroundColor(10, 10, 10);
run("Label...", "format=Text starting=0 interval=1 x=1150 y=1000 font=18 text=[" + ligand1 + "]");
run("Label...", "format=Text starting=0 interval=1 x=1150 y=1018 font=18 text=[" + ligand2 + "]");
run("Label...", "format=Text starting=0 interval=1 x=1150 y=1036 font=18 text=[" + drug1 + "]");
run("AVI... ", "compression=JPEG frame=10 save=" + well_location_path + "/Analysis/" + movie_name + ".avi");

