function MultiStackRegBeacons(beacon, data_dir) {
	run("Image Sequence...", "open=" + data_dir + "Beacon-" + beacon + "/Red_Unreg starting=2 sort");
	run("Image Sequence...", "open=" + data_dir + "Beacon-" + beacon + "/GFP_Unreg starting=2 sort");
	run("MultiStackReg", "stack_1=Red_Unreg action_1=Align file_1=[] stack_2=GFP_Unreg action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body]");
	selectWindow("Red_Unreg");
	run("Image Sequence... ", "format=TIFF use save=" + data_dir + "Beacon-" + beacon + "/Red_Reg");
	close();
	selectWindow("GFP_Unreg");
	run("Image Sequence... ", "format=TIFF use save=" + data_dir + "Beacon-" + beacon + "/GFP_Reg");
	close();
}

data_dir = "/graylab/share/dane/HeLa_1/"

for (beacon = 1; beacon < 2; beacon++) {
        MultiStackRegBeacons(beacon, data_dir);
}
