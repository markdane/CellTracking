run("Image Sequence...", "open=/Users/dane/Documents/CellTracking/Sample2chPNGs/temp_images_files/Beacon-3/Red_Unreg starting=2 sort");
run("Image Sequence...", "open=/Users/dane/Documents/CellTracking/Sample2chPNGs/temp_images_files/Beacon-3/GFP_Unreg starting=2 sort");
run("MultiStackReg", "stack_1=Red_Unreg action_1=Align file_1=[] stack_2=GFP_Unreg action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body]");
selectWindow("Red_Unreg");
run("Image Sequence... ", "format=TIFF use save=/Users/dane/Documents/CellTracking/Sample2chPNGs/temp_images_files/Beacon-3/Red_Reg");
close();
selectWindow("GFP_Unreg");
run("Image Sequence... ", "format=TIFF use save=/Users/dane/Documents/CellTracking/Sample2chPNGs/temp_images_files/Beacon-3/GFP_Reg");
