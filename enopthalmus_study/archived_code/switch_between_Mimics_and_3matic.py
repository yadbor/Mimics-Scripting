# One script is used for both Mimics and 3-matic. 
# For that reason we have to check if we are in Mimics or in 3-matic 
try:
    import trimatic
except:
    in_3matic = False
else:
    in_3matic = True
SHARED_OBJ = "Union"
#If True we are in Mimics
if not in_3matic:
    # import required modules 
    import os
    import subprocess
    # Open Heart.mcs Mimics project
    path = r"C:\MedData\DemoFiles\Heart.mcs"
    mimics.file.open_project(path)
    # Find the masks of interest
    masks_names = ["LA","LV", "Aorta"]
    masks = []
    for m in masks_names:
        mask = mimics.data.masks.find(m)
        if mask:
            print("Mask "+ m +" is present.")
            masks.append(mask)
    # Unite masks
    if len(masks) == 3:
        un1 = mimics.segment.boolean_operations(masks[0],masks[1],"Unite")
        union = mimics.segment.boolean_operations(un1,masks[2],"Unite")
        union.name  = SHARED_OBJ
        # Create the Part of the Union mask
        union_part = mimics.segment.calculate_part(union)
        union_part.name = SHARED_OBJ
        #Export the Union Part in the location of the script
        root_path_of_script = os.path.split(os.path.abspath(__file__))[0]
        path_of_stl = os.path.join(root_path_of_script,union_part.name + ".stl")
        mimics.file.export_part(union_part,path_of_stl)
        with open(os.path.join(os.path.split(__file__)[0],"my_temp.txt"),"w") as f:
            f.write("File is created!\n")
        #Prepare to run 3-matic
        trimatic = mimics.file.get_path_to_3matic()
        command = trimatic 
        args = ("-run_script", __file__, path_of_stl,f.name)
        process = subprocess.Popen((command,) + args, shell=False, stdout=subprocess.PIPE)
        process.wait()
        with open(f.name,"r")as f:
            lines = f.readlines()
        os.remove(f.name)
        for i in range(2):
            mimics.file.import_stl(lines[i+1].strip())
    else:
        print("Please check if a mask is missing! Three masks are required.")
#If True we are in 3-matic
else:
    import sys
    path_of_stl = sys.argv[1]
    f = sys.argv[2]
    trimatic.import_part_stl(path_of_stl)
    part = trimatic.find_parts(SHARED_OBJ)
    if part:
        plane = trimatic.create_plane_fit(part[0])
        cut_parts = trimatic.cut(part[0],plane)
        exp = trimatic.export_stl_ascii(cut_parts,os.path.split(os.path.abspath(__file__))[0])
        with open(f,"a") as f:
            f.write(exp[0]+"\n")
            f.write(exp[1])
        print("To continue please close 3-matic!")
