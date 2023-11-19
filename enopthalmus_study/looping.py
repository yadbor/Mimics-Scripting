
users = ["David", "Ryan", "Rob", "Alan", "Dieter"]
str_users = ';'.join(users)
user = mimics.dialogs.question_box(message = "Choose who is making these measurements", title = "Choose Operator", buttons = str_users)


import os
#change directory to folder where script is stored
os.chdir("C:\\users\\bsantiag\\Documents\\BASILISC Scripts")            
#save file
mimics.file.save_project(filename= "{} Automated Segmentation".format(threshold))


f = r"C:\MedData\DemoFiles\Heart.mcs"
mimics.file.open_project(filename=f)

mimics.file.open_project(
    filename=f,
    read_only_mode=True  # after this call, the project can only
                         # be Saved-As under another file name
    )
    
