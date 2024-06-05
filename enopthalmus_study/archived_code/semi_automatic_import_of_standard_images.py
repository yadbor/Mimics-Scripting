# Tutorial for importing standard images
input_dir=r"C:\MedData\DemoFiles\BMP_Leg"
mimics.dialogs.set_predefined_answer(mimics.dialogs.dialog_id.CHANGE_ORIENTATION, "RAB")
mimics.file.import_standard_images(source_folder=input_dir,xy_resolution=1,z_resolution=1,patient_name="MimMat")    


