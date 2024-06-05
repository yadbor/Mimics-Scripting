#Tutorial for automatic import of dicom images
input_dir=r"C:\MedData\DemoFiles\DICOM_Airway"
mimics.dialogs.set_predefined_answer("ChangeOrientation", "default")
mimics.file.import_dicom_images(source_folder=input_dir)   
mimics.file.anonymize_active_image()

