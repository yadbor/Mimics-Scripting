import mimics
import os
import numpy as np
import nibabel as nib

def get_affine_from_mimics_mask(mask,image):
    """
    Create an affine transformation matrix for the NIfTI image based on the Mimics mask.
    """
    # Get the voxel buffer shape (dimensions)
    voxel_buffer = mask.get_voxel_buffer()
    shape = voxel_buffer.shape

    # Get voxel center coordinates
    p0 = image.get_voxel_center([0, 0, 0])
    x = image.get_voxel_center([shape[0] - 1, 0, 0])
    y = image.get_voxel_center([0, shape[1] - 1, 0])
    z = image.get_voxel_center([0, 0, shape[2] - 1])

    # Calculate the directions and voxel size
    voxel_size = [
        (x[0] - p0[0]) / (shape[0] - 1),
        (y[1] - p0[1]) / (shape[1] - 1),
        (z[2] - p0[2]) / (shape[2] - 1)
    ]

    # Create the affine transformation matrix
    affine = np.eye(4)
    affine[0, 0] = voxel_size[0]
    affine[1, 1] = voxel_size[1]
    affine[2, 2] = voxel_size[2]
    affine[:3, 3] = p0

    return affine

def save_mask_as_nifti(mask, image, output_dir):
    # Get the mask data using get_voxel_buffer
    voxel_buffer = mask.get_voxel_buffer()
    mask_array = np.array(voxel_buffer, dtype=np.uint8)

    # Get the affine transformation matrix for the mask
    affine = get_affine_from_mimics_mask(mask,image)

    # Create the NIfTI image
    nifti_image = nib.Nifti1Image(mask_array, affine)

    # Define the output filename
    output_filename = os.path.join(output_dir, f"{mask.name}.nii.gz")

    # Save the NIfTI image
    nib.save(nifti_image, output_filename)
    print(f"Saved mask '{mask.name}' as NIfTI to '{output_filename}'")

def main():
    # Define the Mimics project file path
    project_file_path = f'C:\MedData\DemoFiles\masked.mcs'

    # Open the Mimics project
    mimics.file.open_project(project_file_path)

    #aquire the first image so we can get the affine transformation matrix
    images = mimics.data.images
    if not images:
        print("No images found in the project.")
        return

    image = images[0]

    # Create the output directory within the project file path
    project_dir = os.path.dirname(project_file_path)
    project_name = os.path.splitext(os.path.basename(project_file_path))[0]
    output_dir = os.path.join(project_dir, f"mimics_export_{project_name}")
    dicom_dir = os.path.join(output_dir, f"DICOM")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(dicom_dir):
        os.makedirs(dicom_dir)
    
    # Loop through all masks in the project
    masks = mimics.data.masks
    for mask in masks:
        save_mask_as_nifti(mask, image, output_dir)
        mask.visible = False

    # Save the Images (now tha the masks are off)
    mimics.file.export_dicom(path=dicom_dir, filename_prefix='dicom')

    # Close the Mimics project
    mimics.file.close_project()

if __name__ == "__main__":
    main()

