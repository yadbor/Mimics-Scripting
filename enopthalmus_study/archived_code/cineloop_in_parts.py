# import modules
import mimics
import logging
import os
import time
import datetime
try:
    import cv2   # First need to install Open CV package for Python. Type pip install opencv-python in your cmd
except ImportError as ie:
    print("================================================================")
    print("=== The 3rd party Python package OpenCV is not installed! ===")
    print("=== To install it, use pip install opencv-python in your cmd! ===")
    print("================================================================")
    raise
##########################################################################################
##########################################################################################
# config
# Here is the configuration area of the script.
# List the names of the anatomy of the Left Heart that you want to visualise in the video.
# In the provided script it is assumed that Left Atrium is LA, Left Ventricle is LV and Aorta is Aorta
LEFT_HEART = ["LA", "LV", "Aorta"]
# List the names of the anatomy of the Right Heart that you want to visualise in the video.
# In the provided script it is assumed that Right Atrium is RA, Right Ventricle is RV and Pulmonary Artery is PA
RIGHT_HEART = ["RA", "RV", "PA"]
# Here you select with part of the heart to visualise in the video.
HEART = LEFT_HEART + RIGHT_HEART
# The name (prefix) of the folder where the screenshots and the video will be saved. The suffix of the folder is the
# name of the Mimics file from where you run the script.
TARGET_FOLDER_NAME = "Output"
# Name of the file of the output video
VIDEO_NAME = "Cineloop_in_parts"
#Frames per second of the video. You can change this value to make the heart beating faster or slower.
FRAMES = 12
# If you want to visualise a logo in your video, the logo should be placed in the same folder as the script.
# If you do not want a logo to be visualised, just leave the script empty. Example: LOGO = ""
LOGO = "mat_logo.jpg"
# To calculate first the 3D Parts from the segmentation masks, set the value of the variable below to True.
# In case you have the Parts already calculated, you can leave it to False
CALCULATE_PARTS_FROM_MASKS = False
##########################################################################################
##########################################################################################
# settings
# Please do NOT modify this value
VIEW_OF_INTEREST = "3D"
IMAGE_FILE_SUFFIX = "jpg"
FONT = cv2.FONT_HERSHEY_DUPLEX
# Here you can change the size of the text that appears in the bottom right side of the video.
# It indicates the name of the cardiac phase that is visualised.
FONT_SIZE = 0.7
# Here you can change the colour (in RGB [0,255]) of the text that appears in the bottom right side of the video.
# It indicates the name of the cardiac phase that is visualised.
FONT_COLOR = (0, 0, 0)
LINE_TYPE = 1
##########################################################################################


class Phase:
    def __init__(self, image_obj):
        self.mimics_obj = image_obj
        self.patient_name = self.get_patient_name()
        self.image_name = self.mimics_obj.name
        self.guid = self.mimics_obj.guid
        self.linked_objs = self.mimics_obj.linked_objects

    def get_patient_name(self) -> str:
        info = self.mimics_obj.get_image_information()
        return info.patient_name
##########################################################################################


def get_correct_phases() -> dict:
    phases = {}
    for im in mimics.data.images:
        size = im.physical_dimensions
        phases[im.guid] = size
    return mimics.data.images

def get_initial_status() -> tuple:
    overlay = mimics.view.is_overlay_enabled()
    active_is = mimics.data.images.get_active()
    visible_objs = []
    for obj in mimics.data.objects:
        if not isinstance(obj, mimics.ImageData) and obj.visible is True:
            visible_objs.append(obj.guid)
    return overlay, active_is, visible_objs

def set_initial_status(overlay, active_is, visible_objs) -> None:
    if not overlay:
        mimics.view.disable_overlay()
    mimics.data.images.set_active(active_is)
    for obj in mimics.data.objects:
        if not isinstance(obj, mimics.ImageData):
            obj.visible = False
            if obj.guid in visible_objs:
                obj.visible = True
    return


def reset_views() -> None:
    # enable Overlay that will make our life easier
    if not mimics.view.is_overlay_enabled():
        mimics.view.enable_overlay()
    # Get current active IS
    currently_active_is = mimics.data.images.get_active()
    # switch to the Standard single IS layout
    mimics.view.set_layout(layout_name=mimics.Layouts.Standard, images=currently_active_is)
    # switch off 3D mask preview if ON
    mimics.view.disable_mask_3d_preview()
    # reset the objects that are visualised
    for obj in mimics.data.objects:
        if not isinstance(obj, mimics.ImageData):
                obj.visible = False
    return None

def calculate_parts() -> None:
    # show the correct objects
    for im in mimics.data.images:
        if im.linked_objects is not None:
            mimics.data.images.set_active(im)
            for objs in im.linked_objects:
                if isinstance(objs, mimics.segment.Mask) and objs.name in HEART:
                    prt = mimics.segment.calculate_part(objs)
                    prt.name = objs.name
                    prt.image = objs.image
                    prt.color = objs.color
                    prt.visible = False
    return None


def prepare_objects(image) -> None:
    # show the correct objects
    phase_visible_objs = []
    for objs in image.linked_objs:
        if isinstance(objs, mimics.Part) and objs.name in HEART:
            objs.visible = True
            phase_visible_objs.append(objs.name)
    logging.info("Phase name:")
    logging.info(image.image_name)
    logging.info("Objects that will be included in the video:")
    logging.info(phase_visible_objs)
    return None


def get_camera() -> tuple:
    view = None
    for v in mimics.data.views:
        if v.type == VIEW_OF_INTEREST:
            cam = v.get_camera()
            settings = cam.get_settings()
            view = v
    return settings, view


def set_camera(settings) -> None:
    for v in mimics.data.views:
        if v.type == VIEW_OF_INTEREST:
            cam = v.get_camera()
            cam.set_settings(settings)
    return None


def capture_frame(view, counter, directory,phase_name) -> bool:
    filename = os.path.join(directory, str(counter) + "." + IMAGE_FILE_SUFFIX)
    try:
        mimics.file.export_view(filename, view)
    except RuntimeError:
        print("It was not possible to export the frame.")
        return False
    img = cv2.imread(filename, 1)
    height, width,_ = img.shape
    txt_size = cv2.getTextSize(phase_name, FONT, FONT_SIZE, LINE_TYPE)
    text_location = (width-txt_size[0][0]-10, height-txt_size[0][1]-10)
    cv2.putText(img, phase_name, text_location, FONT, FONT_SIZE, FONT_COLOR, LINE_TYPE)
    if get_logo(LOGO):
        logo = cv2.imread(get_logo(LOGO), 1)
        rows, cols, channels = logo.shape
    # logo_location = img[rows:-1, cols:-1]
    # print(len(logo_location))
        img[50:rows+50, width-cols-50:-50] = logo
    cv2.imwrite(filename, img)
    return True


def create_dir() -> tuple:
    p = mimics.file.get_project_information()
    path = p.project_path
    name_of_file_with_extension = os.path.basename(path)
    name_of_file = os.path.splitext(name_of_file_with_extension)[0]
    target_folder = TARGET_FOLDER_NAME + "_" + name_of_file
    root = os.path.dirname(__file__)
    new_dir = False
    if not os.path.isdir(os.path.join(root, target_folder)):
        try:
            os.makedirs(os.path.join(root, target_folder))
            new_dir = True
            logging.info("A new folder Output is created")
        except PermissionError:
            logging.info("=== The script does not have the permissions to create a folder here! ===")
            raise
    else:
        logging.info("The Output folder already exists.")
    return os.path.join(root, target_folder), new_dir



def clean_dir(directory) -> None:
    filelist = os.listdir(directory)
    for file in filelist:
        os.remove(os.path.join(directory, file))
    logging.info("The Output folder is cleaned.")
    return None


def get_logo(logo_name)-> str:
    root = os.path.dirname(__file__)
    logo_full_path = os.path.join(root,logo_name)
    if os.path.isfile(logo_full_path):
        return logo_full_path


def create_the_video(directory)-> None:
    image_folder = directory
    video_name = VIDEO_NAME
    images = [img for img in os.listdir(image_folder) if img.endswith("."+IMAGE_FILE_SUFFIX)]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, _ = frame.shape
    fourcc = cv2.VideoWriter_fourcc(*'WMV1')
    video = cv2.VideoWriter(os.path.join(image_folder, video_name + ".avi"), fourcc, FRAMES, (width, height))
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    cv2.destroyAllWindows()
    video.release()
    return None


if __name__ == "__main__":
        log_file = os.path.join(os.path.dirname(__file__), 'report.log')
        logging.basicConfig(filename=log_file, filemode='w', level=logging.INFO)
        with open(log_file, "w")as f:
            f.truncate()
        logging.info('Cineloop in Parts is started. Mimics UI will be frozen for a few seconds...')
        start = time.time()
        now = datetime.datetime.now()
        logging.info("Date and Time:{}".format(now.now()))
        mimics.disable_update_gui()
        phases = get_correct_phases()
        logging.info("The following image sets will be included in the script: ")
        t = []
        for p in phases:
            t.append(p.name)
        logging.info(t)
        output_dir, is_dir_new = create_dir()
        if not is_dir_new:
            clean_dir(output_dir)
        if CALCULATE_PARTS_FROM_MASKS:
            calculate_parts()
            mimics.dialogs.question_box(message='Set the desired camera view in the 3D viewport. Then click OK to continue..', buttons='OK', title='Set the camera', ui_blocking=False)
        cl_phases = []
        overlay, active_is, visible_objs = get_initial_status()
        settings, view = get_camera()
        # print(settings, view)
        for c, p in enumerate(phases):
            ph = Phase(p)
            cl_phases.append(ph)
            reset_views()
            prepare_objects(ph)
            set_camera(settings)
            export = True
            if view:
                if not capture_frame(view, c, output_dir, ph.image_name):
                    export = False
            else:
                logging.info('The 3D view is not visible in Mimics. Please select a layout where the 3D view is visible.')
                export = False
        set_initial_status(overlay, active_is, visible_objs)

        if export:
            create_the_video(output_dir)
            logging.info('The video is exported.')
        mimics.enable_update_gui()
        logging.info("The script is terminated.")
