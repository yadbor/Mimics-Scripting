import mimics
import numpy as np

# parameters
dilateAirPixels = 2  # Number of pixels that the air spaces are dilated to create bone shell
closingPixels = 3  # Number of pixels used for the closing operation

def make_air():
    air = mimics.segment.create_mask()
    air.name = 'Air'
    mimics.segment.threshold(mask=air, threshold_min=0,
                             threshold_max=350)  # thresholds are set in gray values
    return air
    
def make_soft_air():

    airSofterThreshold = mimics.segment.create_mask()  # We make this mask to distract air particles from final segmentation
    airSofterThreshold.name = 'Air softer threshold'
    mimics.segment.threshold(mask=airSofterThreshold, threshold_min=0,
                             threshold_max=400)  
    return airSofterThreshold

def unite(mask_a, mask_b):
	return mimics.segment.boolean_operations(mask_a, mask_b, operation='Unite')

def minus(mask_a, mask_b):
	return mimics.segment.boolean_operations(mask_a, mask_b, operation='Minus')

def intersect(mask_a, mask_b):
	return mimics.segment.boolean_operations(mask_a, mask_b, operation='Intersect')

def morphology_op(mask, pixels, op = 'Dilate', name = None):
	return mimics.segment.morphology_operations(input_mask=mask, 
												operation=op, number_of_pixels=pixels, connectivity=26, 
												target_mask_name=name, limited_to_mask=None)

def make_shells(air, bone, min_width = 170, dilate_pixels, close_pixels):
    
	largest = mimics.data.objects.duplicate(air)
	mimics.segment.keep_largest(largest)
    # How big is background (the largest part of air)?
    # If it is less than min_width in X then add the second largest part.
	# This is found from largest(air - largest(air)) i.e. the largest of what is left after taking away the largest
	largest_X = mimics.measure.get_bounding_box(largest, first_axis=[1, 0, 0], second_axis=[0, 1, 0]).first_vector[0]
    if abs(largest_X) < min_width:
		second_largest = minus(air, largest)
		mimics.segment.keep_largest(second_largest)
		background = unite(largest, second_largest)
	else:
		background = largest

	air_internal = minus(air, background)
	air_dilated = morphology_op(air_internal, dilate_pixels, op = 'Dilate')
	bone_shell = minus(air_dilated, air_internal)

	temporary = unite(bone, bone_shell)
	bone_and_shell = morphology_op(temporary, close_pixels, op = 'Close')
	
    mimics.data.objects.delete(temporary)
    boneAndShell.name = 'Bone and shell'
    boneAndShell.visible = False

	
	background1 = mimics.data.objects.duplicate(air)
    background1.name = 'Background 1'
    background1 = mimics.segment.keep_largest(background1)

    bboxBackground1 = mimics.measure.get_bounding_box(background1, first_axis=[1, 0, 0], second_axis=[0, 1, 0])

    if abs(bboxBackground1.first_vector[0]) < min_width:
        background2 = mimics.segment.boolean_operations(air, background1, operation='Minus')
        background2 = mimics.segment.keep_largest(background2)
        background2.name = 'Background 2'
        background = mimics.segment.boolean_operations(background2, background1, operation='Unite')
    else:
        background = background1

    background.name = 'Background'
    airInSkull = mimics.segment.boolean_operations(air, background, operation='Minus')
    airInSkull.name = 'Air in skull'
    airInSkullDilated = morphology_op(airInSkull, dilate_pixels, op = 'Dilate')
    airInSkullDilated.name = 'Air in skull dilated'

    boneShell = mimics.segment.boolean_operations(airInSkullDilated, airInSkull, operation='Minus')
    boneShell.name = 'Bone shell'

	
    temporary = mimics.segment.boolean_operations(bone, boneShell, operation='Unite')
	boneAndShell = morphology_op(temporary, close_pixels, op = 'Close')

    mimics.data.objects.delete(temporary)
    boneAndShell.name = 'Bone and shell'
    boneAndShell.visible = False

    return bone_and_shell


	largest = mimics.segment.keep_largest(background)
	
	
	
	
#Initial warning
warningMessage = 'This algorithm should only be used for models that are not used for surgery \n ' \
          'on the orbital floor or nasal area. Since detail may be lost in these areas.'
mimics.dialogs.message_box(warningMessage, title='WARNING', ui_blocking=True)

# Create basic skull using the standard threshold to select boundarypoint
print('Starting!')
bone = mimics.segment.create_mask()
bone.name = 'Bone'
mimics.segment.threshold(mask=bone, 
                         threshold_min=1250,
                         threshold_max=4095)  # thresholds are set in gray values
mimics.segment.keep_largest(bone)  # Important to do this otherwise floating pixels will be merged during closing opereation!

whatQuestion = 'What do you want to improve?'
whatAns = mimics.dialogs.question_box(whatQuestion, buttons='orbital floor; cheek; both', title=None, ui_blocking=True)

if whatAns == 'orbital floor':
    sideQuestion = 'Both sides or one side?'
    sideAns = mimics.dialogs.question_box(whatQuestion, buttons='both sides; one side', title=None, ui_blocking=True)
    if sideAns == 'one side':
        point_question = 'Indicate point below the middle of the socket'
    else:
        point_question = 'Indicate point below the middle of the right eye socket'
        
    boundaryPoint = mimics.analyze.indicate_point(message=point_question, show_message_box=True, confirm=False, title=None)

    print('Creating boundary box')
    # First we create a big box because the air needs to be distractable by using 'keep_largest' on the right and left side
    origin = np.array([200, 12, 5])
    vectors = ([500, 0, 0], [0, 72 0], [0, 0, 65])
    bboxBig = mimics.BoundingBox3d(origin=boundaryPoint - offset, first_vector=vectors[0], second_vector=vectors[1], third_vector=vectors[2])

    print('Creating bone shells around air spaces')
    air = make_air()
    air = mimics.segment.crop_mask(air, bboxBig)
    bone_shell = make_shells(air, min_width = 170, dilate_pixels = dilate_air_pixels, close_pixels = close_air_pixles)

    # A smaller box is made to limit the adaptation of the data
    if sideAns == 'one side':
        offset = np.array([25, 10, 5])
        vectors = ([50, 0, 0], [0, 75 0], [0, 0, 65])
    else:
        offset = np.array([25, 10, 5])
        vectors = ([120, 0, 0], [0, 75 0], [0, 0, 65])

    bboxSmall = mimics.BoundingBox3d(origin=boundaryPoint - offset, first_vector=vectors[0], second_vector=vectors[1], third_vector=vectors[2])

    bone_shell_cropped = mimics.segment.crop_mask(bone_shell, bboxSmall)
    bone_shell_cropped.name = 'Bone and shell cropped'
    bone_shell_cropped.visible = False

    air_softer = make_soft_air()

    # Creating the final mask by joining bone and the shell around the air
    finalMaskOrbital = mimics.segment.boolean_operations(bone_shell_cropped, air_softer, operation='Minus')
    mimics.segment.keep_largest(finalMaskOrbital)
    finalMaskOrbital.name = 'Orbit, unite me!'
    
    # Deleting un-needed masks

elif whatAns == 'cheek':
    sideQuestion = 'Both sides or one side?'
    sideAns = mimics.dialogs.question_box(whatQuestion, buttons='both sides; one side', title=None, ui_blocking=True)
    if sideAns == 'one side':
        point_question = 'Indicate point below the middle of the maxillar sinus hole'
    else:
        point_question = 'Indicate point below the middle of the right maxillar sinus hole'
        
    boundaryPoint = mimics.analyze.indicate_point(message=point_question, show_message_box=True, confirm=False, title=None)

    print('Creating boundary box')
    # First we create a big box because the air needs to be distractable by using 'keep_largest' on the right and left side
    ## For 'floor' this is origin = boundaryPoint - [200, 12, 5] and vectors = ([500, 0, 0], [0, 72 0], [0, 0, 65])
    bboxBig = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 200, boundaryPoint[1] - 15, boundaryPoint[2]],
                                   first_vector=[500, 0, 0], second_vector=[0, 30, 0], third_vector=[0, 0, 35])

    print('Creating bone shells around air spaces')
    air = make_air()
    air = mimics.segment.crop_mask(air, bboxBig)
    bone_shell = make_shells(air, min_width = 170, dilate_pixels = dilate_air_pixels, close_pixels = close_air_pixles)

    # A smaller box is made to limit the adaptation of the data
    if sideAns == 'one side':
        offset = np.array([20, 5, 5])
        vectors = ([30, 0, 0], [0, 13 0], [0, 0, 25])
    else:
        offset = np.array([20, 5, 5])
        vectors = ([85, 0, 0], [0, 16 0], [0, 0, 35])

    bboxSmall = mimics.BoundingBox3d(origin=boundaryPoint - offset, first_vector=vectors[0], second_vector=vectors[1], third_vector=vectors[2])

    bone_shell_cropped = mimics.segment.crop_mask(bone_shell, bboxSmall)
    bone_shell_cropped.name = 'Bone and shell cropped'
    bone_shell_cropped.visible = False

    air_softer = make_soft_air()

    finalMaskCheek = mimics.segment.boolean_operations(bone_shell_cropped, air_softer, operation='Minus')
    mimics.segment.keep_largest(finalMaskCheek)
    finalMaskCheek.name = 'Cheek, unite me!'
    
    # Deleting un-needed masks
    
else: # must be 'both'
    # Two boundarytPoints in this case as are doing both, but only on one side
    boundaryPointOrbital = mimics.analyze.indicate_point(message='Indicate point below the middle of the right eye-socket',
                                                         show_message_box=True,
                                                         confirm=False, title='Orbital point')

    boundaryPointCheek = mimics.analyze.indicate_point(message='Indicate point below the middle of the right maxillar sinus hole',
                                                     show_message_box=True,
                                                     confirm=False, title='Cheek point')

    print('Creating boundary box')
    origin = np.array([200, 12, 5])
    vectors = ([500, 0, 0], [0, 72 0], [0, 0, 65])
    bboxBigOrbital = mimics.BoundingBox3d(origin=boundaryPoint - offset, first_vector=vectors[0], second_vector=vectors[1], third_vector=vectors[2])
    
    origin = np.array([200, 15, 2])
    vectors = ([500, 0, 0], [0, 30 0], [0, 0, 35])
    bboxBigCheek = mimics.BoundingBox3d(origin=boundaryPoint - offset, first_vector=vectors[0], second_vector=vectors[1], third_vector=vectors[2])

    print('Creating bone shells around air spaces')
    air = make_air()
    air_orbital = mimics.data.objects.duplicate(air)
    air_orbital = mimics.segment.crop_mask(air_orbital, bboxBigOrbital)

    air_cheek = mimics.data.objects.duplicate(air)
    air_cheek = mimics.segment.crop_mask(air_cheek, bboxBigCheek)

    air=mimics.segment.boolean_operations(airCheek, airOrbital, operation='Unite')
    air.name = 'Air'

    bone_and_shell = make_sheels(air, min_width = 170, dilate_pixels, close_pixels)
    
    boneAndShellOrbital = mimics.data.objects.duplicate(boneAndShell)
    boneAndShellCheek = mimics.data.objects.duplicate(boneAndShell)
    mimics.data.objects.delete(boneAndShell)

    bboxSmallOrbital = mimics.BoundingBox3d(origin=[boundaryPointOrbital[0] - 25, boundaryPointOrbital[1] - 10, boundaryPointOrbital[2] - 5],
                                     first_vector=[120, 0, 0], second_vector=[0, 65, 0], third_vector=[0, 0, 50])
    bboxSmallCheek = mimics.BoundingBox3d(origin=[boundaryPointCheek[0] - 20, boundaryPointCheek[1]-5, boundaryPointCheek[2]-5],
                            first_vector=[85, 0, 0], second_vector=[0, 16, 0], third_vector=[0, 0, 35])

    boneAndShellOrbitalCropped = mimics.segment.crop_mask(boneAndShellOrbital, bboxSmallOrbital)
    boneAndShellCheekCropped = mimics.segment.crop_mask(boneAndShellCheek, bboxSmallCheek)
    
    boneAndShellCropped = mimics.segment.boolean_operations(boneAndShellOrbitalCropped, boneAndShellCheekCropped, operation='Unite')
    boneAndShellCropped.name = 'Bone and shell cropped'
    boneAndShellCropped.visible = False
    
    mimics.data.objects.delete(boneAndShellOrbital)
    mimics.data.objects.delete(boneAndShellCheek)

    air_softer = make_soft_air()
    
    finalMaskOrbital = mimics.segment.boolean_operations(boneAndShellCropped, air_softer, operation='Minus')
    mimics.segment.keep_largest(finalMaskOrbital)
    finalMaskOrbital.name = 'Unite me!'
    # Delete un-needed masks (all but final one?)

print('Done')
