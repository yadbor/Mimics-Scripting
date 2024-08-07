import mimics

#Initial warning
warningMessage = 'This algorithm should only be used for models that are not used for surgery \n ' \
          'on the orbital floor or nasal area. Since detail may be lost in these areas.'
mimics.dialogs.message_box(warningMessage, title='WARNING', ui_blocking=True)

whatQuestion = 'What do you want to improve?'
whatAns = mimics.dialogs.question_box(whatQuestion, buttons='orbital floor; cheek; both', title=None, ui_blocking=True)

if whatAns == 'orbital floor':
    sideQuestion = 'Both sides or one side?'
    sideAns = mimics.dialogs.question_box(whatQuestion, buttons='both sides; one side', title=None, ui_blocking=True)

    # parameters
    dilateAirPixels = 2  # Number of pixels that the air spaces are dilated to create bone shell
    closingPixels = 3  # Number of pixels used for the closing operation

    # Create basic skull using the standard threshold to select boundarypoint
    print('Creating basic skull to select boundarypoint, please go grab a coffee')
    bone = mimics.segment.create_mask()
    bone.name = 'Bone'
    mimics.segment.threshold(mask=bone, threshold_min=1250,
                             threshold_max=4095)  # thresholds are set in gray values
    mimics.segment.keep_largest(
        bone)  # Important to do this otherwise floating pixels will be merged during closing opereation!

    if sideAns == 'one side':
        boundaryPoint = mimics.analyze.indicate_point(message='Indicate point below the middle of the socket',
                                                      show_message_box=True,
                                                      confirm=False, title=None)
    else:
        boundaryPoint = mimics.analyze.indicate_point(message='Indicate point below the middle of the right eye socket',
                                                  show_message_box=True,
                                                  confirm=False, title=None)

    print('Creating boundary box')
    bboxBig = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 200, boundaryPoint[1] - 12, boundaryPoint[2] - 5],
                                   first_vector=[500, 0, 0], second_vector=[0, 72, 0], third_vector=[0, 0, 65])

    print('Creating bone shells around air spaces')
    air = mimics.segment.create_mask()
    air.name = 'Air'
    mimics.segment.threshold(mask=air, threshold_min=0,
                             threshold_max=350)  # thresholds are set in gray values
    air = mimics.segment.crop_mask(air, bboxBig)

    background1 = mimics.data.objects.duplicate(air)
    background1.name = 'Background 1'
    background1 = mimics.segment.keep_largest(background1)

    #This piece of code prevents that too much air is deleted in case background1 already captured all the
    #background air
    bbxBackground1 = mimics.measure.get_bounding_box(background1, first_axis=[1, 0, 0], second_axis=[0, 1, 0])

    if abs(bbxBackground1.first_vector[0])<170:
        background2 = mimics.segment.boolean_operations(air, background1, operation='Minus')
        background2 = mimics.segment.keep_largest(background2)
        background2.name = 'Background 2'
        background = mimics.segment.boolean_operations(background2, background1, operation='Unite')
    else:
        background = background1
    background.name = 'Background'
    airInSkull = mimics.segment.boolean_operations(air, background, operation='Minus')
    airInSkull.name = 'Air in skull'
    airInSkullDilated = mimics.segment.morphology_operations(input_mask=airInSkull, operation='Dilate',
                                                             number_of_pixels=dilateAirPixels,
                                                             connectivity=26,
                                                             limited_to_mask=None)
    airInSkullDilated.name = 'Air in skull dilated'


    boneShell = mimics.segment.boolean_operations(airInSkullDilated, airInSkull, operation='Minus')
    boneShell.name = 'Bone shell'

    temporary = mimics.segment.boolean_operations(bone, boneShell, operation='Unite')
    boneAndShell = mimics.segment.morphology_operations(input_mask=temporary, operation='Close',
                                                        number_of_pixels=closingPixels,
                                                        connectivity=26,
                                                        limited_to_mask=None)
    mimics.data.objects.delete(temporary)
    boneAndShell.name = 'Bone and shell'
    boneAndShell.visible = False
    if sideAns == 'one side':
        bboxSmall = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 25, boundaryPoint[1] - 10, boundaryPoint[2] - 5],
                                         first_vector=[50, 0, 0], second_vector=[0, 75, 0], third_vector=[0, 0, 65])
    else:
        bboxSmall = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 25, boundaryPoint[1] - 10, boundaryPoint[2] - 5],
                                     first_vector=[120, 0, 0], second_vector=[0, 75, 0], third_vector=[0, 0, 65])
    boneAndShellCropped = mimics.segment.crop_mask(boneAndShell, bboxSmall)
    boneAndShellCropped.name = 'Bone and shell cropped'
    boneAndShellCropped.visible = False

    # Creating the final mask by joining bone and the shell around the air
    airSofterThreshold = mimics.segment.create_mask()  # We make this mask to distract air particles from final segmentation
    airSofterThreshold.name = 'Air softer threshold'
    mimics.segment.threshold(mask=airSofterThreshold, threshold_min=0,
                             threshold_max=400)  # thresholds are set in gray values
    finalMaskOrbital = mimics.segment.boolean_operations(boneAndShellCropped, airSofterThreshold, operation='Minus')
    mimics.segment.keep_largest(finalMaskOrbital)
    finalMaskOrbital.name = 'Orbit, unite me!'


    # # Deleting shit
    # mimics.data.objects.delete(airSofterThreshold)

    # # Deleting all the masks except for the final one
    # objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
    #                'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background']
    # for maskName in objectsList:
    #     mask = mimics.data.masks.find(maskName)
    #     if mask is not None:
    #         mask.visible = False

    # # Deleting all the masks except for the final one
    # objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
    #                'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background', 'Air softer threshold']
    # for maskName in objectsList:
    #     mask = mimics.data.masks.find(maskName)
    #     if mask is not None:
    #         mimics.data.objects.delete(mask)


elif whatAns == 'cheek':
    sideQuestion = 'Both sides or one side?'
    sideAns = mimics.dialogs.question_box(whatQuestion, buttons='both sides; one side', title=None, ui_blocking=True)

    # parameters
    dilateAirPixels = 2  # Number of pixels that the air spaces are dilated to create bone shell
    closingPixels = 3  # Number of pixels used for the closing operation

    # Create basic skull to select boundarypoint
    bone = mimics.segment.create_mask()
    bone.name = 'Bone'
    mimics.segment.threshold(mask=bone, threshold_min=1250,
                             threshold_max=4095)  # thresholds are set in gray values
    mimics.segment.keep_largest(bone)
    # part = mimics.segment.calculate_part(bone)
    if sideAns == 'one side':
        boundaryPoint = mimics.analyze.indicate_point(message='Indicate point below the middle of the maxillar sinus hole',
                                                      show_message_box=True,
                                                      confirm=False, title=None)
    else:
        boundaryPoint = mimics.analyze.indicate_point(message='Indicate point below the middle of the right maxillar sinus hole',
                                                  show_message_box=True,
                                                  confirm=False, title=None)


    print('Creating boundary box')
    # First we create a big box because the air needs to be distractable by using 'keep_largest' on the right and left side
    bboxBig = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 200, boundaryPoint[1] - 15, boundaryPoint[2]],
                                   first_vector=[500, 0, 0], second_vector=[0, 30, 0], third_vector=[0, 0, 35])

    print('Creating bone shells around air spaces')
    air = mimics.segment.create_mask()
    air.name = 'Air'
    mimics.segment.threshold(mask=air, threshold_min=0,
                             threshold_max=350)  # thresholds are set in gray values
    air = mimics.segment.crop_mask(air, bboxBig)

    background1 = mimics.data.objects.duplicate(air)
    background1.name = 'Background 1'
    background1 = mimics.segment.keep_largest(background1)

    #This piece of code prevents that too much air is deleted in case background1 already captured all the
    #background air
    bbxBackground1 = mimics.measure.get_bounding_box(background1, first_axis=[1, 0, 0], second_axis=[0, 1, 0])

    if abs(bbxBackground1.first_vector[0])<170:
        background2 = mimics.segment.boolean_operations(air, background1, operation='Minus')
        background2 = mimics.segment.keep_largest(background2)
        background2.name = 'Background 2'
        background = mimics.segment.boolean_operations(background2, background1, operation='Unite')
    else:
        background = background1
    background.name = 'Background'
    airInSkull = mimics.segment.boolean_operations(air, background, operation='Minus')
    airInSkull.name = 'Air in skull'
    airInSkullDilated = mimics.segment.morphology_operations(input_mask=airInSkull, operation='Dilate',
                                                             number_of_pixels=dilateAirPixels,
                                                             connectivity=26,
                                                             limited_to_mask=None)
    airInSkullDilated.name = 'Air in skull dilated'


    boneShell = mimics.segment.boolean_operations(airInSkullDilated, airInSkull, operation='Minus')
    boneShell.name = 'Bone shell'

    temporary = mimics.segment.boolean_operations(bone, boneShell, operation='Unite')
    boneAndShell = mimics.segment.morphology_operations(input_mask=temporary, operation='Close',
                                                        number_of_pixels=closingPixels,
                                                        connectivity=26,
                                                        limited_to_mask=None)
    mimics.data.objects.delete(temporary)
    boneAndShell.visible = False
    boneAndShell.name = 'Bone and shell'

    # A smaller box is made to limit the adaptation of the data
    if sideAns == 'one side':
        bboxSmall = mimics.BoundingBox3d(origin=[boundaryPoint[0] - 20, boundaryPoint[1] - 5, boundaryPoint[2] - 5],
                                         first_vector=[30, 0, 0], second_vector=[0, 13, 0], third_vector=[0, 0, 25])
    else:
        bboxSmall = mimics.BoundingBox3d(
            origin=[boundaryPointCheek[0] - 20, boundaryPointCheek[1] - 5, boundaryPointCheek[2] - 5],
            first_vector=[85, 0, 0], second_vector=[0, 16, 0], third_vector=[0, 0, 35])

    boneAndShellCropped = mimics.segment.crop_mask(boneAndShell, bboxSmall)
    boneAndShellCropped.name = 'Bone and shell cropped'
    boneAndShellCropped.visible = False

    # finalMask = mimics.segment.boolean_operations(bone, boneAndShell, operation='Unite')
    airSofterThreshold = mimics.segment.create_mask()
    airSofterThreshold.name = 'Air softer threshold'
    mimics.segment.threshold(mask=airSofterThreshold, threshold_min=0,
                             threshold_max=400)  # thresholds are set in gray values
    finalMaskCheek = mimics.segment.boolean_operations(boneAndShellCropped, airSofterThreshold, operation='Minus')
    mimics.segment.keep_largest(finalMaskCheek)
    finalMaskCheek.name = 'Cheek, unite me!'

    # Deleting shit
    mimics.data.objects.delete(airSofterThreshold)

    # Deleting all the masks except for the final one
    objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
                   'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background']
    for maskName in objectsList:
        mask = mimics.data.masks.find(maskName)
        if mask is not None:
            mask.visible = False

    # Deleting all the masks except for the final one
    objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
                   'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background', 'Air softer threshold']
    for maskName in objectsList:
        mask = mimics.data.masks.find(maskName)
        if mask is not None:
            mimics.data.objects.delete(mask)

else:
    # parameters
    dilateAirPixels = 2  # Number of pixels that the air spaces are dilated to create bone shell
    closingPixels = 3  # Number of pixels used for the closing operation

    # Create basic skull using the standard threshold to select boundarypoint
    print('Starting!')
    bone = mimics.segment.create_mask()
    bone.name = 'Bone'
    mimics.segment.threshold(mask=bone, threshold_min=1250,
                             threshold_max=4095)  # thresholds are set in gray values
    mimics.segment.keep_largest(
        bone)  # Important to do this otherwise floating pixels will be merged during closing opereation!

    boundaryPointOrbital = mimics.analyze.indicate_point(message='Indicate point below the middle of the right eye-socket',
                                                         show_message_box=True,
                                                         confirm=False, title='Orbital point')
    boundaryPointCheek = mimics.analyze.indicate_point(message='Indicate point below the middle of the right maxillar sinus hole',
                                                     show_message_box=True,
                                                     confirm=False, title='Cheek point')

    print('Creating boundary box')
    bboxBigOrbital = mimics.BoundingBox3d(origin=[boundaryPointOrbital[0] - 200, boundaryPointOrbital[1] - 12, boundaryPointOrbital[2] - 5],
                                          first_vector=[500, 0, 0], second_vector=[0, 65, 0], third_vector=[0, 0, 50])
    bboxBigCheek = mimics.BoundingBox3d(origin=[boundaryPointCheek[0] - 200, boundaryPointCheek[1]-15, boundaryPointCheek[2]],
                                     first_vector=[500, 0, 0], second_vector=[0, 30, 0], third_vector=[0, 0, 35])

    print('Creating bone shells around air spaces')
    airAndBackground = mimics.segment.create_mask()
    mimics.segment.threshold(mask=airAndBackground, threshold_min=0,
                             threshold_max=350)  # thresholds are set in gray values
    airOrbital = mimics.data.objects.duplicate(airAndBackground)
    airOrbital = mimics.segment.crop_mask(airOrbital, bboxBigOrbital)
    airCheek = mimics.data.objects.duplicate(airAndBackground)
    airCheek = mimics.segment.crop_mask(airCheek, bboxBigCheek)
    air=mimics.segment.boolean_operations(airCheek, airOrbital, operation='Unite')
    air.name = 'Air'
    mimics.data.objects.delete(airOrbital)
    mimics.data.objects.delete(airCheek)
    mimics.data.objects.delete(airAndBackground)

    background1 = mimics.data.objects.duplicate(air)
    background1.name = 'Background 1'
    background1 = mimics.segment.keep_largest(background1)

    #This piece of code prevents that too much air is deleted in case background1 already captured all the
    #background air
    bbxBackground1 = mimics.measure.get_bounding_box(background1, first_axis=[1, 0, 0], second_axis=[0, 1, 0])

    if abs(bbxBackground1.first_vector[0])<170:
        background2 = mimics.segment.boolean_operations(air, background1, operation='Minus')
        background2 = mimics.segment.keep_largest(background2)
        background2.name = 'Background 2'
        background = mimics.segment.boolean_operations(background2, background1, operation='Unite')
    else:
        background = background1
    background.name = 'Background'
    airInSkull = mimics.segment.boolean_operations(air, background, operation='Minus')
    airInSkull.name = 'Air in skull'
    airInSkullDilated = mimics.segment.morphology_operations(input_mask=airInSkull, operation='Dilate',
                                                             number_of_pixels=dilateAirPixels,
                                                             connectivity=26,
                                                             limited_to_mask=None)
    airInSkullDilated.name = 'Air in skull dilated'

    boneShell = mimics.segment.boolean_operations(airInSkullDilated, airInSkull, operation='Minus')
    boneShell.name = 'Bone shell'

    temporary = mimics.segment.boolean_operations(bone, boneShell, operation='Unite')
    boneAndShell = mimics.segment.morphology_operations(input_mask=temporary, operation='Close',
                                                        number_of_pixels=closingPixels,
                                                        connectivity=26,
                                                        limited_to_mask=None)
    mimics.data.objects.delete(temporary)
    boneAndShell.name = 'Bone and shell'
    boneAndShell.visible = False

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


    # Creating the final mask by joining bone and the shell around the air
    airSofterThreshold = mimics.segment.create_mask()  # We make this mask to distract air particles from final segmentation
    airSofterThreshold.name = 'Air softer threshold'
    mimics.segment.threshold(mask=airSofterThreshold, threshold_min=0,
                             threshold_max=400)  # thresholds are set in gray values
    finalMaskOrbital = mimics.segment.boolean_operations(boneAndShellCropped, airSofterThreshold, operation='Minus')
    mimics.segment.keep_largest(finalMaskOrbital)
    finalMaskOrbital.name = 'Unite me!'


    # Deleting all the masks except for the final one
    objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
                   'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background']
    for maskName in objectsList:
        mask = mimics.data.masks.find(maskName)
        if mask is not None:
            mask.visible = False

    # Deleting all the masks except for the final one
    objectsList = ['Air', 'Background 1', 'Background 2', 'Fuchsia', 'Air in skull',
                   'Air in skull dilated', 'Bone shell', 'Bone and shell', 'Bone and shell cropped','Background', 'Air softer threshold']
    for maskName in objectsList:
        mask = mimics.data.masks.find(maskName)
        if mask is not None:
            mimics.data.objects.delete(mask)
print('Done')
