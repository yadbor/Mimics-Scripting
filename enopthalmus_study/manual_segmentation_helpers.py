

muscle_nerve = mimics.segment.create_mask()
muscle_nerve.name = "Manual_Muscle+Nerve"
haematoma = mimics.segment.create_mask()
haematoma.name = "Manual_Haematoma"

mimics.segment.activate_multiple_slice_edit(mask=haematoma, operation="Threshold", edit_type="ellipse")

haematoma_trimmed = mimics.segment.boolean_operations(mask_a=haematoma, mask_b=muscle_nerve, operation='Minus')
haematoma_trimmed.name = "Manual Haematoma - Manual_Muscle+Nerve" 

haematoma_in_orbit = mimics.segment.boolean_operations(mask_a=haematoma_trimmed, mask_b=orbit, operation='Intersect')
haematoma_in_orbit.name = "Haematoma in Orbital Volume"

