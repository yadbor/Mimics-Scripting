import mimics
t = [time.process_time()]; t_label =['init']
m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(148), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_bone'
p_bone = mimics.segment.calculate_part(mask=m_bone, quality='High')

# m_temp_air_c = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
# mimics.data.masks[-1].name = 'm_temp_air_c'

# m_temp_air_ck = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
# mimics.segment.keep_largest(mimics.data.masks[-1])
# mimics.data.masks[-1].name = 'm_temp_air_ck'

# Use new 'Air' threshold, based on mimics skin = [-718 : -177] (was [-1024 : -200])
m_air_all = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-700))
mimics.data.masks[-1].name = 'm_air_all'
m_air_cropped = mimics.segment.crop_mask(mimics.data.masks.duplicate(m_air_all), orbit_ROI)
mimics.data.masks[-1].name = 'm_air_cropped'

t.append(time.process_time()); t_label.append('thresholding')

# m_air_external = mimics.segment.keep_largest(mimics.data.masks.duplicate(m_air_all))
# mimics.data.masks[-1].name = 'm_air_external'
# m_air_internal = mimics.segment.boolean_operations(m_air_all, m_air_external, 'Minus')
# mimics.data.masks[-1].name = 'm_air_intenal'

# Use erode to seperate the parts
m_air_crop_erode = mimics.segment.morphology_operations(m_air_cropped, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name='m_air_cropped_largest', limited_to_mask=None)
# Want the part outside of the face. Should be the largest, and could make bbox longer anterior to make 
# that more likely.
m_air_external = mimics.segment.keep_largest(mimics.data.masks[-1])

# m_air_external = mimics.segment.morphology_operations(m_air_external, operation='Dilate', number_of_pixels=2, connectivity=26, target_mask_name='m_air_external', limited_to_mask=None)

# Or, could region grow using the front of the bbox to set a seed point.
# The moost antero-lateral point on left eye is away from origin, on the right it is the origin
lateral_pt = utils.antero_lateral(orbit_ROI, side)
m_air_external = mimics.segment.region_grow(input_mask=m_air_crop_erode, target_mask=None, point=lateral_pt,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')


m_air_internal = mimics.segment.boolean_operations(m_air_cropped, m_air_external, 'Minus')
mimics.data.masks[-1].name = 'm_air_internal'

t.append(time.process_time()); t_label.append('split air mask')

# m_air_ext_cropped = mimics.segment.keep_largest(mimics.data.masks.duplicate(m_air_cropped))
# mimics.data.masks[-1].name = 'm_air_cropped_largest'
# m_air_int_cropped = mimics.segment.boolean_operations(m_air_cropped, m_air_ext_cropped, 'Minus')
# mimics.data.masks[-1].name = 'm_air_int_cropped'

# m_temp_air_kc = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
# mimics.segment.keep_largest(mimics.data.masks[-1])
# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
# mimics.data.masks[-1].name = 'm_temp_air_kc'

# t.append(time.process_time()); t_label.append('thresholding')

# m_morph_ck = mimics.segment.morphology_operations(m_temp_air_ck, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=None, limited_to_mask=None)
# mimics.data.masks[-1].name = 'm_morph_ck'
# m_morph_kc = mimics.segment.morphology_operations(m_temp_air_kc, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=None, limited_to_mask=None)
# mimics.data.masks[-1].name = 'm_morph_kc'

# t.append(time.process_time()); t_label.append('dilate')

# Dilate internal air 2x8 and unite with bone
m_air_int_dilate = mimics.segment.morphology_operations(m_air_internal, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_dilate', limited_to_mask=None)

m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
mimics.data.masks[-1].name = 'm_unite_bone_air'

# m_bool_ck = mimics.segment.boolean_operations(m_morph_ck, m_bone, 'Unite')
# mimics.data.masks[-1].name = 'm_bool_ck'
# m_bool_kc = mimics.segment.boolean_operations(m_morph_kc, m_bone, 'Unite')
# mimics.data.masks[-1].name = 'm_bool_kc'

t.append(time.process_time()); t_label.append('boolean')

# m_fill_ck = mimics.segment.smart_fill_global(mask=m_bool_ck, hole_closing_distance=7)
# mimics.data.masks[-1].name = 'm_fill_ck'
# m_fill_kc = mimics.segment.smart_fill_global(mask=m_bool_kc, hole_closing_distance=7)
# mimics.data.masks[-1].name = 'm_fill_kc'

m_fill = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=7)
mimics.data.masks[-1].name = 'm_fill'

t.append(time.process_time()); t_label.append('smart fill')

# m_part_ck = mimics.segment.calculate_part(mask=m_fill_ck, quality='High')
# mimics.data.parts[-1].name = 'm_part_ck'
# m_part_kc = mimics.segment.calculate_part(mask=m_fill_kc, quality='High')
# mimics.data.parts[-1].name = 'm_part_kc'

m_part = mimics.segment.calculate_part(mask=m_fill, quality='High')
mimics.data.parts[-1].name = 'm_part'

t.append(time.process_time()); t_label.append('calculate part')

# m_smooth_ck = mimics.tools.smooth(m_part_ck, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
# mimics.data.parts[-1].name = 'm_smooth_ck'
# 
m_smooth = mimics.tools.smooth(m_part, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
mimics.data.parts[-1].name = 'm_smooth'

t.append(time.process_time()); t_label.append('smooth')

# m_wrap_ck = mimics.tools.wrap(m_smooth_ck, smallest_detail=0.2, gap_closing_distance=10,
#                               dilate_result=False, protect_thin_walls=True, keep_originals=True)
# mimics.data.parts[-1].name = 'm_wrap_ck'
# m_wrap_kc = mimics.tools.wrap(m_smooth_kc, smallest_detail=0.2, gap_closing_distance=10,
#                               dilate_result=False, protect_thin_walls=True, keep_originals=True)
# mimics.data.parts[-1].name = 'm_wrap_kc'

m_wrap = mimics.tools.wrap(m_smooth, smallest_detail=0.2, gap_closing_distance=10,
                              dilate_result=False, protect_thin_walls=True, keep_originals=True)
mimics.data.parts[-1].name = 'm_wrap'

t.append(time.process_time()); t_label.append('wrap')

# m_subtract_ck = mimics.segment.calculate_mask_from_part(part=m_wrap_ck, target_mask=None)
# mimics.data.masks[-1].name = 'm_subtract_ck'
# m_subtract_kc = mimics.segment.calculate_mask_from_part(part=m_wrap_kc, target_mask=None)
# mimics.data.masks[-1].name = 'm_subtract_kc'

m_subtract = mimics.segment.calculate_mask_from_part(part=m_wrap, target_mask=None)
mimics.data.masks[-1].name = 'm_subtract'

t.append(time.process_time()); t_label.append('mask from part')

m_temp_orbit = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_temp_orbit'

t.append(time.process_time()); t_label.append('make temp orbit')

# m_intersect_ck = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract_ck, operation='Minus')
# mimics.data.masks[-1].name = 'm_intersect_ck'
# m_intersect_kc = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract_kc, operation='Minus')
# mimics.data.masks[-1].name = 'm_intersect_kc'

m_intersect = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract, operation='Minus')
mimics.data.masks[-1].name = 'm_intersect'

t.append(time.process_time()); t_label.append('m_subtract minus temp_orbit')

seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])

# m_orbit_vol_ck = mimics.segment.region_grow(input_mask=m_intersect_ck, target_mask=None, point=seed,
#                                             slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
# mimics.data.masks[-1].name = 'm_orbit_vol_ck'
# m_orbit_vol_kc = mimics.segment.region_grow(input_mask=m_intersect_kc, target_mask=None, point=seed,
#                                             slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
# mimics.data.masks[-1].name = 'm_orbit_vol_kc'

m_orbit_vol = mimics.segment.region_grow(input_mask=m_intersect, target_mask=None, point=seed,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
mimics.data.masks[-1].name = 'm_orbit_vol'

t.append(time.process_time()); t_label.append('region grow')

steps = dict(zip(t_label[1:], [b-a for a,b in utils.pairwise(t)]))
print(steps)

# new_subtracted = mimics.segment.boolean_operations(m_subtract_kc, m_fill_kc, operation='Minus')
# new_morph = mimics.segment.morphology_operations(new_subtracted,operation='Erode',number_of_pixels=3,connectivity=8)
# new_grow=mimics.segment.region_grow(new_morph, None, point=seed, slice_type='Axial', keep_original_mask=True, multiple_layer=True,connectivity='6-connectivity')
# new_dilate = mimics.segment.morphology_operations(new_grow,operation='Dilate',number_of_pixels=3,connectivity=8)
# new_orbit = mimics.segment.boolean_operations(new_dilate, m_bone, operation='Minus')
