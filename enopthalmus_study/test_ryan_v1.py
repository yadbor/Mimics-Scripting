# orbit_ROI
# <mimics.BoundingBox3d((4.85060977935791, -113.2152328491211, 107.19224548339844), (58.89358043670654, 0.0, 0.0), (0.0, 73.4571762084961, 0.0), (0.0, 0.0, 65.81268310546875))>
# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)

import mimics

t = [time.process_time()]; t_label =['init']
m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_bone'
p_bone = mimics.segment.calculate_part(mask=m_bone, quality='High')

m_temp_air_c = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_temp_air_c'

m_temp_air_ck = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.segment.keep_largest(mimics.data.masks[-1])
mimics.data.masks[-1].name = 'm_temp_air_ck'

m_temp_air_kc = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
mimics.segment.keep_largest(mimics.data.masks[-1])
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_temp_air_kc'

t.append(time.process_time()); t_label.append('thresholding')

m_morph_ck = mimics.segment.morphology_operations(m_temp_air_ck, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=None, limited_to_mask=None)
mimics.data.masks[-1].name = 'm_morph_ck'
m_morph_kc = mimics.segment.morphology_operations(m_temp_air_kc, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name=None, limited_to_mask=None)
mimics.data.masks[-1].name = 'm_morph_kc'

t.append(time.process_time()); t_label.append('dilate')

m_bool_ck = mimics.segment.boolean_operations(m_morph_ck, m_bone, 'Unite')
mimics.data.masks[-1].name = 'm_bool_ck'
m_bool_kc = mimics.segment.boolean_operations(m_morph_kc, m_bone, 'Unite')
mimics.data.masks[-1].name = 'm_bool_kc'

t.append(time.process_time()); t_label.append('boolean')

m_fill_ck = mimics.segment.smart_fill_global(mask=m_bool_ck, hole_closing_distance=7)
mimics.data.masks[-1].name = 'm_fill_ck'
m_fill_kc = mimics.segment.smart_fill_global(mask=m_bool_kc, hole_closing_distance=7)
mimics.data.masks[-1].name = 'm_fill_kc'

t.append(time.process_time()); t_label.append('smart fill')

m_part_ck = mimics.segment.calculate_part(mask=m_fill_ck, quality='High')
mimics.data.parts[-1].name = 'm_part_ck'
m_part_kc = mimics.segment.calculate_part(mask=m_fill_kc, quality='High')
mimics.data.parts[-1].name = 'm_part_kc'

t.append(time.process_time()); t_label.append('calculate part')

m_smooth_ck = mimics.tools.smooth(m_part_ck, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
mimics.data.parts[-1].name = 'm_smooth_ck'
m_smooth_kc = mimics.tools.smooth(m_part_kc, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
mimics.data.parts[-1].name = 'm_smooth_kc'

t.append(time.process_time()); t_label.append('smooth')

m_wrap_ck = mimics.tools.wrap(m_smooth_ck, smallest_detail=0.2, gap_closing_distance=10,
                              dilate_result=False, protect_thin_walls=True, keep_originals=True)
mimics.data.parts[-1].name = 'm_wrap_ck'
m_wrap_kc = mimics.tools.wrap(m_smooth_kc, smallest_detail=0.2, gap_closing_distance=10,
                              dilate_result=False, protect_thin_walls=True, keep_originals=True)
mimics.data.parts[-1].name = 'm_wrap_kc'

t.append(time.process_time()); t_label.append('wrap')

m_subtract_ck = mimics.segment.calculate_mask_from_part(part=m_wrap_ck, target_mask=None)
mimics.data.masks[-1].name = 'm_subtract_ck'
m_subtract_kc = mimics.segment.calculate_mask_from_part(part=m_wrap_kc, target_mask=None)
mimics.data.masks[-1].name = 'm_subtract_kc'

t.append(time.process_time()); t_label.append('mask from part')

m_temp_orbit = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_temp_orbit'

t.append(time.process_time()); t_label.append('make temp orbit')

m_intersect_ck = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract_ck, operation='Minus')
mimics.data.masks[-1].name = 'm_intersect_ck'
m_intersect_kc = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract_kc, operation='Minus')
mimics.data.masks[-1].name = 'm_intersect_kc'

t.append(time.process_time()); t_label.append('subtract M_subtract from temp_orbit')

seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])

m_orbit_vol_ck = mimics.segment.region_grow(input_mask=m_intersect_ck, target_mask=None, point=seed,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
mimics.data.masks[-1].name = 'm_orbit_vol_ck'
m_orbit_vol_kc = mimics.segment.region_grow(input_mask=m_intersect_kc, target_mask=None, point=seed,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
mimics.data.masks[-1].name = 'm_orbit_vol_kc'

t.append(time.process_time()); t_label.append('region grow')

step_times = [b-a for a,b in utils.pairwise(t)]
step_labels = t_label[1:]
steps = dict(zip(step_labels, step_times))
print(steps)