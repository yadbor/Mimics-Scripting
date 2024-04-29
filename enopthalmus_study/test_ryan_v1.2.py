import mimics
import orbital_analysis

t = [time.perf_counter()]; t_label =['init']

# Make the anterior blocking mask and the Orbital ROI
with mimics.disabled_gui_update():
    m_anterior = orbital_analysis.make_anterior_mask(rim, globe)
    m_anterior.name = 'm_anterior'
    #mimics.data.masks[-1].name = 'm_anterior'

# Hide all the planes as they get in the way of seeing the rim etc.
for p in mimics.data.planes:
    p.visible = False

t.append(time.perf_counter()); t_label.append('make anterior mask')

orbit_ROI = orbital_analysis.make_orbit_ROI(rim, point)
print(f'orbit ROI is {orbit_ROI}')

t.append(time.perf_counter()); t_label.append('make orbit ROI')

# Make the bone lower threshold split the difference between Bone (CT) [226 : 3071] and Spongial BOne (Adult CT) [148 : 661]
m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(187), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)

# Could use keep largest (on uncropped mask), but that runs the risk of losing disconnected bone fragments.
# Instead, clean up by removing small areas with 'Open' 1 x 8
m_bone = mimics.segment.morphology_operations(m_bone, operation='Open', number_of_pixels=1, connectivity=8, target_mask_name='m_bone', limited_to_mask=None)

# mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_bone'
p_bone = mimics.segment.calculate_part(mask=m_bone, quality='High')

# Use new 'Air' threshold, based on mimics skin = [-718 : -177] (air was [-1024 : -200])
m_air_all = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200))
mimics.data.masks[-1].name = 'm_air_all'
m_air_cropped = mimics.segment.crop_mask(mimics.data.masks.duplicate(m_air_all), orbit_ROI)
mimics.data.masks[-1].name = 'm_air_cropped'

t.append(time.perf_counter()); t_label.append('thresholding')

# Use erode to seperate the parts
m_air_crop_erode = mimics.segment.morphology_operations(m_air_cropped, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name='m_air_crop_erode', limited_to_mask=None)

# Want the part outside of the face. Should be the largest, and could make bbox longer anterior to make 
# that more likely.
# m_air_external = mimics.segment.keep_largest(mimics.data.masks[-1])

# Or, could region grow using the front of the bbox to set a seed point.
# The most antero-lateral point on left eye is away from origin, on the right it is the origin
lateral_pt = utils.antero_lateral(orbit_ROI, side)
m_air_external = mimics.segment.region_grow(input_mask=m_air_crop_erode, target_mask=None, point=lateral_pt,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')

# Dilate bask to previous size, but should only have the external part of the air now
m_air_external = mimics.segment.morphology_operations(m_air_external, operation='Dilate', number_of_pixels=10, connectivity=26, target_mask_name='m_air_external', limited_to_mask=m_air_cropped)
# Could dialate more here, to catch *everything* outside, then subtract a soft tissue mask to clean up the skin surface before subtracting from total air.
# Cover both skin [-718 : -177] and soft tissue [-700 : 225]
m_soft_tissue = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-720), mimics.segment.HU2GV(225), bounding_box=orbit_ROI)
mimics.data.masks[-1].name = 'm_soft_tissue'
m_air_external = mimics.segment.boolean_operations(m_air_external, m_soft_tissue, 'Minus')

# Subtract the external from all the air in the ROI, to leave the internal air
m_air_internal = mimics.segment.boolean_operations(m_air_cropped, m_air_external, 'Minus')
mimics.data.masks[-1].name = 'm_air_internal'

t.append(time.perf_counter()); t_label.append('split air mask')

# Open internal air to get rid of small regions (caused by low deisnity tissue overlapping with our definition of air)
m_air_int_open = mimics.segment.morphology_operations(m_air_internal, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_open', limited_to_mask=None)


# Dilate internal air 2 x 8 and unite with bone (maybe 1 x 8 would be better?)
m_air_int_dilate = mimics.segment.morphology_operations(m_air_int_open, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_dilate', limited_to_mask=None)
m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
mimics.data.masks[-1].name = 'm_unite_bone_air'

t.append(time.perf_counter()); t_label.append('dilate and unite')

m_fill = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=7) # Is this distance too big?
mimics.data.masks[-1].name = 'm_fill'

t.append(time.perf_counter()); t_label.append('smart fill')

m_part = mimics.segment.calculate_part(mask=m_fill, quality='High')
mimics.data.parts[-1].name = 'm_part'

t.append(time.perf_counter()); t_label.append('calculate part')

m_smooth = mimics.tools.smooth(m_part, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
mimics.data.parts[-1].name = 'm_smooth'

t.append(time.perf_counter()); t_label.append('smooth')

m_wrap = mimics.tools.wrap(m_smooth, smallest_detail=0.2, gap_closing_distance=10,
                              dilate_result=False, protect_thin_walls=True, keep_originals=True)
mimics.data.parts[-1].name = 'm_wrap'

t.append(time.perf_counter()); t_label.append('wrap')

# Make the smoothed, wrapped part back into a mask. It should now define the orbit boundary
m_subtract = mimics.segment.calculate_mask_from_part(part=m_wrap, target_mask=None)
mimics.data.masks[-1].name = 'm_subtract'

t.append(time.perf_counter()); t_label.append('mask from part')

# Add the anterior block to the subtract mask, to remove tissue in front of the eorbital rim
m_subtract = mimics.segment.boolean_operations(m_subtract, m_anterior, 'Unite')

# Create a block of everything (shoudl be soft tissue + air?)
m_temp_orbit = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071))
mimics.segment.crop_mask(mimics.data.masks[-1], orbit_ROI)
mimics.data.masks[-1].name = 'm_temp_orbit'

t.append(time.perf_counter()); t_label.append('make temp orbit')
# Subtracat the orbital bounday from above
m_intersect = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_subtract, operation='Minus')
mimics.data.masks[-1].name = 'm_intersect'

t.append(time.perf_counter()); t_label.append('m_subtract minus temp_orbit')

seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])

# Region grow trhe soft rissue - boundary from the back of the globe, to seperate intra-orbital component
m_orbit_vol = mimics.segment.region_grow(input_mask=m_intersect, target_mask=None, point=seed,
                                            slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
mimics.data.masks[-1].name = 'm_orbit_vol'

t.append(time.perf_counter()); t_label.append('region grow')

steps = dict(zip(t_label[1:], [b-a for a,b in utils.pairwise(t)]))
print(steps)

# new_subtracted = mimics.segment.boolean_operations(m_subtract_kc, m_fill_kc, operation='Minus')
# new_morph = mimics.segment.morphology_operations(new_subtracted,operation='Erode',number_of_pixels=3,connectivity=8)
# new_grow=mimics.segment.region_grow(new_morph, None, point=seed, slice_type='Axial', keep_original_mask=True, multiple_layer=True,connectivity='6-connectivity')
# new_dilate = mimics.segment.morphology_operations(new_grow,operation='Dilate',number_of_pixels=3,connectivity=8)
# new_orbit = mimics.segment.boolean_operations(new_dilate, m_bone, operation='Minus')
