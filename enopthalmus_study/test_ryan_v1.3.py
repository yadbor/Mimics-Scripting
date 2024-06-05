import mimics

import os
os.chdir(r'D:\Projects & Research\Enophthalmos Study\Mimics-Scripting\enopthalmus_study')
import time

import utils, orbital_analysis
from utils import Events

def segment_orbit(rim, globe, point, side):

    t = Events()

    # Make the anterior blocking mask and the Orbital ROI
    with mimics.disabled_gui_update():
        m_anterior = orbital_analysis.make_anterior_mask(rim, globe)
        m_anterior.name = 'm_anterior'
        #mimics.data.masks[-1].name = 'm_anterior'
        # Make a bigger one as well to plug gaps
        m_ant_big = mimics.segment.morphology_operations(m_anterior, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_ant_big', limited_to_mask=None)
        # Hide all the planes as they get in the way of seeing the rim etc.
        for p in mimics.data.planes:
            p.visible = False

    t.add('make anterior mask')

    orbit_ROI = orbital_analysis.make_orbit_ROI(rim, point)
    t.add('make orbit ROI')

    # Create the bone mask, limited to the orbit_ROI
    #m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(148), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    # Make the bone lower threshold split the difference between Bone (CT) [226 : 3071] and Spongial BOne (Adult CT) [148 : 661]?
    #m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(187), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)

    # A different approach is to use the bone thresholds and then use a local threshold to add lower density material.
    m_bone = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(226), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    m_bone = mimics.segment.local_threshold(mask=m_bone, threshold_min=mimics.segment.HU2GV(148), threshold_max=mimics.segment.HU2GV(3071), search_distance=2, isolate=True)
    mimics.data.masks[-1].name = 'm_bone'

    # Could use keep largest (on uncropped mask), but that runs the risk of losing disconnected bone fragments.
    # Instead, clean up by removing small areas with 'Open' 1 x 8
    #m_bone = mimics.segment.morphology_operations(m_bone, operation='Open', number_of_pixels=1, connectivity=8, target_mask_name='m_bone', limited_to_mask=None)

    # Better appears to be to use smooth_mask
    m_bone = mimics.segment.smooth_mask(m_bone)
    mimics.data.masks[-1].name = 'm_bone'

    p_bone = mimics.segment.calculate_part(mask=m_bone, quality='High')
    mimics.data.parts[-1].name = 'p_bone'

    # Create the air mask, limited to the orbit_ROI
    # Use new 'Air' threshold, based on mimics skin = [-718 : -177] (air was [-1024 : -200])
    m_air = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(-200), bounding_box=orbit_ROI)
    mimics.data.masks[-1].name = 'm_air'

    t.add('thresholding')

    # Use erode to seperate the parts
    m_air_erode = mimics.segment.morphology_operations(m_air, operation='Erode', number_of_pixels=2, connectivity=26, target_mask_name='m_air_erode', limited_to_mask=None)

    # Want the air outside of the face. Should be the largest, and could make bbox longer anterior to make that more likely.
    # m_air_external = mimics.segment.keep_largest(mimics.data.masks[-1])

    # Or, could region_grow using the front of the bbox to set a seed point.
    # The most antero-lateral point on left eye is away from origin, on the right it is the origin. 
    # The point is offset a small amount (defaults delta = 2 mm) into the box to make sure it is inside the mask.
    lateral_pt = utils.antero_lateral(orbit_ROI, side)
    m_air_external = mimics.segment.region_grow(input_mask=m_air_erode, target_mask=None, point=lateral_pt,
                                                slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    # Dilate bask to previous size, but should only have the external part of the air now
    m_air_external = mimics.segment.morphology_operations(m_air_external, operation='Dilate', number_of_pixels=10, connectivity=26, target_mask_name='m_air_external', limited_to_mask=m_air)

    # Subtract the external from all the air in the ROI, to leave the internal air
    m_air_internal = mimics.segment.boolean_operations(m_air, m_air_external, 'Minus')
    mimics.data.masks[-1].name = 'm_air_internal'

    t.add('split air mask')

    # Open internal air to get rid of small regions caused by low density tissue overlapping with our definition of air.
    m_air_int_open = mimics.segment.morphology_operations(m_air_internal, operation='Open', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_open', limited_to_mask=None)

    # Dilate internal air 2 x 8 and unite with bone (maybe 1 x 8 would be better?)
    m_air_int_dilate = mimics.segment.morphology_operations(m_air_int_open, operation='Dilate', number_of_pixels=2, connectivity=8, target_mask_name='m_air_int_dilate', limited_to_mask=None)
    m_unite_bone_air = mimics.segment.boolean_operations(m_air_int_dilate, m_bone, 'Unite')
    mimics.data.masks[-1].name = 'm_bone+air'

    # Clean up the 'Opened' air mask
    mimics.data.masks.delete(m_air_int_open)

    t.add('dilate and unite')

    # Smartfill the 'orbital boundary' mask to close small holes
    m_fill = mimics.segment.smart_fill_global(mask=m_unite_bone_air, hole_closing_distance=7) # Is this distance too big?
    mimics.data.masks[-1].name = 'm_fill'

    t.add('smart fill')

    # Create a part from the filled mask, smooth it and wrap it to close more gaps in orbital boundary, then convert back to a mask
    m_part = mimics.segment.calculate_part(mask=m_fill, quality='High')
    mimics.data.parts[-1].name = 'm_part'

    t.add('calculate part')

    # Note keep_originals=True for debugging - not used later so could be deleted here
    m_smooth = mimics.tools.smooth(m_part, smooth_factor=0.5, iterations=5, compensate_shrinkage=False, keep_originals=True)
    mimics.data.parts[-1].name = 'm_smooth'

    t.add('smooth')

    p_wrap = mimics.tools.wrap(m_smooth, smallest_detail=0.2, gap_closing_distance=10,
                                dilate_result=False, protect_thin_walls=True, keep_originals=True)
    mimics.data.parts[-1].name = 'p_wrap'

    t.add('wrap')

    # Make the smoothed, wrapped part back into a mask. It should now define the orbit boundary
    m_wrapped = mimics.segment.calculate_mask_from_part(part=p_wrap, target_mask=None)
    mimics.data.masks[-1].name = 'm_wrapped'

    t.add('mask from part')

    # Remove tissue in front of the orbit to isolate the contents. This should be just subtracting the 'anterior' mask,
    # but due to the resolution of the blocks there can still be gaps between the bone and the anterior mask.

    # Make two masks - one expanded to segment the contents and one to size to expand back into
    # Add the anterior block to the subtract mask, to remove tissue in front of the orbital rim
    m_not_orbit = mimics.segment.boolean_operations(m_wrapped, m_anterior, 'Unite')
    mimics.data.masks[-1].name = 'm_subtract'
    # Now make the bigger one
    m_oversize = mimics.segment.morphology_operations(m_not_orbit, operation='Close', number_of_pixels=2, connectivity=26, target_mask_name='m_oversize', limited_to_mask=None)

    # Create a block of everything (should be soft tissue + air?)
    m_temp_orbit = mimics.segment.threshold(mimics.segment.create_mask(), mimics.segment.HU2GV(-1024), mimics.segment.HU2GV(3071), bounding_box=orbit_ROI)
    mimics.data.masks[-1].name = 'm_temp_orbit'

    t.add('make temp orbit')

    # Subtract the to-size orbital boundary from above, leaving the orbital contents + extra tissue areas
    m_vol_to_fill = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_not_orbit, operation='Minus')
    mimics.data.masks[-1].name = 'm_vol_to_fill'

    # Subtract the oversize orbital boundary from above, hopefully seperating the orbit from anything else
    m_segmented = mimics.segment.boolean_operations(mask_a=m_temp_orbit, mask_b=m_oversize, operation='Minus')

    # then shrink what's left to get seperate regions
    m_segmented = mimics.segment.morphology_operations(m_segmented, operation='Open', number_of_pixels=6, connectivity=8, target_mask_name='m_segmented', limited_to_mask=None)
    #mimics.data.masks[-1].name = 'm_segmented'

    t.add('m_subtract minus temp_orbit')

    seed = (globe.center[0], globe.center[1] + globe.radius, globe.center[2])

    # Region grow the soft tissue - boundary from the back of the globe, to seperate intra-orbital component
    m_core = mimics.segment.region_grow(input_mask=m_segmented, target_mask=None, point=seed,
                                                slice_type='Axial', keep_original_mask=True, multiple_layer=True, connectivity='6-connectivity')
    mimics.data.masks[-1].name = 'm_core'

    t.add('region grow core')

    # Dilate the core region back out to fill the volume, limited to the original volume to fill mask
    m_orbit = mimics.segment.morphology_operations(m_core, operation='Dilate', number_of_pixels=10, connectivity=8, target_mask_name='m_orbit', limited_to_mask=m_vol_to_fill)

    # This will extrude anteriorly if there are gaps between the bone & anterior mask, so subtract the elarged anterior mask
    m_orbit = mimics.segment.boolean_operations(mask_a=m_orbit, mask_b=m_ant_big, operation='Minus')
    # Subtract the wrapped bone again for good measure
    m_orbit = mimics.segment.boolean_operations(mask_a=m_orbit, mask_b=m_wrapped, operation='Minus')

    # There may be bits outside fo the orbit (not sure how), so region grow to get onlu connected parts, then smooth the mask
    m_orbit = mimics.segment.region_grow(input_mask=m_orbit, target_mask=None, point=seed,
                                                slice_type='Axial', keep_original_mask=False, multiple_layer=True, connectivity='6-connectivity')
    m_orbit = mimics.segment.smooth_mask(m_orbit)


    print(f'Total time taken {t.elapsed()} seconds')
    print(t.as_dict())

    return m_orbit
