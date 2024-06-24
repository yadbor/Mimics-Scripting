import os  # for scandir() etc
import re  # for regexp matching

root = r'D:\Projects & Research\Enophthalmos Study'

# Each user analysed a seperate folder, listed here
p_list = {
'ryan': {'base': 'Ryan Processing'}, 
'rob':  {'base': 'Rob Processing'}
}
# Gather the file names for each user
for u in p_list.keys():
    entry = p_list[u]
    base = os.path.join(root, entry['base'])
    entry['projects'] = [f.path for f in os.scandir(base) if re.match(r'.*.mcs', f.name)]



for f in os.scandir(base):
    if re.match('.*Bone.*', f.name): 
        print(f.name)
        os.rename(src=f.path, dst=f.path.replace('Bone', 'bone'))
    if re.match('.*Soft.*', f.name):
        os.rename(src=f.path, dst=f.path.replace('Soft', 'soft'))


def to_mask():
    part_to_copy = [p for p in mimics.data.parts if p.selected]
    mimics.segment.calculate_mask_from_part(part=part_to_copy[0])
    mimics.data.masks[-1].name = part_to_copy.name


# Find masks related to Orbital Volume in current image set
img = mimics.data.images.get_active()
correct_orbits = [ov for ov in img.linked_objects if re.match('.*Orbital.*', ov.name) and isinstance(ov,  mimics.segment.Mask)]

def to_mask():
    for part_to_copy in [p for p in mimics.data.parts if p.selected]:
        mimics.segment.calculate_mask_from_part(part=part_to_copy)
        new_name = re.sub(r'(?<=Volume).*', '', part_to_copy.name)
        mimics.data.masks[-1].name = new_name

def tom():
    for part_to_copy in mimics.data.parts:
        if part_to_copy.selected:
            part_images = part_to_copy.image
            mimics.data.images.set_active(part_images)
            mimics.segment.calculate_mask_from_part(part=part_to_copy)
            new_name = re.sub(r'(?<=Volume).*', '', part_to_copy.name)
            mimics.data.masks[-1].name = new_name

def ren():
    for m in mimics.data.masks:
        if m.selected:
            m.name = m.name + '_old'
            mimics.segment.calculate_part(mask=m, quality='High')

    