


f = open(path, 'r', encoding='utf8')
while True:
    try:
        line = f.readline()
    except UnicodeDecodeError:
        f.close()
        encodeing = 'ansi'
        break
    if not line:
        f.close()
        encoding = 'utf8'
        break

# now open your file for actual reading and data handling
with open(path, 'r', encoding=encoding) as f:

# OR

import codecs

last_position = -1

def mixed_decoder(unicode_error):
    global last_position
    string = unicode_error[1]
    position = unicode_error.start
    if position <= last_position:
        position = last_position + 1
    last_position = position
    new_char = string[position].decode("ansi")
    #new_char = u"_"
    return new_char, position + 1

codecs.register_error("mixed", mixed_decoder)
# Bonus: reads as UTF-8 until an error occurs and does not need in-place error handling.
