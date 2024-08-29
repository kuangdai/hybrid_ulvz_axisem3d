def replace_in_file(source, replace_dict, dest=None):
    if dest is None:
        dest = source
    with open(source, 'r') as fs:
        text = fs.read()
    for key, val in replace_dict.items():
        text = text.replace(key, str(val))
    with open(dest, 'w') as fs:
        fs.write(text)
