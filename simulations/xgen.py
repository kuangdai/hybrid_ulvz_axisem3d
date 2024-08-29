import os
import fnmatch
import sys
from pathlib import Path

if __name__ == "__main__":
    run_name = sys.argv[1]

    # python script name
    py_files = fnmatch.filter(os.listdir('.'), '[0-9][0-9]_*.py')
    py_no_name = {py_file[:2]: py_file[3:-3] for py_file in py_files}
    py_no_name_sort = {}
    for i in range(100):
        no = f'{i:02d}'
        if no in py_no_name.keys():
            py_no_name_sort[no] = py_no_name[no]

    # run scripts
    os.system(f'rm -rf outputs/{run_name}')
    for no, name in py_no_name_sort.items():
        os.system(f'python {no}_{name}.py {run_name}')

    # copy scripts for post-processing
    os.system('cp {geodetic,seismogram,merge_rotate}.py %s' % f'outputs/{run_name}/')

    # replace @@ folder name
    for no, name in py_no_name_sort.items():
        os.system(f'mv outputs/{run_name}/@@_{name} '
                  f'outputs/{run_name}/{no}_{name}')

    # replace @@ inside
    files = [item for item in Path(f'outputs/{run_name}').rglob('*') if item.is_file()]
    for fname in files:
        fname = str(fname)
        if fname.endswith('.e') or fname.endswith('.png') or fname.startswith('STATIONS'):
            continue
        with open(fname, 'r') as fs:
            text = fs.read()
        for no, name in py_no_name_sort.items():
            text = text.replace(f'@@_{name}', f'{no}_{name}')
        with open(fname, 'w') as fs:
            fs.write(text)
