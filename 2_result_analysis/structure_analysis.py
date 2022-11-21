from rdkit import Chem
from pymol import cmd, CmdException


res_dic = {'ALA':'A',
           'ARG':'R',
           'ASN':'N',
           'ASP':'D',
           'CYS':'C',
           'GLN':'Q',
           'GLU':'E',
           'GLY':'G',
           'HIS':'H',
           'ILE':'I',
           'LEU':'L',
           'LYS':'K',
           'MET':'M',
           'PHE':'F',
           'PRO':'P',
           'SER':'S',
           'THR':'T',
           'TRP':'W',
           'TYR':'Y',
           'VAL':'V'}
           

# get distance object information in PyMOL
def get_raw_distances(names='', state=1, selection='all', quiet=1):
    from chempy import cpv

    state, quiet = int(state), int(quiet)
    if state < 1:
        state = cmd.get_state()

    valid_names = cmd.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                print(' Error: no such distance object: ' + name)
                raise CmdException

    raw_objects = cmd.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',space=locals())

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError, IndexError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))
    return r


