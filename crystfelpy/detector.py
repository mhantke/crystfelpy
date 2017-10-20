#
# Reading Crytfel geometry files
#
# The code below is based on bits of the CrystFEL C library (Thomas White <taw@physics.org> et al.)
# Author: Max Hantke (hantke@xray.bmc.uu.se)
#

import re, copy
import numpy as np

count = lambda a_list, value: len([element for element in a_list if (element==value)])

def get_detector_geometry(filename):

    reject = False

    rg_tmp = {}
    rgc_tmp = {}
    
    det = {
	"panels"   : {},
	"bad"      : {},
	"mask_good": 0,
	"mask_bad" : 0,
	"rigid_groups" : {},
	"path_dim" : 0,
	"dim_dim" : 0,
	"rigid_group_collections" : {},
	# The default defaults... 
	"defaults" : { "orig_min_fs" : None,
	               "orig_min_ss" : None,
	               "orig_max_fs" : None,
	               "orig_max_ss" : None,
	               "cnx" : None,
	               "cny" : None,
	               "clen" : None,
                       "clen_from" : None,
	               "coffset" : 0.0,
	               "res" : None,
	               "badrow" : '-',
	               "no_index" : 0,
	               "fsx" : 1.0,
	               "fsy" : 0.0,
	               "fsz" : 0.0,
	               "ssx" : 0.0,
	               "ssy" : 1.0,
	               "ssz" : 0.0,
                       "rail_x" : None,  # The actual default rail direction 
                       "rail_y" : None,  # is below
                       "rail_z" : None,
                       "clen_for_centering" : None,
	               "adu_per_eV" : None,
	               "adu_per_photon" : None,
	               "max_adu" : None,
	               "mask" : "",
	               "mask_file" : None,
	               "satmap" : None,
	               "satmap_file" : None,
	               "data" : "",
	               "dim_structure" : [],
	               #"name" : "", # Maybe not needed
        }
    }

    bad_default = {
	"min_x" : np.nan,
	"max_x" : np.nan,
	"min_y" : np.nan,
	"max_y" : np.nan,
	"min_fs": 0,
	"max_fs": 0,
	"min_ss": 0,
	"max_ss": 0,
	"is_fsss": None,
	"panel" : None,
    }

    default_dim_structure = ['ss', 'fs']
        
    with open(filename, 'r') as f:
        for line in f.readlines():
            line_wo_comment = line.split(';')[0].strip()
            if len(line_wo_comment) == 0:
                continue
            if line_wo_comment.count('=') != 1:
                print "ERROR: Syntax error. Every line must contain one equal sign separating key and value pair. (%s)" % line_wo_comment
            name, value = line_wo_comment.split('=')
            name = name.strip()
            value = value.strip()
            name_bits = re.split('[/\\.]+', name)
            n0 = name_bits[0]
            if len(name_bits) < 2:
                # Top level option
                if n0 in ['mask_good', 'mask_bad']:
                    det[n0] = int(value)
                elif n0.startswith('rigid_group_collection'):
                    rgc_name = n0[len('rigid_group_collection'):].strip('_')
                    rgc_tmp[rgc_name] = value.split(',')
                elif n0.startswith('rigid_group'):
                    rg_name = n0[len('rigid_group'):].strip('_')
                    rg_tmp[rg_name] = value.split(',')
                elif n0 in ['photon_energy', 'photon_energy_from', 'photon_energy_scale']:
                    print "WARNING: Ignoring the entry \'%s\'." % n0
                else:
                    unrec = _parse_arg_for_panel(det, det['defaults'], None, n0, value)
                    if unrec:
                        print "ERROR: Unrecongnised top level field \'%s\'" % n0
            else:
                n1 = name_bits[1]
                if n0.startswith('bad'):
                    if n0 not in det['bad']:
                        det['bad'][n0] = copy.deepcopy(bad_default)
                    badregion = det['bad'][n0]
                    unrec = _parse_arg_for_bad(det, badregion, n1, value)
                else:
                    if n0 not in det['panels']:
                        det['panels'][n0] = copy.deepcopy(det['defaults'])
		    panel = det['panels'][n0]
                    unrec = _parse_arg_for_panel(det, panel, n0, n1, value)

    if len(det['panels']) == 0:
        print "ERROR: No panel descriptions in geometry file."
        return

    # Determine and check path dimensions
        
    path_dim_data = np.array([p['data'].count('%') for p in det['panels'].values()])
    if not (path_dim_data[0] == path_dim_data).all():
        print "ERROR: All panels' data entries must have the same number of placeholders"
        reject = True

    path_dim_mask = np.array([p['mask'].count('%') for p in det['panels'].values()])
    if not (path_dim_mask[0] == path_dim_mask).all():
        print "ERROR: All panels' mask entries must have the same number of placeholders"
        reject = True
            
    if path_dim_mask[0] > path_dim_data[0]:
        print "ERROR: Number of placeholders in mask cannot be larger than for data"
        reject = True
            
    det["path_dim"] = path_dim_data[0]
            
    # Determine and check dim dimensions
        
    dim_dim = np.array([count(p['dim_structure'], '%') for p in det['panels'].values()])
    if not (dim_dim[0] == dim_dim).all():
        print "ERROR: All panels' dim entries must have the same number of placeholders"
        reject = True
    if (dim_dim > 1).any():
        print "ERROR: Maximum one placeholder dim coordinate is allowed."
            
    ss_defined = np.array([(count(p['dim_structure'], 'ss')==1) for p in det['panels'].values()]).all()
    if not ss_defined:
        print "ERROR: One slow scan dim must be defined."
        reject = True
    fs_defined = np.array([(count(p['dim_structure'], 'fs')==1) for p in det['panels'].values()]).all()

    if not fs_defined:
        print "ERROR: One fast scan dim must be defined."
        reject = True
            
    det['dim_dim'] = dim_dim

    for p_key, p in det['panels'].items():

        for s in ['min_fs', 'max_fs', 'min_ss', 'max_ss']:
            if p['orig_%s' % s] is None:
		print "ERROR: Please specify the %s coordinate for panel %s." % (s, p_key)
		reject = True

        for c in ['x', 'y']:
            if p['cn%s' % c] is None:
                print "ERROR: Please specify the corner %s coordinate for panel %s" % (c, p_key)
		reject = True
            
        if p['clen'] is None and p['clen_from'] is None:
	    print "ERROR: Please specify the camera length for panel %s" % (p_key)
	    reject = True

        if p['res'] is None:
	    print "ERROR: Please specify the resolution for panel %s" % (p_key)
	    reject = True
            
        if p['adu_per_eV'] is None and p['adu_per_photon'] is None:
            print "ERROR: Please specify either adu_per_eV or adu_per_photon for panel %s" % (p_key)
	    reject = True

        if p['clen_for_centering'] is None and not p['rail_x'] is None:
	    print "ERROR: You must specify clen_for_centering if you specify the rail direction (panel %s)" % (p_key)
	    reject = True

	if not p['mask_file'] is None and p['mask'] is None:
	    print "ERROR: You have specified \'mask_file\' but not \'mask\'. \'mask_file\' will therefore have no effect (panel %s)." % (p_key)
            reject = True

        # The default rail direction
	if p['rail_x'] is None:
	    p['rail_x'] = 0.0
	    p['rail_y'] = 0.0
	    p['rail_z'] = 1.0
            
	if p['clen_for_centering'] is None:
            p['clen_for_centering'] = 0.0

	# It's OK if the badrow direction is '0'
	# It's not a problem if "no_index" is still zero 
	# The default transformation matrix is at least valid 

        # Calculate width and height
	p['w'] = int(p['orig_max_fs'] - p['orig_min_fs'] + 1)
	p['h'] = int(p['orig_max_ss'] - p['orig_min_ss'] + 1)

    for b_key, b in det['bad'].items():
        if b['is_fsss'] is None:
            print "ERROR: Please specify the coordinate ranges for bad region %s" % (b_key)
	    reject = True

    for rg_key, panel_keys in rg_tmp.items():
        for panel_key in panel_keys:
            if panel_key not in det['panels']:
                print "ERROR: Cannot add panel to rigid group \'%s\'. Panel not found: \'%s\'" % (rg_key, panel_keys)
                reject = True
            else:
                if rg_key not in det['rigid_groups']:
                    det['rigid_groups'][rg_key] = {}
                det['rigid_groups'][rg_key][panel_key] = det['panels'][panel_key]
                
    for rgc_key, rg_keys in rgc_tmp.items():
        for rg_key in rg_keys:
            if rg_key not in det['rigid_groups']:
                print "ERROR: Cannot add rigid group to collection \'%s\'. Rigid group not found: %s" % (rgc_key, rg_key)
            else:
                if rgc_key not in det['rigid_group_collections']:
                    det['rigid_group_collections'][rgc_key] = {}
                det['rigid_group_collections'][rgc_key][rg_key] = det['rigid_groups'][rg_key]
        
    # Calculate matrix inverses and other stuff
    for p_key, p in det['panels'].items():
	if p['fsx']*p['ssy'] == p['ssx']*p['fsy']:
	    print "ERROR: Panel \'%s\' transformation singular." % (p_key)
            
        d = float(p['fsx']*p['ssy'] - p['ssx']*p['fsy'])
	p['xfs'] =  p['ssy'] / d
	p['yfs'] = -p['ssx'] / d
	p['xss'] = -p['fsy'] / d
	p['yss'] =  p['fsx'] / d
        
    min_d, max_d = find_min_max_d(det);

    if not reject:
        return det
        
# Test if fs,ss in panel "p" is further {out,in} than {*p_max_d,*p_min_d}, and
# if so update det->furthest_{out,in}_{panel,fs,ss}.
def check_point(panel, fs, ss, min_d, max_d, det):
    
    xs = fs*panel["fsx"] + ss*panel["ssx"]
    ys = fs*panel["fsy"] + ss*panel["ssy"]
    
    rx = (xs + panel["cnx"]) / panel["res"]
    ry = (ys + panel["cny"]) / panel["res"]
    
    d = np.sqrt(rx**2 + ry**2)
    
    if d > max_d:        
	det["furthest_out_panel"] = panel
	det["furthest_out_fs"] = fs
	det["furthest_out_ss"] = ss
        return min_d, d
    elif d < min_d:
	det["furthest_in_panel"] = panel
	det["furthest_in_fs"]    = fs
	det["furthest_in_ss"]    = ss
        return d, max_d
    else:
        return min_d, max_d
        
def find_min_max_d(det):
    min_d = np.inf
    max_d = 0.0
    for p in det['panels'].values():
        min_d, max_d = check_point(panel=p, fs=0,      ss=0,      min_d=min_d, max_d=max_d, det=det)
        min_d, max_d = check_point(panel=p, fs=0,      ss=p['h'], min_d=min_d, max_d=max_d, det=det)
        min_d, max_d = check_point(panel=p, fs=p['w'], ss=0,      min_d=min_d, max_d=max_d, det=det)
        min_d, max_d = check_point(panel=p, fs=p['w'], ss=p['h'], min_d=min_d, max_d=max_d, det=det)
    return min_d, max_d

def dir_conv(value):
    x = 0.
    y = 0.
    z = 0.
    for value_bit in value.split(' '):
        a = value_bit.strip()
        if a.endswith('x'):
            x = float(a[:-1])
        elif a.endswith('y'):
            y = float(a[:-1])
        elif a.endswith('z'):
            z = float(a[:-1])
        else:
            print "ERROR: Invalid direction %s" % value
            return None, None, None
    return x, y, z


def _parse_arg_for_bad(det, badregion, key, value):
    unrec = False
    xy_args = ['min_x', 'max_x', 'min_y', 'max_y']
    fsss_args = ['min_fs', 'max_fs', 'min_ss', 'max_ss']
    if key in xy_args:
        if badregion['is_fsss'] is None:
            badregion['is_fsss'] = False
        elif badregion['is_fsss'] == True:
            print "ERROR: You can't mix x/y and fs/ss in a bad region."
        else:
            badregion[key] = float(value)
    elif key in fsss_args:
        if badregion['is_fsss'] is None:
            badregion['is_fsss'] = True
        elif badregion['is_fsss'] == False:
            print "ERROR: You can't mix x/y and fs/ss in a bad region."
        else:
            badregion[key] = float(value)
    elif key == 'panel':
        badregion[key] = value
    else:
        print "ERROR: Unrecognised field \'%s\'." % (key)
        unrec = True

    return unrec
    
def _parse_arg_for_panel(det, panel, panel_key, key, value):
    unrec = False
    if key in ['min_x', 'max_x', 'min_y', 'max_y', 'min_fs', 'max_fs', 'min_ss', 'max_ss']: 
        panel['orig_%s' % key] = float(value)
    elif key in ['corner_x', 'corner_y']:
        panel['cn%s' % key[-1]] = float(value)
    elif key == 'rail_direction':
        tmp_x, tmp_y, tmp_z = dir_conv(value)
        if tmp_x is not None and tmp_y is not None and tmp_z is not None:
            panel['rail_x'] = tmp_x
            panel['rail_y'] = tmp_y
            panel['rail_z'] = tmp_z
        else:
            print "ERROR: Illegal input for rail_direction \'%s\'" % value
            unrec = True
    elif key in ['clen_for_centering', 'adu_per_eV', 'adu_per_photon']:
        panel[key] = float(value)
    elif key == 'rigid_group':
        if panel not in det['rigid_groups']:
            if panel_key is None:
                print "ERROR: Cannot add unknown panel to rigid group."
                unrec = True
            else:
                det['rigid_groups'][panel_key] = panel
    elif key == 'clen':
        try:
            float(value)
        except ValueError:
            panel['clen'] = None
            panel['clen_from'] = value
        else:
            panel['clen'] = float(value)
            panel['clen_from'] = None
    elif key in ['data', 'mask']:
        if not value.startswith("/"):
	    print "ERROR: Invalid %s location \'%s\'" % (key, value)
        else:
	    panel[key] = value
    elif key in ['mask_file', 'saturation_map', 'saturation_map_file']:
        panel[key] = value
    elif key in ['coffset', 'res', 'max_adu']:
	panel[key] = float(value)
    elif key == 'badrow_direction':
        if value not in ['x', 'y', 'f', 's', '-']:
            print "ERROR: badrow_direction must be \'x\', \'y\', \'f\', \'s\' or \'-\'. Assuming \'-\'."
	    panel['badrow'] = '-'
        else:
            panel['badrow'] = 'f' if value in ['x', 'f'] else 's' 
    elif key == "no_index":
	panel['no_index'] = int(bool(value))
    elif key == 'fs':
        tmp_x, tmp_y, tmp_z = dir_conv(value)
        if tmp_x is not None and tmp_y is not None and tmp_z is not None:
            panel['fsx'] = tmp_x
            panel['fsy'] = tmp_y
            panel['fsz'] = tmp_z
        else:
            print "ERROR: Illegal fast scan direction \'%s\'" % value
            unrec = True
    elif key == 'ss':
        tmp_x, tmp_y, tmp_z = dir_conv(value)
        if tmp_x is not None and tmp_y is not None and tmp_z is not None:
            panel['ssx'] = tmp_x
            panel['ssy'] = tmp_y
            panel['ssz'] = tmp_z
        else:
            print "ERROR: Illegal slow scan direction \'%s\'" % value
            unrec = True
    elif key.startswith('dim'):
        if value not in ['%', 'ss', 'fs']:
            print "ERROR: dim must be \'%\', \'ss\', or \'fs\'."
            unrec = True
        else:
            idim = int(key[3])
            while len(panel['dim_structure']) < (idim+1):
                panel['dim_structure'].append(None)
            panel['dim_structure'][idim] = value
    else:
        print "ERROR: Unrecognised field '%s'" % key
        unrec = True
    return unrec

def count_number_of_pixels(det):
    N = 0
    for p in det['panels'].values():
        N += int(p['h'] * p['w'])
    return N

def calc_pixel_coords(det):
    N = count_number_of_pixels(det)
    X = np.zeros(N, dtype='float64')
    Y = np.zeros(N, dtype='float64')
    i = 0
    for p in det['panels'].values():
        Nss = p['h']
        Nfs = p['w']
        Np = Nss * Nfs
        FS, SS = np.meshgrid(np.arange(Nfs), np.arange(Nss))
        XS = (FS.ravel() * p['fsx'] + SS.ravel() * p['ssx'])
        YS = (FS.ravel() * p['fsy'] + SS.ravel() * p['ssy'])
        X[i:i+Np] = (XS[:] + p['cnx']) / p['res']
        Y[i:i+Np] = (YS[:] + p['cny']) / p['res']
        i += Np
    return X, Y

def draw_mask_image(det, unit_is_pixel=True):
    res = np.array([p['res'] for p in det['panels'].values()])
    if not (res[0] == res).all():
        print "ERROR: You cannot set unit_is_pixel=True if detector contains panels with different res."
        return
    res = res[0]
    Xpix, Ypix = calc_pixel_coords(det)
    Xpix = res * Xpix
    Ypix = res * Ypix
    Xpix = np.asarray(np.round(Xpix - Xpix.min()), dtype='int64')
    Ypix = np.asarray(np.round(Ypix - Ypix.min()), dtype='int64')
    Nx = Xpix.max() + 1
    Ny = Ypix.max() + 1
    I = Ypix * Nx + Xpix 
    M = np.zeros(Ny*Nx, dtype='int32')
    M[I] += 1
    M = M.reshape((Ny, Nx))
    return M
