from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def get_phase_cmap():
    ''' Global color map for cell cycle phases '''
    cmap = {
        'S': 'goldenrod',
        'G1/2': 'dodgerblue',
        'G1': 'dodgerblue',
        'G2': 'lightblue',
        'LQ': 'lightgrey'
    }
    return cmap


def get_cna_cmap():
    ''' Global color map for copy number alterations '''
    cmap = {
        'gain': 'red',  # red
        'loss': 'deepskyblue',  # dark blue
        'neutral': '#CCCCCC',  # grey
        'unaltered': '#CCCCCC'  # grey
    }
    return cmap


def get_rt_cmap(return_colors=False):
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    if return_colors:
        return ListedColormap(color_list), rt_colors
    return ListedColormap(color_list)


def get_acc_cmap(return_colors=False):
    """ Return a colormap for replication accuracy states. False positives are green, false negatives are purple, and correct calls are gray. """
    acc_colors = {0:'#CCCCCC', -1: '#532A44', 1: '#00685E'}
    color_list = []
    for i in [-1, 0, 1]:
        color_list.append(acc_colors[i])
    if return_colors:
        return ListedColormap(color_list), acc_colors
    return ListedColormap(color_list)        


def get_cn_cmap(cn_data):
    color_reference = {0:'#3182BD', 1:'#9ECAE1', 2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33', 6:'#B30000', 7:'#980043', 8:'#DD1C77', 9:'#DF65B0', 10:'#C994C7', 11:'#D4B9DA'}
    min_cn = int(cn_data.min())
    max_cn = int(cn_data.max())
    assert min_cn - cn_data.min() == 0
    assert max_cn - cn_data.max() == 0
    color_list = []
    for cn in range(min_cn, max_cn+1):
        if cn > max(color_reference.keys()):
            cn = max(color_reference.keys())
        color_list.append(color_reference[cn])
    return ListedColormap(color_list)


def get_methods_cmap():
    cmap = {
        'PERT clone': 'olive',
        'PERT comp.': 'yellowgreen',
        'PERT': 'yellowgreen',
        'Kronos': 'lightcoral',
        'laks': 'darksalmon',
        'Laks': 'darksalmon',
        'true': 'steelblue'
    }
    return cmap


def get_chrom_cmap():
    cmap = {
        'autosomes': '#CCCCCC', # grey
        'chrX': 'C1' # orange
    }
    return cmap


def get_bkpt_cmap():
    cmap = {
        'No': '#CCCCCC', # grey
        False: '#CCCCCC', # grey
        'Yes': 'C2', # green
        True: 'C2' # green
    }
    return cmap


def get_clone_cmap():
    cmap = {
        'A': 'cadetblue',
        1: 'cadetblue',
        'B': 'chocolate',
        2: 'chocolate',
        'C': 'olivedrab',
        3: 'olivedrab',
        'D': 'tan',
        4: 'tan',
        'E': 'plum',
        5: 'plum',
        'F': 'indianred',
        6: 'indianred',
        'G': 'lightpink',
        7: 'lightpink',
        'H': 'slategrey',
        8: 'slategrey',
        'I': 'darkseagreen',
        9: 'darkseagreen',
        'J': 'darkkhaki',
        10: 'darkkhaki',
        'K': 'lightsteelblue',
        11: 'lightsteelblue',
        'L': 'darksalmon',
        12: 'darksalmon',
        'M': 'lightgreen',
        13: 'lightgreen',
        'N': 'lightpink',
        14: 'lightpink',
        'O': 'lightgrey',
        15: 'lightgrey',
        'P': 'lightblue',
        16: 'lightblue',
        'Q': 'coral',
        17: 'coral',
        'R': 'lightcyan',
        18: 'lightcyan',
        'S': 'lightgoldenrodyellow',
        19: 'lightgoldenrodyellow',
        'T': 'darkseagreen',
        20: 'darkseagreen',
        'U': 'indigo',
        21: 'indigo'
    }
    return cmap


def get_htert_cmap():
    cmap = {
        'WT': 'C0',
        'SA039': 'C0',
        'TP53-/-': 'C1',
        'SA906a': 'C1',
        'SA906b': 'orange',
        'TP53-/-,BRCA1+/-' : 'C2',
        'SA1292': 'C2',
        'TP53-/-,BRCA1-/-': 'C3',
        'SA1056': 'C3',
        'TP53-/-,BRCA2+/-': 'C4',
        'SA1188': 'C4',
        'TP53-/-,BRCA2-/-': 'C5',
        'SA1054': 'C5',
        'SA1055': 'chocolate',
        'OV2295': 'lightgreen'
    }
    return cmap


def get_facs_cmap():
    ''' Global color map for FACS isolated cell types '''
    cmap = {
        # FACS sorted cell lines
        'GM18507': 'violet', 'SA928': 'violet', 1: 'violet',
        'T47D': 'gold', 'SA1044': 'gold', 2: 'gold',
    }
    return cmap


def get_metacohort_cmaps():
    cell_type_cdict = {
        'hTERT': 'lightsteelblue', 0: 'lightsteelblue',
        'HGSOC': 'teal', 1: 'teal',
        'TNBC': 'salmon', 2: 'salmon',
        'OV2295': 'lightgreen', 3: 'lightgreen',
        'T47D': 'gold', 4: 'gold',
        'GM18507': 'violet', 5: 'violet',
    }
    cell_type_cmap = LinearSegmentedColormap.from_list('cell_type_cmap', list(cell_type_cdict.values()), N=len(cell_type_cdict))

    signature_cdict = {
        'FBI': 'plum', 0: 'plum',
        'HRD': 'cyan', 1: 'cyan',
        'TD': 'coral', 2: 'coral',
    }
    signature_cmap = LinearSegmentedColormap.from_list('signature_cmap', list(signature_cdict.values()), N=len(signature_cdict))

    condition_cdict = {
        'Line': 'tan', 0: 'tan',
        'PDX': 'lightskyblue', 1: 'lightskyblue',
    }
    condition_cmap = LinearSegmentedColormap.from_list('condition_cmap', list(condition_cdict.values()), N=len(condition_cdict))

    ploidy_cdict = {2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33'}
    ploidy_cmap = LinearSegmentedColormap.from_list('ploidy_cmap', list(ploidy_cdict.values()), N=len(ploidy_cdict))

    return cell_type_cmap, signature_cmap, condition_cmap, ploidy_cmap


def get_rx_cmap():
    cmap = {
        'Rx-': '#CCCCCC', # grey
        'U': '#CCCCCC',
        'untreated': '#CCCCCC',
        'Rx+': 'C0',
        'T': 'C0',
        'treated': 'C0'
    }
    return cmap
