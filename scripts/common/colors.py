from matplotlib.colors import ListedColormap

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


def get_cell_line_cmap():
    ''' Global color map for cell lines '''
    cmap = {
        # orange for GM18507
        'GM18507': 'C1',
        1: 'C1',
        # blue for T47D
        'T47D': 'C0',
        2: 'C0'
    }
    return cmap


def get_cna_cmap():
    ''' Global color map for copy number alterations '''
    cmap = {
        'gain': 'red',  # red
        'loss': 'deepskyblue',  # dark blue
        'neutral': '#CCCCCC'  # grey
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
        'PERT clone': 'lightblue',  # POWDERKEG BLUE
        'PERT comp.': 'dodgerblue',  # TRUE BLUE
        'PERT': 'dodgerblue',  # TRUE BLUE
        'Kronos': 'goldenrod',  # GOLD
        'laks': 'gold',  # YELLOW
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
        'A': 'C0',
        1: 'C0',
        'B': 'C1',
        2: 'C1',
        'C': 'C2',
        3: 'C2',
        'D': 'C3',
        4: 'C3',
        'E': 'C4',
        5: 'C4',
        'F': 'C5',
        6: 'C5',
        'G': 'C6',
        7: 'C6',
        'H': 'C7',
        8: 'C7',
        'I': 'C8',
        9: 'C8',
        'J': 'C9',
        10: 'C9',
        'K': 'lightsteelblue',
        11: 'lightsteelblue',
        'L': 'darksalmon',
        12: 'darksalmon',
        'M': 'lightgreen',
        13: 'lightgreen'
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
        'OV2295': 'C6'
    }
    return cmap