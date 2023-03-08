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
        'gain': '#BA0021',  # red
        'loss': '#003263',  # dark blue
        'neutral': '#C4CED4'  # silver
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
        'autosomes': 'lightsteelblue',
        'chrX': 'salmon'
    }
    return cmap
