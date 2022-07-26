import numpy as np
import pandas as pd
from argparse import ArgumentParser

import torch
from torch.distributions import constraints

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.util import ignore_jit_warnings


def get_args():
    p = ArgumentParser()

    p.add_argument('-si', '--df_s', help='True somatic CN profiles for S-phase cells')
    p.add_argument('-gi', '--df_g', help='True somatic CN profiles for S-phase cells')
    p.add_argument('-n', '--num_reads', type=int, help='number of reads per cell')
    p.add_argument('-nbr', '--nb_r', type=int, help='negative binomial rate for overdispersion')
    p.add_argument('-a', '--a', type=int, help='amplitude of sigmoid curve when generating rt noise')
    p.add_argument('-b', '--betas', type=float, nargs='+', help='list of beta coefficients for gc bias')
    p.add_argument('-rt', '--rt_col', type=str, help='column in cn input containing appropriate replication timing values')
    p.add_argument('-gc', '--gc_col', type=str, help='column in cn input containing GC content of each bin')
    p.add_argument('-so', '--s_out', help='simulated S-phase cells')
    p.add_argument('-go', '--g_out', help='simulated G1/2-phase cells')

    return p.parse_args()


def make_gc_features(x, poly_degree):
    """Builds features i.e. a matrix with columns [x, x^2, x^3, x^4]."""
    x = x.unsqueeze(1)
    return torch.cat([x ** i for i in reversed(range(0, poly_degree+1))], 1)


def model_s(gc_profile, cn0=None, rt0=None, num_cells=None, num_loci=None, cn_prior=None, u_guess=70., nb_r_guess=10000., poly_degree=4):
    with ignore_jit_warnings():
        if cn0 is not None:
            num_loci, num_cells = cn0.shape
        assert num_cells is not None
        assert num_loci is not None

    # controls the consistency of replicating on time
    a = pyro.sample('expose_a', dist.Gamma(torch.tensor([2.]), torch.tensor([0.2])))
    
    # negative binomial dispersion
    nb_r = pyro.param('expose_nb_r', torch.tensor([nb_r_guess]), constraint=constraints.positive)
    
    # define cell and loci plates
    loci_plate = pyro.plate('num_loci', num_loci, dim=-2)
    cell_plate = pyro.plate('num_cells', num_cells, dim=-1)

    betas = pyro.sample('expose_betas', dist.Normal(0., 1.).expand([poly_degree+1]).to_event(1))

    if rt0 is not None:
        # fix rt as constant when input into model
        rt = rt0
    else:
        with loci_plate:
            # bulk replication timing profile
            rt = pyro.sample('expose_rt', dist.Beta(torch.tensor([1.]), torch.tensor([1.])))

    # fix cn as constant when input into model
    if cn0 is not None:
        cn = cn0

    with cell_plate:

        # per cell replication time
        time = pyro.sample('expose_time', dist.Beta(torch.tensor([1.]), torch.tensor([1.])))

        # per cell reads per copy per bin
        u = pyro.sample('expose_u', dist.Normal(torch.tensor([u_guess]), torch.tensor([u_guess/10.])))
        
        with loci_plate:

            if cn0 is None:
                if cn_prior is None:
                    cn_prior = torch.ones(num_loci, num_cells, 13)
                # sample cn probabilities of each bin from Dirichlet
                cn_prob = pyro.sample('expose_cn_prob', dist.Dirichlet(cn_prior))
                # sample cn state from categorical based on cn_prob
                cn = pyro.sample('cn', dist.Categorical(cn_prob), infer={"enumerate": "parallel"})

            # per cell per bin late or early 
            t_diff = time.reshape(-1, num_cells) - rt.reshape(num_loci, -1)

            # probability of having been replicated
            p_rep = 1 / (1 + torch.exp(-a * t_diff))

            # binary replicated indicator
            rep = pyro.sample('rep', dist.Bernoulli(p_rep), infer={"enumerate": "parallel"})

            # copy number accounting for replication
            rep_cn = cn * (1. + rep)

            # copy number accounting for gc bias
            gc_features = make_gc_features(gc_profile, poly_degree)
            gc_rate = torch.exp(torch.sum(betas * gc_features, 1))
            biased_cn = rep_cn * gc_rate.reshape(-1, 1)

            # expected reads per bin per cell
            expected_reads = (u * biased_cn)

            nb_p = expected_reads / (expected_reads + nb_r)
            
            reads = pyro.sample('reads', dist.NegativeBinomial(nb_r, probs=nb_p), obs=None)

    return reads


def model_g1(gc_profile, cn=None, num_cells=None, num_loci=None, u_guess=70., poly_degree=4):
    with ignore_jit_warnings():
        if cn is not None:
            num_loci, num_cells = cn.shape
        assert num_cells is not None
        assert num_loci is not None
    
    # negative binomial dispersion
    nb_r = pyro.param('expose_nb_r', torch.tensor([10000.0]), constraint=constraints.positive)

    # gc bias params
    betas = pyro.sample('expose_betas', dist.Normal(0., 1.).expand([poly_degree+1]).to_event(1))

    with pyro.plate('num_cells', num_cells):

        # per cell reads per copy per bin
        u = pyro.sample('expose_u', dist.Normal(torch.tensor([u_guess]), torch.tensor([u_guess/10.])))

        with pyro.plate('num_loci', num_loci):

            # copy number accounting for gc bias
            gc_features = make_gc_features(gc_profile, poly_degree)
            gc_rate = torch.exp(torch.sum(betas * gc_features, 1))
            biased_cn = cn * gc_rate.reshape(-1, 1)

            # expected reads per bin per cell
            expected_reads = (u * biased_cn)

            nb_p = expected_reads / (expected_reads + nb_r)

            reads = pyro.sample('reads', dist.NegativeBinomial(nb_r, probs=nb_p), obs=None)

    return reads


def convert_rt_units(rt):
    # make sure rt units range from 0-1
    return (rt - rt.min()) / (rt.max() - rt.min())


def simulate_s_cells(gc_profile, cn, rt, argv):
    pyro.clear_param_store()

    num_loci, num_cells = cn.shape
    gc_profile = torch.tensor(gc_profile.values)
    cn = torch.tensor(cn.values)
    rt_profile = torch.tensor(convert_rt_units(rt.values))

    u_guess = float(argv.num_reads) / (1.5 * torch.mean(cn))

    print('num_loci', num_loci)
    print('num_cells', num_cells)
    print('gc_profile', gc_profile.shape)
    print('cn', cn.shape)
    print('rt_profile', rt_profile.shape)
    print('u_guess', u_guess)

    conditioned_model = poutine.condition(
        model_s,
        data={
            'expose_a': torch.tensor([argv.a]),
            'expose_nb_r': torch.tensor([argv.nb_r]),
            'expose_betas': torch.tensor(argv.betas),
            'expose_rt': rt_profile
        })

    model_trace = pyro.poutine.trace(conditioned_model)

    samples = model_trace.get_trace(gc_profile, cn0=cn, u_guess=u_guess, poly_degree=len(argv.betas)-1)

    t = samples.nodes['expose_time']['value']
    u = samples.nodes['expose_u']['value']

    t_diff = t.reshape(-1, num_cells) - rt_profile.reshape(num_loci, -1)
    p_rep = 1 / (1 + torch.exp(-argv.a * t_diff))

    rep = samples.nodes['rep']['value']

    rep_cn = cn * (1. + rep)

    reads = samples.nodes['reads']['value']

    # normalize read count
    reads_norm = (reads / torch.sum(reads, 0)) * argv.num_reads
    reads_norm = reads_norm.type(torch.int64)

    return reads_norm, reads, rep, p_rep, t


def simulate_g_cells(gc_profile, cn, argv):
    pyro.clear_param_store()

    num_loci, num_cells = cn.shape
    gc_profile = torch.tensor(gc_profile.values)
    cn = torch.tensor(cn.values)

    u_guess = float(argv.num_reads) / (1. * torch.mean(cn))

    conditioned_model = poutine.condition(
        model_g1,
        data={
            'expose_nb_r': torch.tensor([argv.nb_r]),
            'expose_betas': torch.tensor(argv.betas),
        })

    model_trace = pyro.poutine.trace(conditioned_model)

    samples = model_trace.get_trace(gc_profile, cn=cn, u_guess=u_guess, poly_degree=len(argv.betas)-1)

    u = samples.nodes['expose_u']['value']

    reads = samples.nodes['reads']['value']

    # normalize read count
    reads_norm = (reads / torch.sum(reads, 0)) * argv.num_reads
    reads_norm = reads_norm.type(torch.int64)

    return reads_norm, reads


def main():
    argv = get_args()

    df_s = pd.read_csv(argv.df_s, sep='\t')
    df_g = pd.read_csv(argv.df_g, sep='\t')

    df_s.chr = df_s.chr.astype(str)
    df_g.chr = df_g.chr.astype(str)

    cn_s = pd.pivot_table(df_s, index=['chr', 'start'], columns='cell_id', values='true_G1_state')
    cn_g = pd.pivot_table(df_g, index=['chr', 'start'], columns='cell_id', values='true_G1_state')
    gc_profile = df_s[['chr', 'start', argv.gc_col]].drop_duplicates()[argv.gc_col]
    rt_profile = df_s[['chr', 'start', argv.rt_col]].drop_duplicates()[argv.rt_col]

    # S-phase: condition each model based in argv parameters and simulate read count
    reads_norm, reads, rep, p_rep, t = simulate_s_cells(gc_profile, cn_s, rt_profile, argv)

    reads_norm_df = pd.DataFrame(reads_norm.numpy(), columns=cn_s.columns, index=cn_s.index)
    reads_df = pd.DataFrame(reads.numpy(), columns=cn_s.columns, index=cn_s.index)
    rep_df = pd.DataFrame(rep.numpy(), columns=cn_s.columns, index=cn_s.index)
    p_rep_df = pd.DataFrame(p_rep.numpy(), columns=cn_s.columns, index=cn_s.index)
    t_df = pd.DataFrame(t.numpy(), columns=['true_t'], index=cn_s.columns)

    print('reads_norm_df\n', reads_norm_df.head())
    print('reads_df\n', reads_df.head())
    print('rep_df\n', rep_df.head())
    print('p_rep_df\n', p_rep_df.head())
    print('t_df\n', t_df.head())

    # merge normalized read count
    reads_norm_df = reads_norm_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_reads_norm')
    reads_norm_df.chr = reads_norm_df.chr.astype(str)
    df_s = pd.merge(df_s, reads_norm_df)
    # merge raw reads before normalizing total read count
    reads_df = reads_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_reads_raw')
    reads_df.chr = reads_df.chr.astype(str)
    df_s = pd.merge(df_s, reads_df)
    # merge true replication states
    rep_df = rep_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_rep')
    rep_df.chr = rep_df.chr.astype(str)
    df_s = pd.merge(df_s, rep_df)
    # merge probability of each bin being replicated
    p_rep_df = p_rep_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_p_rep')
    p_rep_df.chr = p_rep_df.chr.astype(str)
    df_s = pd.merge(df_s, p_rep_df)
    # merge s-phase times
    df_s = pd.merge(df_s, t_df.reset_index(), on='cell_id')

    # G1-phase: condition each model based in argv parameters and simulate read count
    reads_norm_g, reads_g = simulate_g_cells(gc_profile, cn_g, argv)
    
    reads_norm_g_df = pd.DataFrame(reads_norm_g.numpy(), columns=cn_g.columns, index=cn_g.index)
    reads_g_df = pd.DataFrame(reads_g.numpy(), columns=cn_g.columns, index=cn_g.index)
    
    # merge normalized read count
    reads_norm_g_df = reads_norm_g_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_reads_norm')
    reads_norm_g_df.chr = reads_norm_g_df.chr.astype(str)
    df_g = pd.merge(df_g, reads_norm_g_df)
    # merge raw reads before normalizing total read count
    reads_g_df = reads_g_df.reset_index().melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='true_reads_raw')
    reads_g_df.chr = reads_g_df.chr.astype(str)
    df_g = pd.merge(df_g, reads_g_df)

    df_s.to_csv(argv.s_out, sep='\t', index=False)
    df_g.to_csv(argv.g_out, sep='\t', index=False)



if __name__ == '__main__':
    main()