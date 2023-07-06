import os
import logging
import torch
import numpy as np
import pandas as pd
from argparse import ArgumentParser

import pyro
import pyro.distributions as dist
import pyro.distributions.transforms as T
from torch import nn
from pyro.nn import PyroModule

# for CI testing
smoke_test = ('CI' in os.environ)
pyro.set_rng_seed(1)
torch.manual_seed(1)
np.random.seed(1)

assert issubclass(PyroModule[nn.Linear], nn.Linear)
assert issubclass(PyroModule[nn.Linear], PyroModule)


def get_args():
    p = ArgumentParser()

    p.add_argument('-cr', '--clone_rt', type=str, help='matrix of RT samples for each clone')
    p.add_argument('-f', '--features', type=str, help='clone feature metadata (e.g. signature, ploidy, cell type)')
    p.add_argument('--remove_x', action='store_true', default=False, help='filter out all loci with chr=="X"')
    p.add_argument('-bip', '--beta_importance_posteriors', type=str, help='beta importance posteriors')
    p.add_argument('-rtp', '--rt_profile_posteriors', type=str, help='rt profile posteriors')
    p.add_argument('-lc', '--loss_curve', type=str, help='loss curve')

    return p.parse_args()


def preprocess_data(clone_features, clone_rt):
    # set loci as index and create a tensor of RT values
    rt_data = torch.tensor(clone_rt.values).type(torch.float64)

    # convert cell types to a tensor of floats
    cell_types = torch.tensor(clone_features[[
        'type_HGSOC',
        'type_TNBC',
        'type_hTERT',
        'type_T47D',
        'type_GM18507',
        'type_OV2295'
    ]].values).type(torch.float64)

    # convert signatures to a tensor of floats
    signatures = torch.tensor(clone_features[[
        'signature_FBI',
        'signature_HRD',
        'signature_TD'
    ]].values).type(torch.float64)

    # binarize the clone_features 'ploidy' column into 0s where ploidy<= 2 and 1s where ploidy>2
    clone_features['wgd'] = clone_features['ploidy'].apply(lambda x: 1 if x>2.0 else 0)
    clone_features['ngd'] = clone_features['ploidy'].apply(lambda x: 1 if x<=2.0 else 0)
    # convert wgd to a tensor of floats
    wgd = torch.tensor(clone_features[[
        'wgd', 
        'ngd'
    ]].values).type(torch.float64)

    return rt_data, cell_types, signatures, wgd


def model(cell_types, signatures, wgd, rt_data=None, num_loci=None, num_profiles=None):
    """
    :param cell_types: tensor of cell type labels for each profile
    :param signatures: tensor of signature labels for each profile
    :param wgd: tensor of WGD labels for each profile
    :param rt_data: tensor of RT values for each locus in each profile
    :param num_loci: number of loci in the RT data (when rt_data is None)
    :param num_profiles: number of profiles in the RT data (when rt_data is None)
    """
    assert (rt_data is not None) or (num_loci is not None)
    if rt_data is not None:
        num_loci, num_profiles = rt_data.shape
    assert num_loci is not None
    assert num_profiles is not None
    assert cell_types.shape == (num_profiles, 6)  # HGSOC, TNBC, hTERT, T47D, GM18507, OV2295
    assert signatures.shape == (num_profiles, 3)  # FBI, HRD, TD 
    assert wgd.shape == (num_profiles, 2)  # WGD, NGD

    # define profile and loci plates
    loci_plate = pyro.plate('num_loci', num_loci, dim=-2)
    profile_plate = pyro.plate('num_profiles', num_profiles, dim=-1)

    # define importance terms for all the variables in the regression model
    # these should all be on the positive real domain so that positive rt values correspond to early rt in rt_data
    beta_global_rt_importance = pyro.sample('beta_global_rt_importance', dist.Gamma(1., 1.))
    beta_cell_type_importance = pyro.sample('beta_cell_type_importance', dist.Gamma(1., 1.))
    beta_signature_importance = pyro.sample('beta_signature_importance', dist.Gamma(1., 1.))
    beta_wgd_importance = pyro.sample('beta_wgd_importance', dist.Gamma(1., 1.))

    # global replication timing
    global_rt = pyro.sample('global_rt', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))

    # cell type specific replication timing
    ct_hgsoc = pyro.sample('ct_hgsoc', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ct_tnbc = pyro.sample('ct_tnbc', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ct_htert = pyro.sample('ct_htert', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ct_t47d = pyro.sample('ct_t47d', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ct_gm18507 = pyro.sample('ct_gm18507', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ct_ov2295 = pyro.sample('ct_ov2295', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))

    # get signature-specific replication timing
    sig_fbi = pyro.sample('sig_fbi', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    sig_hrd = pyro.sample('sig_hrd', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    sig_td = pyro.sample('sig_td', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))

    # whole genome doubled specific replication timing
    wgd_rt = pyro.sample('wgd_rt', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))
    ngd_rt = pyro.sample('ngd_rt', dist.Normal(torch.zeros(num_loci), torch.ones(num_loci)).to_event(1))

    # normalize all the replication timing distributions
    global_rt_norm = T.Normalize()(global_rt)
    ct_hgsoc_norm = T.Normalize()(ct_hgsoc)
    ct_tnbc_norm = T.Normalize()(ct_tnbc)
    ct_htert_norm = T.Normalize()(ct_htert)
    ct_t47d_norm = T.Normalize()(ct_t47d)
    ct_gm18507_norm = T.Normalize()(ct_gm18507)
    ct_ov2295_norm = T.Normalize()(ct_ov2295)
    sig_fbi_norm = T.Normalize()(sig_fbi)
    sig_hrd_norm = T.Normalize()(sig_hrd)
    sig_td_norm = T.Normalize()(sig_td)
    wgd_rt_norm = T.Normalize()(wgd_rt)
    ngd_rt_norm = T.Normalize()(ngd_rt)

    # random variable to represent the standard deviation when sampling 
    sigma = pyro.sample("sigma", dist.Uniform(0., 1.))

    with profile_plate:

        # add a term that captures the average RT value for each profile
        # this should compensate for samples which are drastically shifted to early or late replication
        # use a Beta(1, 1) prior to keep the distribution uniform on [0, 1]
        # theta = pyro.sample('expose_theta', dist.Normal(torch.tensor([0.]), torch.tensor([1.])))
        
        # multiply the beta of each term by its normalized RT profile
        # start with the global RT profile
        term_global_rt = beta_global_rt_importance * global_rt_norm.reshape(-1, 1)
        
        # cell type RT profiles
        term_hgsoc = beta_cell_type_importance * (ct_hgsoc_norm.reshape(-1, 1) * cell_types[:, 0])
        term_tnbc = beta_cell_type_importance * (ct_tnbc_norm.reshape(-1, 1) * cell_types[:, 1])
        term_htert = beta_cell_type_importance * (ct_htert_norm.reshape(-1, 1) * cell_types[:, 2])
        term_t47d = beta_cell_type_importance * (ct_t47d_norm.reshape(-1, 1) * cell_types[:, 3])
        term_gm18507 = beta_cell_type_importance * (ct_gm18507_norm.reshape(-1, 1) * cell_types[:, 4])
        term_ov2295 = beta_cell_type_importance * (ct_ov2295_norm.reshape(-1, 1) * cell_types[:, 5])

        # signature RT profiles
        term_fbi = beta_signature_importance * (sig_fbi_norm.reshape(-1, 1) * signatures[:, 0])
        term_hrd = beta_signature_importance * (sig_hrd_norm.reshape(-1, 1) * signatures[:, 1])
        term_td = beta_signature_importance * (sig_td_norm.reshape(-1, 1) * signatures[:, 2])

        # WGD RT profile
        term_wgd = beta_wgd_importance * (wgd_rt_norm.reshape(-1, 1) * wgd[:, 0])
        term_ngd = beta_wgd_importance * (ngd_rt_norm.reshape(-1, 1) * wgd[:, 1])

        # sum all terms together to generate a mean on the real domain
        mean = term_global_rt + term_hgsoc + term_tnbc + term_htert + term_t47d + term_gm18507 + term_ov2295 + term_fbi + term_hrd + term_td + term_wgd + term_ngd

        # sigmoid transform the mean so it lies on the 0-1 domain
        new_mean = torch.sigmoid(mean)

        with loci_plate:
            # sample the observed RT data using the transformed mean and sigma
            return pyro.sample("obs", dist.Normal(new_mean, sigma), obs=rt_data)

def main():
    argv = get_args()

    # load the input data
    clone_features = pd.read_csv(argv.features)
    clone_rt = pd.read_csv(argv.clone_rt)

    # Optional: filter out all loci with chr=='X'
    if argv.remove_x:
        clone_rt = clone_rt.query('chr!="X"')
    
    # set clone_rt as the index
    clone_rt = clone_rt.set_index(['chr', 'start', 'end'])

    # preprocess the data into torch tensors
    rt_data, cell_types, signatures, wgd = preprocess_data(clone_features, clone_rt)

    # define the guide function according to the model using AutoNormal
    # initialize the SVI object and the optimizer
    auto_guide = pyro.infer.autoguide.AutoNormal(model)
    adam = pyro.optim.Adam({"lr": 0.02})
    elbo = pyro.infer.Trace_ELBO()
    svi = pyro.infer.SVI(model, auto_guide, adam, elbo)

    # train the model
    losses = []
    for step in range(1000 if not smoke_test else 2):  # Consider running for more steps.
        loss = svi.step(cell_types, signatures, wgd, rt_data, None, None)
        losses.append(loss)
        if step % 1 == 0:
            logging.info("Elbo loss: {}".format(loss))
    
    losses = []
    loss_diffs = []
    frac_loss_diffs = []
    total_loss_diff = np.nan
    converged = False
    for step in range(1000 if not smoke_test else 2):
        loss = svi.step(cell_types, signatures, wgd, rt_data, None, None)
        losses.append(loss)
        
        # check for convergence after 250 steps
        # see if the ratio between the total loss difference (first and last) and 
        # the the average of the 10 most recent loss differences is less than 1e-3
        if step > 250:
            total_loss_diff = abs(losses[-1] - losses[0])
            temp_loss_diff = []
            # loop over the past 10 losses and find their difference from the current loss
            for i in range(10):
                temp_loss_diff.append(abs(losses[-1] - losses[-i-2]))
            # take the mean of the 5 differences, the smaller this is, the closer the losses are to converging
            mean_temp_loss_diff = np.mean(temp_loss_diff)
            loss_diffs.append(mean_temp_loss_diff)
            frac_loss_diff = mean_temp_loss_diff / total_loss_diff
            frac_loss_diffs.append(frac_loss_diff)
            if frac_loss_diff < 1e-3:
                converged = True
        if converged:
            print('Converged at step {}'.format(step))
            break

    # create a dataframe for the iteration and loss values
    loss_df = pd.DataFrame({
        'iteration': np.arange(len(losses)),
        'loss': losses
    })

    # draw multiple samples in parallel from our trained guide
    # these will form our posterior distributions
    with pyro.plate("samples", 800, dim=-3):
        samples = auto_guide(cell_types, signatures, wgd)
    
    # create a dataframe for the posterior RT profiles
    # create an output df with the same index as clone_df
    output_df = pd.DataFrame(index=clone_rt.index)
    # add a column for the each parallel sample of expose_rho
    for i in range(samples['global_rt'].shape[0]):
        temp_a = samples['global_rt'][i, 0, 0, :]
        temp_a_norm = T.Normalize()(temp_a)
        output_df['global_rt_{}'.format(i)] = temp_a_norm.detach().numpy()
    # do the same for all the other parameters which are along the expose_loci plate
    for name, value in samples.items():
        if name.startswith('ct_') or name.startswith('sig_') or name.startswith('wgd_') or name.startswith('ngd_'):
            for i in range(value.shape[0]):
                temp_rt = value[i, 0, 0, :]
                temp_rt_norm = T.Normalize()(temp_rt)
                output_df['{}_{}'.format(name, i)] = temp_rt_norm.detach().numpy()
    output_df.reset_index(inplace=True)

    # create an output dataframe for the beta importance terms
    beta_importance_posteriors = pd.DataFrame({
        'beta_cell_type_importance': np.abs(samples['beta_cell_type_importance'][:, 0, 0].detach().numpy()),
        'beta_signature_importance': np.abs(samples['beta_signature_importance'][:, 0, 0].detach().numpy()),
        'beta_wgd_importance': np.abs(samples['beta_wgd_importance'][:, 0, 0].detach().numpy()),
        'beta_global_rt_importance': np.abs(samples['beta_global_rt_importance'][:, 0, 0].detach().numpy())
    })

    # save the loss curve
    loss_df.to_csv(argv.loss_curve, index=False)

    # save the beta importance posteriors
    beta_importance_posteriors.to_csv(argv.beta_importance_posteriors, index=False)

    # save the RT profile posteriors
    output_df.to_csv(argv.rt_profile_posteriors, index=False)


if __name__ == '__main__':
    main()
