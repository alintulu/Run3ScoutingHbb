from __future__ import print_function, division
import sys, os
import csv, json
import numpy as np
from scipy.interpolate import interp1d
import scipy.stats
import pickle
import ROOT
import pandas as pd

import rhalphalib as rl
from rhalphalib import AffineMorphTemplate, MorphHistW2

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

eps=0.0000001
do_muon_CR = False

# Tell me if my sample is too small to care about
def badtemp_ma(hvalues, mask=None):
    # Need minimum size & more than 1 non-zero bins                                                                                                                                   
    tot = np.sum(hvalues[mask])

    count_nonzeros = np.sum(hvalues[mask] > 0)
    if (tot < eps) or (count_nonzeros < 2):
        return True
    else:
        return False

# Turn an reg distribution into a single bin (for muon CR)
def one_bin(template):
    try:
        h_vals, h_edges, h_key, h_variances = template
        return (np.array([np.sum(h_vals)]), np.array([0., 1.]), "onebin", np.array([np.sum(h_variances)]))
    except:
        h_vals, h_edges, h_key = template
        return (np.array([np.sum(h_vals)]), np.array([0., 1.]), "onebin")

def shape_to_num(var, nom, clip=2):
    nom_rate = np.sum(nom)
    var_rate = np.sum(var)

    if abs(var_rate/nom_rate) > clip:
        var_rate = clip*nom_rate

    if var_rate < 0:
        var_rate = 0

    return var_rate/nom_rate

# Read the histogram
def get_template(sName, passed, ptbin, cat, obs, syst, muon=False):
    """
    Read reg template from root file
    """

    f = ROOT.TFile.Open('signalregion.root')

    name = cat+'fail_'
    if passed:
        name = cat+'pass_'
    if cat == 'ggf_':
        name += 'pt'+str(ptbin)+'_'

    name += sName+'_'+syst

    print(name)

    h = f.Get(name)

    sumw = []
    sumw2 = []

    for i in range(1,h.GetNbinsX()+1):
        sumw += [h.GetBinContent(i)]
        sumw2 += [h.GetBinError(i)*h.GetBinError(i)]

    return (np.array(sumw), obs.binning, obs.name, np.array(sumw2))

# Plots of the MC transfer factor
def plot_mctf(tf_MCtempl, regbins, name):
    """
    Plot the MC pass / fail TF as function of (pt,rho) and (pt,reg)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # arrays for plotting pt vs reg                    
    pts = np.linspace(300,1200,15)
    ptpts, regpts = np.meshgrid(pts[:-1] + 0.5 * np.diff(pts), regbins[:-1] + 0.5 * np.diff(regbins), indexing='ij')
    ptpts_scaled = (ptpts - 300.) / (1200. - 300.)
    rhopts = 2*np.log(regpts/ptpts)

    rhopts_scaled = (rhopts - (-6)) / ((-1.7) - (-6))
    validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

    ptpts = ptpts[validbins]
    regpts = regpts[validbins]
    ptpts_scaled = ptpts_scaled[validbins]
    rhopts_scaled = rhopts_scaled[validbins]

    tf_MCtempl_vals = tf_MCtempl(ptpts_scaled, rhopts_scaled, nominal=True)
    df = pd.DataFrame([])
    df['reg'] = regpts.reshape(-1)
    df['pt'] = ptpts.reshape(-1)
    df['MCTF'] = tf_MCtempl_vals.reshape(-1)

    fig, ax = plt.subplots()
    h = ax.hist2d(x=df["reg"],y=df["pt"],weights=df["MCTF"], bins=(regbins,pts))
    plt.xlabel("$m_{reg}$ [GeV]")
    plt.ylabel("$p_{T}$ [GeV]")
    cb = fig.colorbar(h[3],ax=ax)
    cb.set_label("Ratio")
    fig.savefig("plots/MCTF_regpt_"+name+".png",bbox_inches="tight")
    fig.savefig("plots/MCTF_regpt_"+name+".pdf",bbox_inches="tight")
    plt.clf()

    # arrays for plotting pt vs rho                                          
    rhos = np.linspace(-6,-1.7,23)
    ptpts, rhopts = np.meshgrid(pts[:-1] + 0.5*np.diff(pts), rhos[:-1] + 0.5 * np.diff(rhos), indexing='ij')
    ptpts_scaled = (ptpts - 300.) / (1200. - 300.)
    rhopts_scaled = (rhopts - (-6)) / ((-1.7) - (-6))
    validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

    ptpts = ptpts[validbins]
    rhopts = rhopts[validbins]
    ptpts_scaled = ptpts_scaled[validbins]
    rhopts_scaled = rhopts_scaled[validbins]

    tf_MCtempl_vals = tf_MCtempl(ptpts_scaled, rhopts_scaled, nominal=True)

    df = pd.DataFrame([])
    df['rho'] = rhopts.reshape(-1)
    df['pt'] = ptpts.reshape(-1)
    df['MCTF'] = tf_MCtempl_vals.reshape(-1)

    fig, ax = plt.subplots()
    h = ax.hist2d(x=df["rho"],y=df["pt"],weights=df["MCTF"],bins=(rhos,pts))
    plt.xlabel("rho")
    plt.ylabel("$p_{T}$ [GeV]")
    cb = fig.colorbar(h[3],ax=ax)
    cb.set_label("Ratio")
    fig.savefig("plots/MCTF_rhopt_"+name+".png",bbox_inches="tight")
    fig.savefig("plots/MCTF_rhopt_"+name+".pdf",bbox_inches="tight")

    return

def example_rhalphabet(tmpdir,
                       throwPoisson = True,
                       fast=0,
                       year="2022"):
    """ 
    Create the data cards!
    """

    # Systematics
    sys_dict = {}

    sys_dict['JES'] = rl.NuisanceParameter('CMS_scale_j_{}'.format(year), 'lnN')
    sys_dict['JER'] = rl.NuisanceParameter('CMS_res_j_{}'.format(year), 'lnN')
    exp_systs = ['JES','JER']

    # Simple lumi systematics                                                                                                                                                            
    sys_lumi_uncor = rl.NuisanceParameter('CMS_lumi_13TeV_2016', 'lnN')

    # define bins    
    ptbins = {}
    ptbins['ggf'] = np.array([350, 400, 450, 500, 1200])
    #ptbins['ggf'] = np.array([300, 350, 400, 450, 500, 550, 600, 675, 800, 1200])
    #ptbins['ggf'] = np.array([300, 1200])

    npt = {}
    npt['ggf'] = len(ptbins['ggf']) - 1

    regbins = np.linspace(40, 201, 24)
    reg = rl.Observable('reg', regbins)

    validbins = {}

    cats = ['ggf']
    ncat = len(cats)

    # Build qcd MC pass+fail model and fit to polynomial
    tf_params = {}
    for cat in cats:

        fitfailed_qcd = 0

        # here we derive these all at once with 2D array                            
        ptpts, regpts = np.meshgrid(ptbins[cat][:-1] + 0.3 * np.diff(ptbins[cat]), regbins[:-1] + 0.5 * np.diff(regbins), indexing='ij')
        rhopts = 2*np.log(regpts/ptpts)
        ptscaled = (ptpts - 300.) / (1200. - 300.)
        rhoscaled = (rhopts - (-6)) / ((-1.7) - (-6))
        validbins[cat] = (rhoscaled >= 0) & (rhoscaled <= 1)
        rhoscaled[~validbins[cat]] = 1  # we will mask these out later   

        n_max_fail = 100
        while fitfailed_qcd < n_max_fail:
        
            qcdmodel = rl.Model('qcdmodel_'+cat)
            qcdpass, qcdfail = 0., 0.

            for ptbin in range(npt[cat]):

                failCh = rl.Channel('ptbin%d%s%s%s' % (ptbin, cat, 'fail',year))
                passCh = rl.Channel('ptbin%d%s%s%s' % (ptbin, cat, 'pass',year))
                qcdmodel.addChannel(failCh)
                qcdmodel.addChannel(passCh)

                # QCD templates from file                           
                failTempl = get_template('QCD', 0, ptbin+1, cat+'_', obs=reg, syst='nominal')
                passTempl = get_template('QCD', 1, ptbin+1, cat+'_', obs=reg, syst='nominal')

                failCh.setObservation(failTempl, read_sumw2=True)
                passCh.setObservation(passTempl, read_sumw2=True)

                qcdfail += sum([val for val in failCh.getObservation()[0]])
                qcdpass += sum([val for val in passCh.getObservation()[0]])

            qcdeff = qcdpass / qcdfail
            print('Inclusive P/F from Monte Carlo = ' + str(qcdeff))

            # initial values                                                                 
            print('Initial fit values read from file initial_vals*')
            with open('initial_vals_'+cat+'.json') as f:
                initial_vals = np.array(json.load(f)['initial_vals'])
            print(initial_vals)

            tf_MCtempl = rl.BasisPoly("tf_MCtempl_"+cat+year,
                                      (initial_vals.shape[0]-1,initial_vals.shape[1]-1),
                                      ['pt', 'rho'], 
                                      init_params=initial_vals,
                                      basis='Bernstein',
                                      limits=(-20, 20), coefficient_transform=None)

            tf_MCtempl_params = qcdeff * tf_MCtempl(ptscaled, rhoscaled)

            for ptbin in range(npt[cat]):

                failCh = qcdmodel['ptbin%d%sfail%s' % (ptbin, cat, year)]
                passCh = qcdmodel['ptbin%d%spass%s' % (ptbin, cat, year)]
                failObs = failCh.getObservation()
                passObs = passCh.getObservation()
                
                qcdparams = np.array([rl.IndependentParameter('qcdparam_'+cat+'_ptbin%d_regbin%d' % (ptbin, i), 0) for i in range(reg.nbins)])
                sigmascale = 10.
                scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
                
                fail_qcd = rl.ParametericSample('ptbin%d%sfail%s_qcd' % (ptbin, cat, year), rl.Sample.BACKGROUND, reg, scaledparams[0])
                failCh.addSample(fail_qcd)
                pass_qcd = rl.TransferFactorSample('ptbin%d%spass%s_qcd' % (ptbin, cat, year), rl.Sample.BACKGROUND, tf_MCtempl_params[ptbin, :], fail_qcd)
                passCh.addSample(pass_qcd)
                
                failCh.mask = validbins[cat][ptbin]
                passCh.mask = validbins[cat][ptbin]

            qcdfit_ws = ROOT.RooWorkspace('w')

            simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
            qcdfit = simpdf.fitTo(obs,
                                  ROOT.RooFit.Extended(True),
                                  ROOT.RooFit.SumW2Error(True),
                                  ROOT.RooFit.Strategy(2),
                                  ROOT.RooFit.Save(),
                                  ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                                  ROOT.RooFit.PrintLevel(1),
                              )
            qcdfit_ws.add(qcdfit)
            qcdfit_ws.writeToFile(os.path.join(str(tmpdir), 'testModel_qcdfit_'+cat+'_'+year+'.root'))

            # Set parameters to fitted values  
            allparams = dict(zip(qcdfit.nameArray(), qcdfit.valueArray()))
            pvalues = []
            for i, p in enumerate(tf_MCtempl.parameters.reshape(-1)):
                p.value = allparams[p.name]
                pvalues += [p.value]
            
            if qcdfit.status() != 0:
                print('Could not fit qcd')
                fitfailed_qcd += 1

                new_values = np.array(pvalues).reshape(tf_MCtempl.parameters.shape)
                with open("initial_vals_"+cat+".json", "w") as outfile:
                    json.dump({"initial_vals":new_values.tolist()},outfile)

            else:
                break

        if fitfailed_qcd >=n_max_fail:
            #raise RuntimeError(f'Could not fit qcd after {n_max_fail} tries')
            print(f'Could not fit qcd after {n_max_fail} tries')

        print("Fitted qcd for category " + cat)

        # Plot the MC P/F transfer factor                                                   
        plot_mctf(tf_MCtempl, regbins, cat)                           

        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', qcdfit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

        # initial values                                                                                                                                         
        with open('initial_vals_data_'+cat+'.json') as f:
            initial_vals_data = np.array(json.load(f)['initial_vals'])

        tf_dataResidual = rl.BasisPoly("tf_dataResidual_"+year+cat,
                                       (initial_vals_data.shape[0]-1,initial_vals_data.shape[1]-1), 
                                       ['pt', 'rho'], 
                                       init_params=initial_vals_data,
                                       basis='Bernstein',
                                       limits=(-20,20), coefficient_transform=None)

        tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
        tf_params[cat] = qcdeff * tf_MCtempl_params_final * tf_dataResidual_params

    # build actual fit model now
    model = rl.Model('testModel_'+year)

    # exclude QCD from MC samps
    samps = ['ZJetsqq', 'ZJetsbb', 'TTbar', 'W', 'VV', 'VBF', 'ZH', 'WH', 'ttH', 'ggF']
    sigs = ['ggF']

    cols = ['bin','region','samp','syst','up','val']
    df = pd.DataFrame(columns=cols)

    for cat in cats:
        for ptbin in range(npt[cat]):
            for region in ['pass', 'fail']:

                print("Bin: " + cat + " bin " + str(ptbin) + " " + region)

                # drop bins outside rho validity                                                
                mask = validbins[cat][ptbin]
                failCh.mask = validbins[cat][ptbin]
                passCh.mask = validbins[cat][ptbin]

                ch = rl.Channel('ptbin%d%s%s%s' % (ptbin, cat, region, year))
                model.addChannel(ch)

                isPass = region == 'pass'
                templates = {}
            
                for sName in samps:

                    templates[sName] = get_template(sName, isPass, ptbin+1, cat+'_', obs=reg, syst='nominal')
                    nominal = templates[sName][0]

                    if(badtemp_ma(nominal)):
                        print("Sample {} is too small, skipping".format(sName))
                        continue

                    # expectations
                    templ = templates[sName]
                    
                    if sName in sigs:
                        stype = rl.Sample.SIGNAL
                    else:
                        stype = rl.Sample.BACKGROUND
                
                    sample = rl.TemplateSample(ch.name + '_' + sName, stype, templ)
                    sample.autoMCStats(lnN=True)

                    # lumi systematic
                    sample.setParamEffect(sys_lumi_uncor, 1.022)

                    # experimental systematics
                    for sys in exp_systs:
                        syst_up = get_template(sName, isPass, ptbin+1, cat+'_', obs=reg, syst=sys+'Up')[0]
                        syst_do = get_template(sName, isPass, ptbin+1, cat+'_', obs=reg, syst=sys+'Down')[0]

                        eff_up = shape_to_num(syst_up,nominal)
                        eff_do = shape_to_num(syst_do,nominal)

                        sample.setParamEffect(sys_dict[sys], eff_up, eff_do)
                        
                        df = df.append(pd.DataFrame([[cat+' '+str(ptbin+1),region,sName,sys,1,eff_up-1]],columns=cols))
                        df = df.append(pd.DataFrame([[cat+' '+str(ptbin+1),region,sName,sys,0,eff_do-1]],columns=cols))

                    ch.addSample(sample)

                data_obs = get_template('data', isPass, ptbin+1, cat+'_', obs=reg, syst='nominal')
                ch.setObservation(data_obs, read_sumw2=True)

    for cat in cats:
        for ptbin in range(npt[cat]):

                failCh = model['ptbin%d%sfail%s' % (ptbin, cat, year)]
                passCh = model['ptbin%d%spass%s' % (ptbin, cat, year)]

                qcdparams = np.array([rl.IndependentParameter('qcdparam_'+cat+'_ptbin%d_regbin%d' % (ptbin, i), 0) for i in range(reg.nbins)])
                initial_qcd = failCh.getObservation()[0].astype(float)  # was integer, and numpy complained about subtracting float from it

                for sample in failCh:
                    initial_qcd -= sample.getExpectation(nominal=True)

                if np.any(initial_qcd < 0.):
                    raise ValueError('initial_qcd negative for some bins..', initial_qcd)

                sigmascale = 10  # to scale the deviation from initial                      
                scaledparams = initial_qcd * (1 + sigmascale/np.maximum(1., np.sqrt(initial_qcd)))**qcdparams
                fail_qcd = rl.ParametericSample('ptbin%d%sfail%s_qcd' % (ptbin, cat, year), rl.Sample.BACKGROUND, reg, scaledparams)
                failCh.addSample(fail_qcd)
                pass_qcd = rl.TransferFactorSample('ptbin%d%spass%s_qcd' % (ptbin, cat, year), rl.Sample.BACKGROUND, tf_params[cat][ptbin, :], fail_qcd)
                passCh.addSample(pass_qcd)


    df.to_csv('output/systematics.csv')

    with open(os.path.join(str(tmpdir), 'testModel_'+year+'.pkl'), 'wb') as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), 'testModel_'+year))

if __name__ == '__main__':

    print("Starting to run...")

    for directory in ['output', 'plots']:
        if not os.path.exists(directory):
            os.mkdir(directory)

    example_rhalphabet('output')

