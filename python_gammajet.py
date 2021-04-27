# generating the CSV
def getIsoCutValue(centrality):
    if centrality < 10:
        return 6.5
    elif centrality < 30:
        return 5.0
    elif centrality < 50:
        return 3.5
    else:
        return 1.5


# this only considers the highest pT cluster in each event
def main(ntuplefilenames, csvfilename):
    clusterptrange = (15, 40)
    jetptmin = 5

    start = time.time()
    with open(csvfilename, 'wb') as csvfile:
        fieldnames = ['centrality_v0m', 'icluster',
                      'cluster_pt', 'cluster_Lambda', 'cluster_phi', 'cluster_eta', 'cluster_iso_tpc_02_sub',
                      'jet_ak04tpc_pt_raw', 'jet_ak04tpc_phi', 'jet_ak04tpc_eta']
        csvwriter = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        csvwriter.writeheader()

        # need to normalize by number of clusters, so keep track
        # it'll get 1 added to it before the rows start being written
        # so may as well zero index
        clustercounter = -1

        for (ifile, ntuplefilename) in enumerate(ntuplefilenames):
            print('{0} Processing {1} ({2}/{3})'.format(datetime.datetime.now(), os.path.basename(ntuplefilename), ifile + 1, len(ntuplefilenames)))
            rootfile = ROOT.TFile.Open(ntuplefilename, 'READ')
            tree = rootfile.Get('AliAnalysisTaskNTGJ/_tree_event')
            nevents = tree.GetEntries()

            # get relevant branches
            acluster_pt = rnp.tree2array(tree, branches='cluster_pt')
            acluster_eta = rnp.tree2array(tree, branches='cluster_eta')
            acluster_phi = rnp.tree2array(tree, branches='cluster_phi')
            acluster_e_cross = rnp.tree2array(tree, branches='cluster_e_cross')
            acluster_e_max = rnp.tree2array(tree, branches='cluster_e_max')
            acluster_ncell = rnp.tree2array(tree, branches='cluster_ncell')
            acluster_tof = rnp.tree2array(tree, branches='cluster_tof')
            acluster_distance_to_bad_channel = rnp.tree2array(tree, branches='cluster_distance_to_bad_channel')
            acluster_lambda_square = rnp.tree2array(tree, branches='cluster_lambda_square')
            acluster_iso_tpc_02 = rnp.tree2array(tree, branches='cluster_iso_tpc_02')
            acluster_iso_tpc_02_ue = rnp.tree2array(tree, branches='cluster_iso_tpc_02_ue')
            aue_estimate_tpc_const = rnp.tree2array(tree, branches='ue_estimate_tpc_const')

            ajet_ak04tpc_pt_raw = rnp.tree2array(tree, branches='jet_ak04tpc_pt_raw')
            ajet_ak04tpc_phi = rnp.tree2array(tree, branches='jet_ak04tpc_phi')
            ajet_ak04tpc_eta = rnp.tree2array(tree, branches='jet_ak04tpc_eta')

            acentrality_v0m = rnp.tree2array(tree, branches='centrality_v0m')

            # loop through events
            for ievent in range(nevents):
                centrality_v0m = acentrality_v0m[ievent]
                ue_estimate_tpc_const = aue_estimate_tpc_const[ievent]

                clusterptmax = 0
                iclusterptmax = -1

                # skip events outside of the centrality range of interest
                if centrality_v0m < 0 or centrality_v0m > 90:
                    continue

                for icluster in range(len(acluster_pt[ievent])):
                    # first, check if the cluster passes all the cluster cuts
                    if abs(acluster_eta[ievent][icluster]) > 0.67:
                        continue
                    if (acluster_e_cross[ievent][icluster] / acluster_e_max[ievent][icluster]) < 0.05:
                        continue
                    if acluster_ncell[ievent][icluster] < 2:
                        continue
                    if acluster_distance_to_bad_channel[ievent][icluster] < 1.0:
                        continue
                    if abs(acluster_tof[ievent][icluster]) > 20:
                        continue

                    # check if the cluster passes the isolation cut
                    cluster_iso_tpc_02_sub = acluster_iso_tpc_02[ievent][icluster] + acluster_iso_tpc_02_ue[ievent][icluster] - (ue_estimate_tpc_const * 0.2 * 0.2 * np.pi)
                    if cluster_iso_tpc_02_sub > getIsoCutValue(centrality_v0m):
                        continue

                    # check if the cluster is the highest-pt cluster yet found
                    cluster_pt = acluster_pt[ievent][icluster]
                    if cluster_pt > clusterptmax:
                        clusterptmax = cluster_pt
                        iclusterptmax = icluster

                # if the max-pt cluster is out of the ptrange, skip this event
                if clusterptmax < clusterptrange[0] or clusterptmax > clusterptrange[1]:
                    continue

                cluster_lambda0 = acluster_lambda_square[ievent][iclusterptmax][0]
                cluster_phi = acluster_phi[ievent][iclusterptmax]
                cluster_eta = acluster_eta[ievent][iclusterptmax]
                cluster_iso = acluster_iso_tpc_02[ievent][iclusterptmax] + acluster_iso_tpc_02[ievent][iclusterptmax] - (ue_estimate_tpc_const * 0.2 * 0.2 * np.pi)
                clustercounter = clustercounter + 1

                # loop through jets
                for ijet in range(len(ajet_ak04tpc_pt_raw[ievent])):
                    # ignore jets below pT threshold
                    jet_pt = ajet_ak04tpc_pt_raw[ievent][ijet]
                    if jet_pt < jetptmin:
                        continue

                    jet_phi = ajet_ak04tpc_phi[ievent][ijet]
                    jet_eta = ajet_ak04tpc_eta[ievent][ijet]

                    # fill csv
                    row = {}
                    row['centrality_v0m'] = centrality_v0m
                    row['icluster'] = clustercounter
                    row['cluster_pt'] = clusterptmax
                    row['cluster_Lambda'] = cluster_lambda0
                    row['cluster_phi'] = cluster_phi
                    row['cluster_eta'] = cluster_eta
                    row['cluster_iso_tpc_02_sub'] = cluster_iso
                    row['jet_ak04tpc_pt_raw'] = jet_pt
                    row['jet_ak04tpc_phi'] = jet_phi
                    row['jet_ak04tpc_eta'] = jet_eta

                    csvwriter.writerow(row)

            rootfile.Close()


# CSV as input
def getdeltaphi(jet_phi, cluster_phi):
    dphi = jet_phi - cluster_phi
    if dphi < 0:
        dphi = dphi + 2 * np.pi
    if dphi > np.pi:
        dphi = 2 * np.pi - dphi

    return dphi

def getDfFromCsv(filename):
    start = time.time()
    df = pd.read_csv(filename, '\t')
    df['delta_phi'] = df.apply(lambda x: getdeltaphi(x['jet_ak04tpc_phi'], x['cluster_phi']), axis=1)
    df['ptratio'] = df['jet_ak04tpc_pt_raw'] / df['cluster_pt']
    df['weights'] = 1
    df = df.query('abs(jet_ak04tpc_eta) < 0.5')
    end = time.time()
    print('Processed {0} in {1:.0f} seconds'.format(os.path.basename(filename), end - start))
    return df


def purities15(centrality, pt):
    if 0 < centrality and centrality < 10:
        if 15 < pt and pt < 20: return 0.306
        if 20 < pt and pt < 25: return 0.378
        if 25 < pt and pt < 40: return 0.253
    if 10 < centrality and centrality < 30:
        if 15 < pt and pt < 20: return 0.350
        if 20 < pt and pt < 25: return 0.285
        if 25 < pt and pt < 40: return 0.477
    if 30 < centrality and centrality < 50:
        if 15 < pt and pt < 20: return 0.317
        if 20 < pt and pt < 25: return 0.339
        if 25 < pt and pt < 40: return 0.508
    if 50 < centrality and centrality < 90:
        if 15 < pt and pt < 20: return 0.249
        if 20 < pt and pt < 25: return 0.333
        if 25 < pt and pt < 40: return 0.539
    return 1


def purities18(centrality, pt):
    if 0 < centrality and centrality < 10:
        if 15 < pt and pt < 20: return 0.426
        if 20 < pt and pt < 25: return 0.404
        if 25 < pt and pt < 40: return 0.382
    if 10 < centrality and centrality < 30:
        if 15 < pt and pt < 20: return 0.350
        if 20 < pt and pt < 25: return 0.337
        if 25 < pt and pt < 40: return 0.320
    if 30 < centrality and centrality < 50:
        if 15 < pt and pt < 20: return 0.336
        if 20 < pt and pt < 25: return 0.444
        if 25 < pt and pt < 40: return 0.505
    if 50 < centrality and centrality < 90:
        if 15 < pt and pt < 20: return 0.401
        if 20 < pt and pt < 25: return 0.441
        if 25 < pt and pt < 40: return 0.517
    return 1


def getPurityWeight(centrality, pt, sigma2long, runperiod):
    if runperiod == '2015':
        p = purities15(centrality, pt)
    elif runperiod == '2018':
        p = purities18(centrality, pt)

    if sigma2long < 0.3:
        return 1 / p
    elif 0.4 < sigma2long and sigma2long < 1.5:
        return (1 - p) / p
    else:
        return 1

# this is currently slow
df15['purityweights'] = df15.apply(lambda x: getPurityWeight(x['centrality_v0m'], x['cluster_pt'], x['cluster_Lambda'], '2015'), axis=1)    
df18['purityweights'] = df18.apply(lambda x: getPurityWeight(x['centrality_v0m'], x['cluster_pt'], x['cluster_Lambda'], '2018'), axis=1)
