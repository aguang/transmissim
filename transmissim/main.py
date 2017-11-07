import yaml
from transmissim import simulate

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params

    with open(params, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    for section in cfg:
        print(section)

    print(cfg['seed'])
    print(cfg['modules'])
    seed = random.randint(1,4294967295)
    if cfg['main']['seed']:
        print(cfg['main']['seed'])
        seed = int(cfg['main']['seed'])

    sim_network = cfg['modules']['network']
    sim_phylogeny = cfg['modules']['phylogeny']
    sim_sequence = cfg['modules']['sequence']

    # assert sim_network, sim_phylogeny, sim_sequence are booleans

    viral_tree_program = cfg['programs']['viraltreeprogram']
    reads = cfg['programs']['readsprogram']

    phylogeny_out = cfg['output']['phylogenyout']
    sequence_out = cfg['output']['sequenceout']

    sim_contact = cfg['network']['contact']
    sim_transmission_network = cfg['network']['transmission']
    R0 = cfg['network']['R0']
    number_of_hosts = cfg['network']['numberofhosts']
    duration = cfg['network']['duration']
    rate_import_case = cfg['network']['rateimportcase']

    if(sim_network && sim_phylogeny && sim_sequence):
        print("Module option: Simulate network, phylogeny, and sequences.")
        if(sim_contact != 0):
            print("Contact network simulation not implemented yet.")
        else:
            network = simulate.outbreaker(cluster_R0 = R0, n_hosts = number_of_hosts,
                cluster_duration = duration, rate_import_case = rate_import_case)

    elif(sim_network && sim_phylogeny && !sim_sequence):
        # run just network & phylogeny
        print("Not implemented yet.")

    elif(sim_network && !sim_phylogeny && sim_sequence):
        print("Network to Sequence simulation not implemented yet.")

    elif(!sim_network && sim_phylogeny && sim_sequence):
        # just simulate transmission tree & sequences
        print("Not implemented yet.")

    elif(!sim_network && !sim_phylogeny && sim_sequence):
        # just simulate genome sequences
        print("Not implemented yet.")

    elif(!sim_network && sim_phylogeny && !sim_sequence):
        # just simulate phylogeny
        print("Not implemented yet.")

    elif(sim_network && !sim_phylogeny && !sim_sequence):
        # just simulate network
        print("Not implemented yet.")

    else:
        print("All modules 0, nothing to simulate.")

        print(analysis_start == "Pipeline: all")
        if analysis_start == "all":
            network = outbreaker(cluster_R0, n_hosts, cluster_duration, rate_import_case)
            vt = transmission_tree(network, cluster_duration, ancestral_duration, tree_out, simphy_path, seed, birth_rate, death_rate)
            if analysis_end != "TT":
                t = vt.write(format=5)
                sequence(t, root_file)
                print("sequence finished")
                if analysis_end != "GS":
                    reads(art, sequencing_system,reads_out,read_length,coverage,mean_fragment_length,sd_fragment_length)

        if analysis_start == "TT":
            print("Pipeline: Transmission Tree -> Reads")
            t=''
            with open(full_tree, 'r') as f:
                t=f.readline()
            sequence(t, root_file)
            if analysis_end != "GS":
                reads(art, sequencing_system, reads_out, read_length, coverage,mean_fragment_length,sd_fragment_length)

        if analysis_start == "GS":
            reads(art, sequencing_system, reads_out, read_length, coverage,mean_fragment_length,sd_fragment_length)