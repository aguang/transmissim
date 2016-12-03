#!/usr/bin/env python
from rpy2.robjects.packages import importr
from transmission import binary_tree

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params
    if params:
        options = [line.strip().split('//')[0].strip() for line in open(params)]
        
        analysis = options[0]

        if analysis == "all" or "TN->R":

            outbreaker = importr('outbreaker')

            base = importr('base')
            w = base.rep(0.8, 350)
            duration = 350
            test = outbreaker.simOutbreak(R0 = 2, infec_curve=w, n_hosts=200, duration=duration, rate_import_case=0)
            full_tree = binary_tree(test)
            print(full_tree)

            print("done")
