import yaml
import simulate

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='simulation')

    ds = ' [%(default)s]'
    parser.add_argument('-p', '--params', help='parameters file for simulation')
    opts = parser.parse_args()

    params = opts.params

    with open(params, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    simulate.main(cfg)