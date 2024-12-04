#!/usr/bin/env python3
from adecoi import __version__

import adecoi.utils.custom_logger as custom_logger
from adecoi.utils.log_colours import green,cyan,red

import os
import sys
import argparse
import snakemake

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):
 
    parser = argparse.ArgumentParser(add_help=True)

    i_group = parser.add_argument_group('Run options')
    i_group.add_argument('-c',"--config", action="store",help="Input config file in yaml format, all command line arguments can be passed via the config file.", dest="config")
    i_group.add_argument('--reference_name', action="store",
                         help="Optional input for the reference name in the focal alignment. Default: Reference",
                         dest="reference_name")
    i_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads. Default: 1")
    i_group.add_argument("-h","--help",action="store_true",dest="help")

    """
    Exit with help menu if no args supplied
    """

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    config = setup_config_dict(cwd,args.config)
    
    # ready to run? either verbose snakemake or quiet mode

    status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,
                                    config=config, cores=config[KEY_THREADS],lock=False,quiet=True,log_handler=config[KEY_LOG_API]
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1



if __name__ == '__main__':
    main()
