#!/usr/bin/env python
'''
A script to retrieve the limits

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import glob
import pwd
import argparse
import errno
import socket
import signal
import logging
import math
import ROOT
import subprocess
from multiprocessing import Pool
from socket import gethostname

scratchDir = 'data' if 'uwlogin' in gethostname() else 'nfs_scratch'

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def limitsWrapper(args):
    getLimits(*args)

def runCommand(command):
    return subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

def getCommands(**kwargs):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    combineCommands = [
        'combine -M AsymptoticLimits -m {a} {datacard} -n "HToAAH{h}A{a}_{tag}"',
        'combine -M Significance -m {a} {datacard} -n "HToAAH{h}A{a}_{tag}Observed"',
        'combine -M Significance -m {a} {datacard} -n "HToAAH{h}A{a}_{tag}APriori" -t -1 --expectSignal=1',
        'combine -M Significance -m {a} {datacard} -n "HToAAH{h}A{a}_{tag}APosteriori" -t -1 --expectSignal=1 --toysFreq',
    ]

    return combineCommands


def runLimit(tag,h,a,**kwargs):
    dryrun = kwargs.get('dryrun',False)

    datacard = 'datacards_shape/MuMuTauTau/mmmt_{}_HToAAH{}A{}.txt'.format(tag,h,'X' if 'parametric' in tag else a)

    combineCommands = getCommands(tag,h,a,**kwargs)
    for cc in combineCommands:
        cc.format(datacard=datacard,h=h,a=a,tag=tag)
        if dryrun:
            logging.info(cc)
        else:
            out = runCommand(cc)
            print out

def submitLimit(tag,h,amasses,**kwargs):
    dryrun = kwargs.get('dryrun',False)
    jobName = kwargs.get('jobName',None)
    pointsPerJob = kwargs.get('pointsPerJob',10)

    a = '${A}'

    datacard = 'datacards_shape/MuMuTauTau/mmmt_{}_HToAAH{}A{}.txt'.format(tag,h,'X' if 'parametric' in tag else '${A}')

    combineCommands = getCommands(**kwargs)

    sample_dir = '/{}/{}/condor_projects/{}/{}/{}'.format(scratchDir,pwd.getpwuid(os.getuid())[0], jobName, tag, h)
    python_mkdir(sample_dir)

    # create submit dir
    submit_dir = '{}/submit'.format(sample_dir)
    if os.path.exists(submit_dir):
        logging.warning('Submission directory exists for {0}.'.format(jobName))
        return

    # make dag
    dag_dir = '{}/dags/'.format(sample_dir)
    python_mkdir(dag_dir)

    # create file list
    input_name = '{}/amasses.txt'.format(dag_dir)
    with open(input_name,'w') as f:
        for ai in amasses:
            f.write('{}\n'.format(ai))

    # output dir
    output_dir = 'srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/{}/{}/{}/{}'.format(pwd.getpwuid(os.getuid())[0], jobName, tag, h)

    # create bash script
    bash_name = '{}/script.sh'.format(sample_dir)
    bashScript = '#!/bin/bash\n'
    bashScript += 'ls\n'
    bashScript += 'printenv\n'
    bashScript += 'cp -r $CMSSW_VERSION/src/datacards_shape .\n'
    bashScript += 'while read A; do\n'
    for cc in combineCommands:
        bashScript += cc.format(datacard=datacard,h=h,a=a,tag=tag)+'\n'
    bashScript += 'done < $INPUT\n'
    with open(bash_name,'w') as file:
        file.write(bashScript)
    os.system('chmod +x {0}'.format(bash_name))

    # create farmout command
    farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --fwklite --job-generates-output-name'
    farmoutString += ' --input-file-list={} --input-files-per-job={} --assume-input-files-exist'.format(input_name,pointsPerJob)
    farmoutString += ' --submit-dir={} --output-dag-file={}/dag --output-dir={}'.format(submit_dir, dag_dir, output_dir)
    farmoutString += ' --extra-usercode-files="src/datacards_shape/MuMuTauTau" {} {}'.format(jobName, bash_name)

    if not dryrun:
        logging.info('Submitting {}/{}/{}'.format(jobName,tag,h))
        os.system(farmoutString)
    else:
        print farmoutString





def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Process limits')

    # limit information
    parser.add_argument('tag', nargs='?',type=str,default='',help='MuMuTauTau tag')
    parser.add_argument('--mh', nargs='?',type=int,default=125,help='Higgs mass')
    parser.add_argument('--ma', nargs='?',type=int,default=15,help='Pseudoscalar mass')
    # job submission
    parser.add_argument('--jobName', nargs='?',type=str,default='',help='Jobname for submission')
    parser.add_argument('--submit',action='store_true',help='Submit Full CLs')
    parser.add_argument('--dryrun',action='store_true',help='Dryrun for submission')
    # logging
    parser.add_argument('-j',type=int,default=1,help='Number of cores')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')

    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    if args.submit:
        amasses = range(5,23,2)
        if 'parametric' in args.tag: amasses = [x*0.1 for x in range(50,211,1)]
        submitLimit(args.tag,args.mh,amasses,dryrun=args.dryrun,jobName=args.jobName)
    else:
        runLimit(args.tag,args.mh,args.ma,dryrun=args.dryrun)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)

