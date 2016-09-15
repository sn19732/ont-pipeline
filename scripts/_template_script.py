#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
# Script meant to be run from pipeline repo, append lib to path:
sys.path.append('lib/')

def parse_arguments():
    """
    Parse arguments.
    """
    parser = argparse.ArgumentParser(
        description='Awesome script.')
    parser.add_argument(
        '-i', metavar='input', type=str, help="Input.", required=True)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
