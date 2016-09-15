#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
# Script meant to be run from pipeline repo, append that to path:
sys.path.append('./')


def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Awesome script.')
    parser.add_argument(
        '-i', metavar='input', type=str, help="Input.")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_arguments()
