#!/usr/bin/env python3
"""
Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
Name: fla_format_code.py
Purpose: Formatting of .c, .cpp, .h and .hh files as per standards.
Script usage: 
    python/python3 scripts/fla_format_code.py <path to folder or specific file>
Examples:
    python scripts/fla_format_code.py test/main/src
    python scripts/fla_format_code.py test/main/src/test_getrf.c
"""

# Global Python Imports
import os
import sys


def fla_format_code():
    
    """
    Formatting of .c, .cpp, .h and .hh files as per standards.
    """
    
    try:
        input_file_or_folder = sys.argv[1]

        if os.path.isfile(input_file_or_folder):
            file_path = input_file_or_folder
            new_file_path = input_file_or_folder + "_new"
            cmd = "clang-format -style=file %s > %s" \
                  % (file_path, new_file_path)
            os.system(cmd)
            os.remove(file_path)
            os.rename(new_file_path, file_path)
        elif os.path.isdir(input_file_or_folder):
            for root, dirs, files in os.walk(input_file_or_folder):
                for file in files:
                    if file.endswith(".c") or file.endswith(".cpp") or \
                       file.endswith(".h") or file.endswith(".hh"):
                        file_path = os.path.abspath(os.path.join(root, file))
                        new_file_path = \
                            os.path.abspath(os.path.join(root, file + "_new"))
                        cmd = "clang-format -style=file %s > %s" \
                              % (file_path, new_file_path)
                        os.system(cmd)
                        os.remove(file_path)
                        os.rename(new_file_path, file_path)
    except Exception as e:
        print(e)


if __name__ == '__main__':
    fla_format_code()

