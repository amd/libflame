#!/usr/bin/env python3
"""
Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
Name: fla_format_code.py
Purpose: Formatting of .c, .cpp, .h and .hh files as per standards.
Script usage:
    python/python3 scripts/fla_format_code.py \
        <specific file/folder or list of multiple files/folders>
Examples:
    python3 scripts/fla_format_code.py test/main/src
    python3 scripts/fla_format_code.py test/main/src src/map
    python3 scripts/fla_format_code.py test/main/src/test_getrf.c
    python3 scripts/fla_format_code.py test/main/src/test_getrf.c src/map/lapack2flamec/FLA_gesdd.c
"""

# Global Python Imports
import os
import sys


def fla_format_code():
    
    """
    Formatting of .c, .cpp, .h and .hh files as per standards.
    """
    
    try:
        input_file_or_folder = sys.argv
        if isinstance(input_file_or_folder, list):
            for element in input_file_or_folder:
                if "fla_format_code.py" in element:
                    pass
                else:
                    if os.path.isfile(element):
                        file_path = element
                        new_file_path = element + "_new"
                        cmd = "clang-format -style=file %s > %s" \
                              % (file_path, new_file_path)
                        os.system(cmd)
                        os.remove(file_path)
                        os.rename(new_file_path, file_path)
                    elif os.path.isdir(element):
                        for root, dirs, files in os.walk(element):
                            for file in files:
                                if file.endswith(".c") or \
                                   file.endswith(".cpp") or \
                                   file.endswith(".h") or \
                                   file.endswith(".hh"):
                                    file_path = \
                                        os.path.abspath(os.path.join(root,
                                                                     file))
                                    new_file_path = os.path.\
                                        abspath(os.path.join(root,
                                                             file + "_new"))
                                    cmd = "clang-format -style=file %s > %s" \
                                          % (file_path, new_file_path)
                                    os.system(cmd)
                                    os.remove(file_path)
                                    os.rename(new_file_path, file_path)
    except Exception as e:
        print(e)


if __name__ == '__main__':
    fla_format_code()
