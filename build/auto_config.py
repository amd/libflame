"""
Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All Rights Reserved.
Name: auto_config.py
Purpose: To check and recognize the zen architecture family
"""

#Global Imports
import platform
import re
import subprocess
import sys


def config_check():                                                             #Function to check the system info

    try:
        global model, family, vendor, stepping                                  #Global Variables
        if 'Windows' in platform.system():
            result = subprocess.Popen('wmic cpu get caption', shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()     #Execute wmic shell command with subprocess
            result = result[0].decode('utf-8')

            result = result.replace('\n', '')                                   #Replace the newline character with empty char
            parse_string = result.split(" ")                                    #Parse the string into list of string
            parse_string = [data for data in parse_string if data.strip()]      #Strip the empty strings from list

            vendor = parse_string[1]
            family = int(parse_string[3])
            model = int(parse_string[5])
            stepping = int(parse_string[7])
        elif 'Linux' in platform.system():
            result = subprocess.Popen('lscpu', shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()     #Execute lscpu command with subprocess

            stepping = int(re.findall(r'\WStepping:.*', result[0].
                                      decode('utf-8'), re.MULTILINE)[0].
                           strip('\n').split(' ')[-1])
            family = int(re.findall(r'\WCPU family:.*', result[0].
                                    decode('utf-8'), re.MULTILINE)[0].
                         strip('\n').split(' ')[-1])
            model = int(re.findall(r'\WModel:.*', result[0].decode('utf-8'),
                                   re.MULTILINE)[0].strip('\n').split(' ')[-1])
            vendor = re.findall(r'\WModel name:.*', result[0].decode('utf-8'),
                                re.MULTILINE)[0].strip('\n').split('Model '
                                                                   'name:')[-1]

        """
        AMD family numbers
        Zen/Zen+/Zen2 (0x17), Zen3/Zen4 (0x19) and zen5 (0x1a) family numbers
        """
        zen_family = [23, 25, 26]

        """
        Bulldozer / Piledriver / Steamroller / Excavator family number
        """
        amd_family = 21

        """
        AMD CPUID model numbers
        """
        zen_model = [48, 255]
        zen2_model = [0, 255]
        zen3_model = [(0, 15), (32, 95)]
        zen4_model = [(16, 31), (96, 175)]
        zen5_model = [(0, 15), (16, 31)]
        excavator_model = [96, 127]
        steamroller_model = [48, 63]
        piledriver_model = [2, 16, 31]
        bulldozer_model = [0, 1]

        if vendor.count("Intel64"):                                             #Check the CPU configuration Intel64/AMD64
            return
        elif 'AMD' in vendor:                                                   #.count("AMD64"):
            """
            Extracting AMD family name
            """
            if family == zen_family[0]:
                if zen_model[0] <= model <= zen_model[1]:
                    family="zen2"
                elif zen2_model[0] <= model <= zen2_model[1]:
                    family="zen"
                else:
                    print("Unknown model number")
            elif family == zen_family[1]:
                if (zen3_model[0][0] <= model <= zen3_model[0][1]) or (
                        zen3_model[1][0] <= model <= zen3_model[1][1]):
                    family="zen3"
                elif (zen4_model[0][0] <= model <= zen4_model[0][1]) or (
                        zen4_model[1][0] <= model <= zen4_model[1][1]):
                    family="zen4"
                else:
                    print("Unknown model number zen4")
            elif family == zen_family[2]:
                if ((zen5_model[0][0]) <= model <= (zen5_model[0][1])) or \
                    ((zen5_model[1][0]) <= model <= (zen5_model[1][1])):
                    family="zen5"
                else:
                    print("Unknown model number for zen5")
            elif family == amd_family:
                if excavator_model[0] <= model <= excavator_model[1]:           #Check for specific models of excavator family
                    family="excavator"
                elif steamroller_model[0] <= model <= steamroller_model[1]:     #Check for specific models of steamroller family
                    family="steamroller"
                elif model == piledriver_model[0] or \
                        (piledriver_model[1] <= model <= piledriver_model[2]):  #Check for specific models of piledriver family
                    family="piledriver"
                elif model == bulldozer_model[0] or \
                        model == bulldozer_model[1]:
                    family="bulldozer"
                else:
                    print("Unknown model number")
            else:
                print("Unknown family")
        else:
            print("UNKNOWN CPU")
        return family
    except Exception as e:
        print("Exception due to %s" % e)

config = config_check()                                                         #Function call for config family names
print(config)
