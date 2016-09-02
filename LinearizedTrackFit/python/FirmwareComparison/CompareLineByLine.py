import sys# , getopt
import argparse
from ConvertFirmwareOutputToInt import *

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--firmware', required=True, help="Firmware output file", dest="firmware_output_file_name")
    parser.add_argument('-e', '--emulator', required=True, help="Emulator output file", dest="emulator_output_file_name")
    args = parser.parse_args()

    # print 'Firmware output file is "', args.firmware_output_file_name
    # print 'Emulator output file is "', args.emulator_output_file_name

    firmware_output_file_name_int = convert_firmare_output_to_int(args.firmware_output_file_name)

    firmware_lines = open(firmware_output_file_name_int).readlines()
    emulator_lines = open(args.emulator_output_file_name).readlines()

    total_tracks = 0
    for line in range(len(emulator_lines)):
        if firmware_lines[line] != emulator_lines[line]:
            print "Found a difference at line "+str(total_tracks)+":"
            print "emulator =", emulator_lines[line]
            print "firmware =", firmware_lines[line]
            return -1
        total_tracks += 1

    print
    print "Compared "+str(total_tracks)+" tracks."
    print "Firmware and emulator match exactly."


if __name__ == "__main__":
    main(sys.argv[1:])
