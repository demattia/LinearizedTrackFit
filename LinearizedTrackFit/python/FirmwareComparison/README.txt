# Scripts to Compare Emulator and Firmware

The output of the firmware can be obtained by running the test bench on the input file produced by the emulator. The two results can be compared by running the script "CompareLineByLine.py". This script will first convert the output of the firmware (in binary format) to signed integer and then compare it line by line to the output produced by the emulator. The first 43 lines of the firmware output are skipped to account for the latency. The lines contain all track parameters and chi2 terms. If a difference is found the script stops and prints the line number and the contents from both files. If no difference is found the script prints the total number of compared lines (tracks) and a statement for an exact match.

For specific usage instructions, run the script with the -h or --help option.
