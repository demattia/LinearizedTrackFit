def twos_comp(val, bits):
    """compute the 2's compliment of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val                         # return positive value as is


def convert_firmare_output_to_int(firmware_output_file_name):
    latency = 43
    firmware_output_file_name_int = firmware_output_file_name.rstrip(".txt")+"_Int.txt"
    output_file = open(firmware_output_file_name_int, "w")
    for line in open(firmware_output_file_name):
        split_line = line.split(" ")
        clock = int(split_line[0])
        if clock < latency:
            continue
        output_line = ""
        for i in range(1, 13):
            if split_line[i].find("X") != -1:
                continue
            output_line += str(twos_comp(int(split_line[i], 2), len(split_line[i]))) + " "
        # print clock, output_line.rstrip(" ")
        output_file.write(output_line+"\n")

    output_file.close()
    return firmware_output_file_name_int
