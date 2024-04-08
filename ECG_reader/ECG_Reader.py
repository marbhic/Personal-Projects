#On my honor,as an Aggie, I have neither given nor received unauthorized aid on this academic work. An Aggie does not
#lie, cheat or steal, or tolerate those who do


import numpy as np
# the purpose of this function is to read the file and extract the header, ecg values, and time values
def read_ecg(filename):
    # opens the file based on what is stored in filename variable
    variables = open(filename,'r')

    # the file is read line by line
    # the header is the first line and is stored as its own variable
    data = variables.readlines()
    header = data[0]

    # the other lines in the files are the ecg values and are appended to the voltage list
    voltage = []
    x = 1
    while x != len(data):
        y = data[x]
        voltage.append(float(y))
        x +=1

    # header is split to find the time step
    header = header.split()
    index = header.index("Step")
    time_step = float(header[index+1])

    # the stop time for the ecg is determined by multiplying the number of ecg values by the time step
    stop = float(len(voltage) * time_step)
    # the time values are found by creating an array with a range from 0 to the calculated stop time with
    # the time step serving as the time gaps
    time = np.arange(0,stop,time_step)
    variables.close()
    # function returns header, ecg values, and time values when called on
    return header, voltage, time