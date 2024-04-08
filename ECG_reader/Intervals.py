#On my honor,as an Aggie, I have neither given nor received unauthorized aid on this academic work. An Aggie does not
#lie, cheat or steal, or tolerate those who do


from scipy.signal import find_peaks

# this function finds the time between P waves
def PP_interval (voltage,time):
    # findpeaks is used to find the loations of the P waves
    P_peaks, _ = find_peaks(voltage, height=(.25, .4))
    # the location of the P waves is used to find the time values of the P values
    Ptime_val = []
    for i in range(0,len(P_peaks)):
        y = P_peaks[i]
        x = time[y]
        Ptime_val.append(x)
    # the time values are used to find the time intervals between the P waves
    PPinterval = []
    for i in range(1,len(Ptime_val)):
        interval = Ptime_val[i] - Ptime_val[i-1]
        PPinterval.append(interval)
    return PPinterval

# this function finds the time between R waves
def RR_interval (voltage,time):
    # findpeaks is used to find the location of the R waves
    R_peaks, _ = find_peaks(voltage, height=1)
    # the locations are used to find the time values of the R waves
    Rtime_val = []
    for i in range(0,len(R_peaks)):
        y = R_peaks[i]
        x = time[y]
        Rtime_val.append(x)
    # the time values are used to find the time intervals between the R waves
    RRinterval = []
    for i in range(1,len(Rtime_val)):
        interval = Rtime_val[i] - Rtime_val[i-1]
        RRinterval.append(interval)
    return RRinterval

# this function finds the interval between P and R waves
def PR_interval (voltage,time):
    # findpeaks is used to find the locations of the P and R waves
    P_peaks, _ = find_peaks(voltage, height=(.25, .4))
    R_peaks, _ = find_peaks(voltage, height=1)
    # if there were no P or R waves, the PR interval can't be found
    if len(P_peaks)==0 or len(R_peaks) == 0:
        PRinterval = 0
    elif len(P_peaks) != len(R_peaks):
        PRinterval = 0

    else:
        # the locations of the P and R waves are used to find the time values of the corresponding waves
        Ptime_val = []
        Rtime_val = []
        for i in range(0,len(P_peaks)):
            y = P_peaks[i]
            x = time[y]
            v = R_peaks[i]
            z = time[v]
            Ptime_val.append(x)
            Rtime_val.append(z)
        PRinterval = []
        #the time values are use to find the PR interval
        for i in range(0,len(Ptime_val)):
            interval = Rtime_val[i] - Ptime_val[i]
            PRinterval.append(interval)
    return PRinterval

def QRS_interval(voltage,time):
    #the QRS interval is found by finding the width of the wave

    # the inverse of ecg is used to find the Q and S peaks
    voltage1 = []
    x = 1
    while x != len(voltage):
        y = voltage[x] * -1
        voltage1.append(float(y))
        x += 1
    # findpeaks is used to find the locations of the Q waves
    Q_peaks, _ = find_peaks(voltage1, height=(.100, .2))
    # the zero value before each Q wave is looked for
    Qvoltage_val = []
    Qzero_val = []
    for i in range(0, len(Q_peaks)):
        # the location of the Q wave
        zero_loc = Q_peaks[i]
        # while loop checks if zero_loc provides a zero value
        # if the value is zero the corresponding time and ecg value are appended
        # if the value is not zero, the value of zero_loc decreases by one
        y = 1
        while y != 0:
            r = voltage[zero_loc]
            if r == -0.0:
                z = time[zero_loc]
                Qzero_val.append(z)
                p = voltage[zero_loc]
                Qvoltage_val.append(p)
                y = 0
            elif r != 0:
                zero_loc = zero_loc - 1
                y += 1
            else:
                break

    # findpeaks is used to find the locations of the S waves
    S_peaks, _ = find_peaks(voltage1, height=.2)

    # the zero value before each Q wave is looked for
    Szero_val = []
    Svoltage_val = []
    for i in range(0, len(S_peaks)):
        # while loop checks if zero_loc provides a zero value

        # if the value is not zero, the value of zero_loc decreases by one
        t = S_peaks[i]
        y = 1
        while y != 0:
            r = voltage[t]
            # if the value is zero the corresponding time and ecg value are appended
            if r == -0.0:
                z = time[t]
                Szero_val.append(z)
                p = voltage[t]
                Svoltage_val.append(p)
                y = 0
            elif r != 0:
                t = t + 1
                y += 1
            else:
                break
    # the zero time values for the Q and S waves are used to find the QRS interval
    QRS_interval = []
    for i in range(0, len(Szero_val)):
        interval = abs(Qzero_val[i] - Szero_val[i])
        QRS_interval.append(interval)
    return QRS_interval