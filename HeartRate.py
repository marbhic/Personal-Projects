#On my honor,as an Aggie, I have neither given nor received unauthorized aid on this academic work. An Aggie does not
#lie, cheat or steal, or tolerate those who do


from scipy.signal import find_peaks
# the purpose of this function is to calculate heart rate
def findHeartRate (voltage,time):
    # the R waves locations are first found using findpeaks
    peaks, _ = find_peaks(voltage, height=.600)

    # the locations are used to find the corresponding ecg and time values and are then appended
    # to their appropriate lists
    peak_val = []
    volt_val = []
    x = 0
    while x != len(peaks):
        r = peaks[x]
        y = time[r]
        z = voltage[r]
        peak_val.append(y)
        volt_val.append(z)
        x += 1

    # the RR interval is calculated by finding the time between each R wave
    # to find the heart rate, 60 is divided by the RR interval
    x = 1
    time_diff = []
    while x != len(peak_val):
        diff = peak_val[x] - peak_val[x - 1]
        bpm = 60 / diff
        x += 1
        time_diff.append(bpm)
    # the average of all the bpm values is used to determine the heart rate
    bpm = (sum(time_diff) / len(time_diff))
    bpm = round(bpm,2)
    return bpm