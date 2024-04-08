#On my honor,as an Aggie, I have neither given nor received unauthorized aid on this academic work. An Aggie does not
#lie, cheat or steal, or tolerate those who do


from scipy.signal import find_peaks

# this function diagnoses AV block first degree
def AV_firstdegree(PRinterval):
    # if there are no corresponding PR intervals for calculations, AV first degree is not the diagnosis
    if PRinterval == 0:
        value = False

    # AV first degree is diagnosed through a PR interval greater than .2 sec
    # the average of the PR intervals are used to determine if the interval exceeds .2 sec
    else:
        PRinterval = sum(PRinterval)/len(PRinterval)
        if PRinterval > .2:
            value = True
        else:
            value = False
    return value

# this function diagnoses AV block second degree
def AV_seconddegree(PPinterval,RRinterval):
    # AV block second degree is diagnosed when there are missing QRS complexes

    # if there are no P waves or R waves diagnosis should be false
    if len(PPinterval) == 0 or len(RRinterval) == 0:
        diagnosis = False

    # if there are QRS complexes missing, the number of P intervals should be considerably bigger
    elif (len(PPinterval) - len(RRinterval)) > 2:
        diagnosis = True
    else:
        diagnosis = False
    return diagnosis

# this function diagnoses tachycardia
def sinus_tachycardia (bpm):
    # tachycardia is diagnosed when the heart rate is over 100 bpm
    if bpm > 100:
        diagnosis = True
    else:
        diagnosis = False
    return diagnosis

# this function diagnoses bradycardia
def sinus_bradycardia (bpm):
    # bradycardia is diagnosed when the heart rate is less than 60 bpm
    if bpm < 60:
        diagnosis = True
    else:
        diagnosis = False
    return diagnosis

# this function diagnoses arrythmia
def arrythmia(PPinterval,RRinterval):
    # if any value of the RR interval is greater than 1, the diagnosis is false
    x = 0
    r = []
    # while loop looks for the values over 1 and append them to the r list
    while x != len(RRinterval):
        y = RRinterval[x]
        if y > 1:
            r.append(y)
            x += 1
        else:
            x += 1
    if len(r) > 0:
        diagnosis = False
    else:
        # if there are no PP intervals or RR interval, the diagnosis is false
        if len(PPinterval)==0 or len(RRinterval)==0:
            diagnosis = False
        else:
            # one way of diagnosing is finding the difference between the greatest PP interval and the smallest PP interval
            # the interval must be greater than .16
            max_interval = max(PPinterval)
            min_interval = min(PPinterval)
            PPinterval = max_interval - min_interval

            # another way of diagnosing is finding the difference between the greatest RR interval and the smallest RR interval
            # the interval must be greater than .16
            max_interval1 = max(RRinterval)
            min_interval1 = min(RRinterval)
            RRinterval = max_interval1 - min_interval1

            #if any of the two intervals are greater than .16, the diagnosis is true
            if PPinterval > .16 or RRinterval > .16:
                diagnosis = True
            else:
                diagnosis = False
    return diagnosis

# this function diagnoses bundle branch block
def bundle_branch(QRS_interval,voltage):
    # to diagnose bundle branch block, the length of the QRS complex must be greater than .12 sec
    # the average of the QRS complexes is used to determine if it is greater than .12 sec
    QRS_interval = sum(QRS_interval)/len(QRS_interval)
    P_peaks, _ = find_peaks(voltage, height=(.25, .4))

    # if there are no P waves, diagnosis is false
    if len(P_peaks) == 0:
        value = False
    #if the average is greater than .12, the diagnosis is true
    elif QRS_interval>.120:
        value = True
    else:
        value = False
    return value

def atrial_fibrillation(voltage):
    # findpeaks is used to first find the P waves
    P_peaks, _ = find_peaks(voltage, height=(.25, .4))
    # when diagnosing atrial fibrillation, there are usually no P waves
    if len(P_peaks) == 0:
        diagnosis = True
    else:
        diagnosis = False
    return diagnosis