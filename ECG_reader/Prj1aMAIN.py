#On my honor,as an Aggie, I have neither given nor received unauthorized aid on this academic work. An Aggie does not
#lie, cheat or steal, or tolerate those who do

import matplotlib.pyplot as plt
# importing functions from other files
from ECG_Reader import read_ecg
from HeartRate import findHeartRate
from Diseases import *
from Intervals import *

# ask user for filename
filename = input("What is the filename: ")

# calls on read_ecg function from ECG_reader file
# function returns the header, ecg, and time values from the file
header, voltage, time = read_ecg(filename)

# calls on findHeartRate function from HeartRate file
# function calculates heart rate by first searching for R waves
# the function then divides 60 by every RR interval, and finding the average of those values
bpm = findHeartRate(voltage,time)

# all the interval functions from Intervals files are called on
# functions use findpeaks to find each wave and then calculate the time between them
QRS_interval = QRS_interval(voltage,time)
PPinterval = PP_interval(voltage,time)
RRinterval = RR_interval(voltage,time)
PRinterval = PR_interval(voltage,time)

# all functions from the Diseases file are called on
# function uses the required values to come with a diagnosis
# if one of the diagnosis is true, its name is stored in the variable diagnosis
if bundle_branch(QRS_interval,voltage) == True:
    diagnosis = "Bundle Branch Block"
elif arrythmia (PPinterval,RRinterval) == True:
    diagnosis = "Arrythmia"
elif atrial_fibrillation(voltage)== True:
    diagnosis = "Atrial Fibrillation"
elif AV_firstdegree(PRinterval)== True:
    diagnosis = "AV Block First Degree"
elif AV_seconddegree(PPinterval,RRinterval)== True:
    diagnosis = "AV Block Second Degree"
elif sinus_bradycardia (bpm) == True:
    diagnosis = "Sinus Bradycardia"
elif sinus_tachycardia(bpm) == True:
    diagnosis = "Sinus Tachycardia"
else:
    diagnosis = "Normal"

# time and voltage are graphed to create the ecg trace with the diagnosis in the title
# patient number and heart rate are printed
plt.plot(time,voltage)
plt.title(diagnosis)
plt.xlabel('Time (s)')
plt.ylabel('Voltage (mV)')
plt.title(diagnosis)
plt.show()
print("Heart Rate:",bpm,"bpm")
# patient number is found by looking for the index after the first number in the header
patient = header.index("Number")
print("Patient Number:", header[patient+1])


