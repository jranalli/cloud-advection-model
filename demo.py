import numpy as np
import matplotlib.pyplot as plt

import camsmoothing

'''
Demonstrates the use of the CAM based on a known uniform 1-d plant spatial 
distribution.
'''

# Read the data
filedata = np.genfromtxt('livermore.csv', delimiter=',')
time = (filedata[:, 0] - filedata[0, 0]) * 24 * 60 * 60  # Remove time offset, convert to seconds
data = filedata[:, 1]

# Extract the data time step
dt = time[1] - time[0]

# Define the plant. This represents a one-dimensional plant that is uniformly
# distributed over 5 km. The reference is placed at a position of 500 m from
# the leading edge of the plant. The cloud speed is 20 m/s.
dx = 1  # 1 meter plant spatial resolution
plant = np.zeros(int(50000/dx))
plant[0:int(5000/dx)] = 1
U = 20
refpos = 500
smoothdata, (freq, filt) = camsmoothing.cam(dt, data, dx, plant, U, ref_pos=refpos)

# Generate a plot of the results
plt.figure(figsize=[4,6])
plt.subplot(211)
plt.plot(time/60/60, data)
plt.plot(time/60/60, smoothdata)
plt.xlim([0,24])
plt.xticks(np.arange(0, 25, 3))
plt.title('Smoothed Irradiance')
plt.legend(['Ref. GHI', 'Smooth GHI'])
plt.ylabel('GHI')
plt.xlabel('Time (hr)')

plt.subplot(212)
plt.semilogx(freq[0:len(filt) // 2], np.abs(filt[0:len(filt) // 2]))
plt.title('Plant Transfer Function')
plt.ylabel('Transfer Function Magnitude')
plt.xlabel('Frequency (hz)')
plt.tight_layout()
plt.show()
