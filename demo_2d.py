import numpy as np
import matplotlib.pyplot as plt

import camsmoothing

'''
Demonstrates the use of the CAM by projecting a 2-d plant point cloud onto a
1-d plant spatial distribution. In this case we will treat each individual
point as an infinitesimally small generator. Each will be projected onto the
line passing through the reference in the cloud motion direction. 

Other methods of representing a 2-d plant in 1-d could be envisioned and are
an area of future research.
'''


def main():
    # Read the data
    filedata = np.genfromtxt('livermore.csv', delimiter=',')
    time = (filedata[:, 0]-filedata[0, 0])*24*60*60  # Remove time offset, convert to seconds
    data = filedata[:, 1]

    # Extract the data time step
    dt = time[1]-time[0]

    # Generate a 500 random point cloud over a 5km x 5km area to represent the
    # plant.
    np.random.seed(42)  # Set a fixed seed
    east = np.random.randint(0, 5000, 500)
    north = np.random.randint(0, 5000, 500)
    refxy = [np.random.randint(0, 5000), np.random.randint(0, 5000)]
    wind_vec = [1, 1]
    projpts = project_points(zip(east, north), refxy, wind_vec)

    # Shift so the plant distances so that the plant begins at x = zero
    refpos = -np.min(projpts)
    projpts = np.array(projpts)-np.min(projpts)

    # Generate the plant by placing a 1 at every position along the projected
    # line where a measurement is located.
    dx = 1
    plant = np.zeros(int(np.max(projpts)/dx))
    plant = np.pad(plant, (0, int(30000/dx)))
    inds = np.round(projpts / dx).astype(np.int32)
    plant[inds] = 1  # Uniform generation for each point
    U = 20  # Wind speed

    # Compute the smoothed time series
    smoothdata, (freq, filt) = camsmoothing.cam(dt, data, dx, plant, U, ref_pos=refpos)

    # Generate figures
    plt.figure(figsize=[12, 8])
    plt.subplot(221)
    plt.scatter(east, north)
    plt.scatter(refxy[0], refxy[1])
    plt.arrow(refxy[0], refxy[1], 500 * wind_vec[0], 500 * wind_vec[1],
            length_includes_head=True, width=25, head_width=100,
            color="#000000")
    plt.xlabel("Easting (m)")
    plt.ylabel("Northing (m)")
    plt.title("2D Plant Layout")
    plt.legend(["Wind Dir", "Plant", "Reference"])
    plt.gca().set_aspect('equal')

    plt.subplot(222)
    plt.plot(plant)
    plt.xlim([-200, np.max(projpts)*1.05])
    plt.xlabel("Projected Distance (m)")
    plt.ylabel("Plant Generation Density (m)")
    plt.title("Projected Plant Distribution")

    plt.subplot(223)
    plt.plot(time/60/60, data)
    plt.plot(time/60/60, smoothdata)
    plt.xlim([0,24])
    plt.xticks(np.arange(0, 25, 3))
    plt.title('Smoothed Irradiance')
    plt.legend(['Ref. GHI', 'Smooth GHI'])
    plt.ylabel('GHI')
    plt.xlabel('Time (hr)')

    plt.subplot(224)
    plt.semilogx(freq[0:len(filt) // 2], np.abs(filt[0:len(filt) // 2]))
    plt.title('Plant Transfer Function')
    plt.ylabel('Transfer Function Magnitude')
    plt.xlabel('Frequency (hz)')
    plt.tight_layout()
    plt.show()


# Helper function for projecting all the points
def project_points(pos, refpos, wind_dir):
    """
    Project a set of points onto a line passing through a reference in a given
    direction.

    :param pos: A list of [E,N] point positions.
    :param refpos: The reference position
    :param wind_dir: An (E,N) vector representing the wind direction.
    :return: a list of distances on the projected line
    """

    # Compute the wind unit vector
    wind_unit = np.array(wind_dir)/(np.linalg.norm(wind_dir))

    dists = []
    # Compute the distance along the line for each point
    for pt in pos:
        (dx, dy) = (pt[0] - refpos[0], pt[1] - refpos[1])
        dists.append(dx * wind_unit[0] + dy * wind_unit[1])  # Dot product

    return dists


if __name__ == "__main__":
    main()
