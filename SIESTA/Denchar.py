"""
This file if for manupulation of denchar out files:
    We really dont need this since VESTA can visualize cube files which denchar can write
    Since finding that out I have to eddited this document and will load denchar files straight into VESTA
"""

def condition_data(*args, **kwargs):
    import glob,os
    filenames=[]
    # os.chdir("/mydir")
    for file in glob.glob("*MOD*"):
            filenames.append(file)
    print(filenames)

    data = []
    for i,filename in enumerate(filenames):
        print(f"Working on file: {filename}")
        with open(filename, "r") as file:
            for line in file:
                line = line.strip().split()
                line = [float(x) for x in line]
                if len(line) != 3:
                    continue
                else:
                    if i == 0:
                        data.append(line)
                    else:
                        for j in range(len(data)):
                            if line[0] == data[j][0] and line[1] == data[j][1]:
                                data[j][2] += line[2]
                                break

    with open (f"conditioned_data", "w+") as file_write:
        print(f"writing to data")
        for i in range(len(data)):
            file_write.write(f"{data[i][0]}  {data[i][1]}  {data[i][2]}\n")


class plot_WFSX():
    def __init__(self, filenames, *args, **kwargs):
        """
        filesnames: (list of strings) Will add all releveant Z values
        kwargs:
            levels: (list of ints) we can set levels to show numerical values in countour plots 
        """
        import pandas as pd
        import numpy as np
        # contour_data = pd.read_csv(filename, delim_whitespace=True, header=None)
        # contour_data.columns=["X", "Y", "Z"]
        for i,filename in enumerate(filenames):
            print(i)
            print(filename)
            contour_data = pd.read_csv(filename, delim_whitespace=True, header=None)
            contour_data.columns=["X", "Y", "Z"]

            Z = contour_data.pivot_table(index='X', columns='Y', values='Z').T.values
            X_unique = np.sort(contour_data.X.unique())
            Y_unique = np.sort(contour_data.Y.unique())
            X, Y = np.meshgrid(X_unique, Y_unique)

            import matplotlib.pyplot as plt
            from matplotlib import rcParams

            # Initialize plot objects
            rcParams['figure.figsize'] = 6, 6 # sets plot size
            fig = plt.figure()
            ax = fig.add_subplot(111)

            # Define levels in z-axis where we want lines to appear
            levels = kwargs.get("levels", [])
            # levels = np.array([-0.4,-0.2,0,0.2,0.4])

            # Generate a color mapping of the levels we've specified
            import matplotlib.cm as cm # matplotlib's color map library
            cpf = ax.contourf(X,Y,Z, len(levels), cmap=cm.RdBu)

            # Set all level lines to black
            line_colors = ['black' for l in cpf.levels]

            # Make plot and customize axes
            # cp = ax.contour(X, Y, Z, levels=levels, colors=line_colors)
            # cp = ax.contour(X, Y, Z, levels=levels)
            # ax.clabel(cp, fontsize=10, colors=line_colors)
            # plt.xticks([0,0.5,1])
            # plt.yticks([0,0.5,1])
            ax.set_xlabel('X-axis')
            plt.title("$\rho/cm^2$")
            _ = ax.set_ylabel('Y-axis')
            plt.savefig(f'{filename}.pdf') # uncomment to save vector/high-res version