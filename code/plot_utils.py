import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import subprocess
import PIL

def initialize_3d(fig_size=1080, ratio=1, subplots=(1, 1), facecolor='#F8F9E6',
                  elev=10, azim=30, l=4):

    fig_w, fig_h = fig_size * ratio, fig_size
    dpi = 300
    fig_width = fig_w / dpi
    fig_height = fig_h / dpi
    fig_size = fig_width * fig_height
    fs = np.sqrt(fig_size)
        
    fig = plt.figure(
        figsize=(fig_width, fig_height),
        dpi=dpi,  # Default dpi, will adjust later for saving
        layout='none',
    )
    fig.patch.set_facecolor(facecolor)

    gs = mpl.gridspec.GridSpec(subplots[0], subplots[1], figure=fig)
    axs = [[None] * subplots[1] for _ in range(subplots[0])]

    for i in range(subplots[0]):
        for j in range(subplots[1]):
            axs[i][j] = fig.add_subplot(gs[i, j], projection='3d')
            ax = axs[i][j]

            ax.set_facecolor(facecolor)


            ax.view_init(elev=elev, azim=azim)

            ax.set_xlim(-l[0], l[0])
            ax.set_ylim(-l[1], l[1])
            ax.set_zlim(-l[2], l[2])
            ax.grid(False)
            for k, axis in enumerate([ax.xaxis, ax.yaxis, ax.zaxis]):
                axis.line.set_linewidth(0 * fs)
                axis.set_visible(False)
                axis.pane.fill = False
                axis.pane.set_edgecolor('none')

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])

            dz = 2*l[2]
            dx = 2*l[0]
            dy = 2*l[1]

            ax.set_box_aspect([1,dy/dx,dz/dx])

    return fig, axs, fs, gs

def plot_vel_distr(v_phi, v_R, v_z, pos, v_lim=500):

    fig_size, rat = 540, 1.5
    subplots=(2,2)
    fig, axs, fs, gs = initialize_figure_old(fig_size=fig_size, ratio=rat, ts=2, sw=0.25, subplots=subplots, wmerge=[0], minor=True, wr=[1,1])

    bins = np.linspace(-50, v_lim)

    vels = [v_phi, v_R, v_z]
    v_names = ['$v_\phi$', '$v_R$', '$v_z$']
    for i, vel in enumerate(vels):
        counts_phi, bins = np.histogram(vel, bins=bins)
        max_count_phi = counts_phi.max()
        normalized_counts = counts_phi / max_count_phi
        axs[0][0].bar(bins[:-1], normalized_counts, width=np.diff(bins), align='edge', alpha=0.85, label=v_names[i])

    # plot legend
    axs[0][0].legend(loc='upper right', fontsize=fs*2)

    # contour plot
    xs = pos[0, :].value
    ys = pos[1, :].value
    zs = pos[2, :].value

    # 2d histogram from xs and ys
    H, xedges, yedges = np.histogram2d(xs, ys, bins=50)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H == 0, H)  # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    axs[0][1].pcolormesh(xedges, yedges, Hmasked)
    axs[0][1].set_aspect('equal')

    # 2d histogram from xs and zs
    H, xedges, yedges = np.histogram2d(xs, zs, bins=50)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H == 0, H)  # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    axs[1][1].pcolormesh(xedges, yedges, Hmasked)
    axs[1][1].set_aspect('equal')

    # set y axes to the right
    axs[0][1].yaxis.tick_right()
    axs[1][1].yaxis.tick_right()

    # set x y labels
    axs[0][1].set_xlabel('$x$', fontsize=2*fs)
    axs[0][1].set_ylabel('$y$', fontsize=2*fs, labelpad=0)
    axs[1][1].set_xlabel('$x$', fontsize=2*fs)
    axs[1][1].set_ylabel('$z$', fontsize=2*fs, labelpad=0)

    axs[0][1].xaxis.set_label_position("top")
    axs[0][1].yaxis.set_label_position("right")
    axs[1][1].yaxis.set_label_position("right")




    return 

def png_to_mp4(fold, title='video', fps=36, digit_format='04d', res=(1920, 1080)):

    # Get a list of all .mp4 files in the directory
    files = [f for f in os.listdir(fold) if f.endswith('.png')]
    files.sort()

    im = PIL.Image.open(fold+files[0])
    resx, resy = im.size

    resx = int(3*resx)//4
    resy = int(3*resy)//4

    if resx % 2 == 1:
        resx += 1
    if resy % 2 == 1:
        resy += 1

    name = os.path.splitext(files[0])[0]
    basename = name.split('_')[0]

    ffmpeg_path = 'ffmpeg'  
    framerate = fps
    abs_path = os.path.abspath(fold)
    parent_folder = os.path.dirname(abs_path)+'\\'
    output_file = parent_folder + "{}.mp4".format(title,framerate)
    crf = 18
    bitrate = "5000k"  # 5 Mbps, adjust as needed
    preset = "slow"  # Adjust as per your patience; slower will be better quality but take longer
    tune = "film"  # Especially if your content is more graphical/animated

    command = f'{ffmpeg_path} -y -r {framerate} -i {fold}{basename}_%{digit_format}.png -c:v libx264 -crf {crf} -preset {preset} -tune {tune} -pix_fmt yuv420p -vf scale={resx}:{resy} {output_file}'

    subprocess.run(command, shell=True)

def png_to_gif(fold, title='video', outfold=None, fps=36, 
               digit_format='04d', quality=500):

    # Get a list of all .mp4 files in the directory
    files = [f for f in os.listdir(fold) if f.endswith('.png')]
    files.sort()

    name = os.path.splitext(files[0])[0]
    basename = name.split('_')[0]

    ffmpeg_path = 'ffmpeg'  
    framerate = fps

    if outfold==None:
        abs_path = os.path.abspath(fold)
        parent_folder = os.path.dirname(abs_path)+'\\'
    else:
        parent_folder = outfold
        if not os.path.exists(parent_folder):
            os.makedirs(parent_folder)
            
    output_file = parent_folder + "{}.gif".format(title,framerate)

    # Create a palette for better GIF quality
    palette_file = parent_folder + "palette.png"
    palette_command = f'{ffmpeg_path} -i {fold}{basename}_%{digit_format}.png -vf "fps={framerate},scale={quality}:-1:flags=lanczos,palettegen" -y {palette_file}'
    subprocess.run(palette_command, shell=True)
    print(palette_file)

    # Use the palette to create the GIF
    gif_command = f'{ffmpeg_path} -r {framerate} -i {fold}{basename}_%04d.png -i {palette_file} -lavfi "fps={framerate},scale={quality}:-1:flags=lanczos [x]; [x][1:v] paletteuse" -y {output_file}'
    subprocess.run(gif_command, shell=True)

def initialize_figure_old(
    fig_size=540,
    ratio=1.5,
    fig_w=None, fig_h=None,
    subplots=(1, 1), grid=True, 
    lw=0.015, ts=2, theme=None,
    pad=0.5,
    color='#222222',
    dpi=300,
    sw=0.15,
    wr=None, hr=None, hmerge=None, wmerge=None,
    ylabel='bottom',
    layout='constrained',
    hspace=None, wspace=None,
    tick_direction='inout',
    minor=True,
    top_bool=True,
    projection=None
):

    """
    This is a generic function to initialize a figure with a given size and
    number of subplots. It's just for the mass density and velocity distribution plot. It's a bit too complex for that purpose, but I already had it written.   
    """

    if fig_w is None:
        fig_w = fig_size * ratio
        fig_h = fig_size

    fig_width = fig_w / dpi
    fig_height = fig_h / dpi
    fig_size = fig_width * fig_height
    fs = np.sqrt(fig_size)
    fig = plt.figure(
        figsize=(fig_width, fig_height),
        dpi=dpi,  # Default dpi, will adjust later for saving
        layout=layout,
    )

    if wr is None:
        wr_ = [1] * subplots[1]
    else:
        wr_ = wr
    if hr is None:
        hr_ = [1] * subplots[0]
    else:
        hr_ = hr
    

    gs = mpl.gridspec.GridSpec(subplots[0], subplots[1], figure=fig, width_ratios=wr_, height_ratios=hr_, hspace=hspace, wspace=wspace)


    ax = [[None] * subplots[1] for _ in range(subplots[0])]

    if theme == "dark":
        fig.patch.set_facecolor(color)
        plt.rcParams.update({"text.color": "white"})

    for i in range(subplots[0]):
        for j in range(subplots[1]):
            
            if hmerge is not None:
                if i in hmerge:
                    ax[i][j] = fig.add_subplot(gs[i, :])
                else:
                    ax[i][j] = fig.add_subplot(gs[i, j])
            elif wmerge is not None:
                if j in wmerge:
                    ax[i][j] = fig.add_subplot(gs[:, j])
                else:
                    ax[i][j] = fig.add_subplot(gs[i, j])
            else:
                ax[i][j] = fig.add_subplot(gs[i, j], projection=projection)

            if theme == "dark":
                ax[i][j].set_facecolor(color)
                ax[i][j].tick_params(colors="white")
                ax[i][j].spines["bottom"].set_color("white")
                ax[i][j].spines["top"].set_color("white")
                ax[i][j].spines["left"].set_color("white")
                ax[i][j].spines["right"].set_color("white")
                ax[i][j].xaxis.label.set_color("white")
                ax[i][j].yaxis.label.set_color("white")

            #ax[i][j].xaxis.set_tick_params(which="minor", bottom=False)

            if grid:
                ax[i][j].grid(
                    which="major",
                    linewidth=fs * lw,
                    color="white" if theme == "dark" else "black",
                )
            for spine in ax[i][j].spines.values():
                spine.set_linewidth(fs * sw)

            if ylabel == 'bottom':
                labeltop_bool = False
                labelbottom_bool = True
            elif ylabel == 'top':
                labeltop_bool = True
                labelbottom_bool = False
                ax[i][j].xaxis.set_label_position('top')

            else:
                labeltop_bool = True
                labelbottom_bool = True
                ax[i][j].xaxis.set_label_position('both')

            
            ax[i][j].tick_params(
                axis="both",
                which="major",
                labelsize=ts * fs,
                size=fs * sw*2,
                width=fs * sw,
                pad= pad * fs,
                top=top_bool,
                labelbottom=labelbottom_bool,
                labeltop=labeltop_bool,
                right=top_bool,
                direction=tick_direction
            )

            if minor:
                ax[i][j].minorticks_on()
                ax[i][j].tick_params(
                axis="both",
                which="major",
                labelsize=ts * fs,
                size=fs * sw*4,
                width=fs * sw,
                pad= pad * fs,
                top=top_bool,
                labelbottom=labelbottom_bool,
                labeltop=labeltop_bool,
                right=top_bool,
                direction=tick_direction
                )
                ax[i][j].tick_params(axis='both', which="minor", 
                direction=tick_direction,
                top=top_bool,
                right=top_bool,
                size=fs * sw*2.5, width=fs * sw,)

    if hmerge is not None:
        for k in hmerge:
            for l in range(1, subplots[1]):
                fig.delaxes(ax[k][l])

    if wmerge is not None:
        for k in wmerge:
            for l in range(1, subplots[0]):
                fig.delaxes(ax[l][k])
            
    
    return fig, ax, fs, gs