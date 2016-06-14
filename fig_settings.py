## Figure size:
# fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
# inches_per_pt = 1.0/72.27               # Convert pt to inches
# golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
# fig_width = fig_width_pt*inches_per_pt  # width in inches
# fig_height =fig_width*golden_mean       # height in inches
# fig_size = [fig_width,fig_height]
# rc("figure",figsize = fig_size)
text_color = "#000000"#"#424242"
ms_size = 50

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':18})
#rc('xtick', labelsize=16) 
#rc('ytick', labelsize=16) 
rc('text', usetex=True,color=text_color)
rc('axes',edgecolor=text_color)
#rc('text',color=text_color)
rc('xtick',color=text_color)
rc('ytick',color=text_color)
rc('lines',markeredgewidth=0.2)
rc('lines',markersize=ms_size*2)

myformat = "pdf"
## color code
color_MA_perturb = "#FF0040"
color_MM_marker = "#088A08"
color_MM_perturb = "#2E64FE"

myblue ="#5C85AD"
mygreen ="#01DF74"
myred ="#DB4D4D"
