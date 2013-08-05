#!/usr/bin/env python
'''
Tool used for visualization and analysis of the COS FUV detectors.

Auther: Justin Ely
Version: 1.0 ish
'''

import pyfits
import pylab
import sys
import numpy
import glob
import os

from ttag_funcs import ttag_image
from astroraf.math.utils import gauss_kern, blur_image
from astroraf.plotting import init_plots

from matplotlib.patches import Rectangle
from Tkinter import Tk, IntVar, StringVar, Frame, Button, Menu, Label, Radiobutton, OptionMenu, W, N, S, E, Checkbutton, Entry, END
import tkFileDialog, Tkconstants
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
 
data_dir = '/user/ely/COS/COS_layout/'

def find_defects(data, thresh=.5, which='light'):
        seg=self.segment.get()
        xlim_a = (1176, 15219)
        ylim_a = (414, 558)
        xlim_b = (870, 14984)
        ylim_b = (470, 608)
        ylen = 10
        xlen = 6

        if seg == 'A':
            infile = os.path.join( data_dir, '12676-all-A-2.fits' )
            lx = xlim_a[0]
            ux = xlim_a[1]
            ly = ylim_a[0]
            uy = ylim_a[1]
        elif seg == 'B':
            infile = os.path.join( data_dir, '12676-LV-B-2.fits' )
            lx = xlim_b[0]
            ux = xlim_b[1]
            ly = ylim_b[0]
            uy = ylim_b[1]

        if which=='dark':
            ylen=10
            xlen=6
            if seg=='A':
                dark=pyfits.getdata('all_pha_0.fits', 1) 
                lx=1172
                ly=340
                ux=15150
                uy=720
            elif seg=='B':
                dark=pyfits.getdata('all_pha_0.fits', 2)
                lx=1040
                ly=402
                ux=14945
                uy=754
            print 'Blurring array'
            array=blur_image(dark, 10, 6)
        else:
            array=pyfits.getdata(infile, 1)

        coords=[]
        for ystep in range(ly, uy-ylen):
            print ystep
            for xstep in range(lx, ux-10):
                sub_array=array[ystep:ystep+ylen, xstep:xstep+xlen]
                med=numpy.median(sub_array)
                std=numpy.std(sub_array)
                if which=='dark':
                    index=numpy.where(sub_array>((med)*thresh))
                else:
                    index=numpy.where(sub_array<((med)*thresh))
                if len(index[0])>0:
                    for i in range(len(index[0])):
                        coords.append((index[0][i]+ystep, index[1][i]+xstep))

        coord_set=set(coords)
        final_coords=[]
        length=len(coord_set)
        for i, item in enumerate(coord_set):
            print i, '/', length
            if (coords.count(item)>=30 and self.has_neighbor(item, coord_set)):
                print 'True'
                final_coords.append(item)
        if len(final_coords)>0:
            xs=[line[1] for line in final_coords]
            ys=[line[0] for line in final_coords]

def find_contiguous(input_array):
    """DocString
    """
    try: np
    except : import numpy as np
    regions = []
    padding = 1
 
    shape = input_array.shape
    valid_index = np.where( input_array > 0 )
    found_array = np.zeros( (shape) )
    for i in xrange( shape[0] ):
        #print i
        for j in xrange( shape[1] ):
            if input_array[i, j] and not found_array[i, j]:
                points = find_extent(input_array, j, i)
		print i, j, len(points)
                if len(points) <= 4: continue
                xs = numpy.array( [item[0] for item in points] )
                ys = numpy.array( [item[1] for item in points] )
                lx = xs.min() - padding
                ly = ys.min() - padding
                dx = ( xs.max() + padding - lx ) +1
                dy = ( ys.max() + padding - ly ) +1
		
                regions.append( (lx, dx, ly, dy) )
		found_array[ly:ly+dy, lx:lx+dx] = 1

    return list( set( regions) )
                    

def get_neighbors(input_array, x, y):
    neighbors = []
    for i in [-2, -1, 0, 1, 2]:
        for j in [-2, -1, 0, 1, 2]:
            if i != j:
                if input_array[y+i, x+j]:
                    neighbors.append( (x+j,  y+i) )

    return neighbors


def find_extent(input_array, initial_x, initial_y):
    MORE_TO_CHECK = True
    all_points = [ (initial_x, initial_y) ]

    while MORE_TO_CHECK:
        xs = [ item[0] for item in all_points ]
        ys = [ item[1] for item in all_points ]

        starting_size = len(all_points)
        for x, y in zip(xs, ys):
            new_points = get_neighbors(input_array, x, y)
            for point in new_points:
                if not point in all_points:
                    all_points.append( point )

        ending_size = len( all_points )
        if starting_size == ending_size: MORE_TO_CHECK = False

    return all_points

class App:
    def __init__(self,  parent):
        #----Variables----#
        self.show_gain=IntVar()
        self.show_illumination=IntVar()
        self.show_spectrum=IntVar()
        self.show_dark=IntVar()
        self.show_dq=IntVar()
        self.show_mydq=IntVar()

        self.xmin=StringVar()
        self.xmax=StringVar()
        self.ymin=StringVar()
        self.ymax=StringVar()
        self.vmin=StringVar()
        self.vmax=StringVar()
        self.vmin.set(0)
        self.vmax.set(0)

        self.segment=StringVar()
        self.segment.set('A')
        self.extract=StringVar()
        self.extract.set('None')
        self.extract_offset=IntVar()
        self.extract_offset.set(0)
        self.cmap=StringVar()
        self.cmap.set('gist_yarg')
        self.dq_show=StringVar()
        self.dq_show.set('184')
        self.infile=StringVar()
        self.draw=StringVar()
        self.draw.set('Modal Gain')
        self.degrade_loc=StringVar()
        self.degrade_years=StringVar()
        self.degrade_width=StringVar()
        self.degrade_xshift=StringVar()
        self.degrade_array=IntVar()
	self.N_degraded=IntVar()
	self.N_degraded.set(0)
        self.grid_limits=IntVar()
        self.grid_limits.set(1)

	self.data_file = None
        self.Againmap=None
        self.Bgainmap=None
        self.Aexp=None
        self.Bexp=None
        self.Adark=None
        self.Bdark=None

        self.last=None

        #-----Open files to speed things up----#
        self.myParent = parent
        self.myContainer = Frame(parent)
        self.myParent.bind('<Return>', self.update_plot)
        self.myParent.title('COS Layout Planner')



        #------------Master draw button-----------#
        self.redraw = Button(self.myContainer, text="Draw Map", command=self.draw_all)
        self.redraw.grid(row=1, column=0, sticky=N)
        
        #-----------Define Clear Plot Button--------#
        self.clear_button = Button(self.myContainer, text="Clear Plot", command=self.clear_axis)
        self.clear_button.grid(row = 2, column = 0)

        #-----Save button---------#
        #self.savefig=Button(self.myContainer, text='Save Figure', command=self.save_fig)
        #self.savefig.grid(row=2, column=0)
        #----Now in menu----------#

        #-----------Define Quit Button--------#
        #self.quitbutton = Button(self.myContainer, text="QUIT", fg="red", command=self.myContainer.quit)
        #self.quitbutton.grid(row = 1, column = 1)
        #---------Now in Menu-----------------#

        #----------Menu-------------#
        menubar = Menu( parent )

        # create a pulldown menu, and add it to the menu bar
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open", command=self.open_file)
        filemenu.add_command(label="Save", command=self.save_current)
        filemenu.add_separator()

        filemenu.add_command(label="Exit", command=self.exit)
        menubar.add_cascade(label="File", menu=filemenu)
        
        # create more pulldown menus
        toolmenu = Menu(menubar, tearoff=0)
        toolmenu.add_command(label="Open Toolbox", command=self.dq_tools)
        toolmenu.add_command(label="Find Low", command=self.find_low)
        #toolmenu.add_command(label="Paste", command=self.hello)
        menubar.add_cascade(label="Tools", menu=toolmenu)
        
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Show Help", command=self.help_file)
        menubar.add_cascade(label="Help", menu=helpmenu)
        #submenu
        cmapmenu=Menu(toolmenu, tearoff=0)
        cmapmenu.add_radiobutton(label='Grey', variable=self.cmap, value='gist_yarg')
        cmapmenu.add_radiobutton(label='Prism', variable=self.cmap, value='prism')
        cmapmenu.add_radiobutton(label='Jet', variable=self.cmap, value='jet')
        cmapmenu.add_radiobutton(label='Ncar', variable=self.cmap, value='gist_ncar')
        #cmapmenu.add_radiobutton(label='autumn', variable=self.cmap, value='autumn', command=self.set_cmap(self.cmap))
        #cmapmenu.add_radiobutton(label='spring', variable=self.cmap, value='spring', command=self.set_cmap(self.cmap))
        toolmenu.insert_cascade(index=2, label='Select Cmap', menu=cmapmenu)

        self.menubar=menubar
        # display the menu
        parent.config(menu=self.menubar)

        
        #----Segment Button---#
        self.segment_selector_label= Label(self.myContainer, text = 'Segment Selector')
        self.segment_selector_label.config(fg = 'blue')
        self.segment_selector_label.grid(row=0, column=1, padx=3, pady=3)
        self.segment_selector=Radiobutton(self.myContainer, text="A", variable=self.segment, value='A', indicatoron=1).grid(row=1, column=1, padx=3, pady=3)
        self.segment_selector=Radiobutton(self.myContainer, text="B", variable=self.segment, value='B', indicatoron=1).grid(row=2, column=1, padx=3, pady=3)
        
        #-----------Define drawing boxes------------#
        self.drawing_label= Label(self.myContainer, text = 'Map to draw')
        self.drawing_label.config(fg = 'blue')
        self.drawing_label.grid(row=0, column=2, padx=3, pady=3)

        self.draw_what=OptionMenu(self.myContainer, self.draw, 'Data File', 'Modal Gain', 'Summed Exposure', 'Summed Dark')
        self.draw_what.grid(row =1, column = 2, sticky = W+E, padx = 3, pady = 3)

        #self.toggle_gain = Checkbutton(self.myContainer, text="Modal Gain Map", variable=self.show_gain)
        #self.toggle_gain.grid(row =1, column = 3, sticky = W, padx = 3, pady = 3)

        #self.degrade_label= Label(self.myContainer, text = 'Degrade Gain')
        #self.degrade_label.config(fg = 'blue')
        #self.degrade_label.grid(row=2, column=2, padx=3, pady=3)
        self.toggle_degrade=Checkbutton(self.myContainer, text="Degrade Gain Map", variable=self.degrade_array)
        self.toggle_degrade.grid(row=2, column=2, sticky=W+E, padx=3, pady=3)

        self.degrade_loc_entry=Entry(self.myContainer, textvariable=self.degrade_loc)
        self.degrade_loc_entry.grid(row=3, column=2, sticky=W+E, padx=3, pady=3)
        self.degrade_loc_label=Label(self.myContainer, text='Y Loc')
        self.degrade_loc_label.config(fg = 'blue')
        self.degrade_loc_label.grid(row=3, column=1, sticky=E, padx=3, pady=3)

        self.degrade_years_entry=Entry(self.myContainer, textvariable=self.degrade_years)
        self.degrade_years_entry.grid(row=4, column=2, sticky=W+E, padx=3, pady=3)
        self.degrade_years_label=Label(self.myContainer, text='Years')
        self.degrade_years_label.config(fg = 'blue')
        self.degrade_years_label.grid(row=4, column=1, sticky=E, padx=3, pady=3)

        self.degrade_width_entry=Entry(self.myContainer, textvariable=self.degrade_width)
        self.degrade_width_entry.grid(row=5, column=2, sticky=W+E, padx=3, pady=3)
        self.degrade_width_label=Label(self.myContainer, text='Width')
        self.degrade_width_label.config(fg = 'blue')
        self.degrade_width_label.grid(row=5, column=1,sticky=E, padx=3, pady=3)

        self.degrade_xshift_entry=Entry(self.myContainer, textvariable=self.degrade_xshift)
        self.degrade_xshift_entry.grid(row=6, column=2, sticky=W+E, padx=3, pady=3)
        self.degrade_xshift_label=Label(self.myContainer, text='X shift')
        self.degrade_xshift_label.config(fg = 'blue')
        self.degrade_xshift_label.grid(row=6, column=1, sticky=E, padx=3, pady=3)



        #self.degrade_years_entry.grid_remove()
        #self.degrade_loc_entry.grid_remove()
        #self.toggle_degrade.grid_remove()

        #self.toggle_degrade.grid_remove()

        #self.toggle_limits=Radiobutton(self.myContainer, text="Set Plot Limits", indicatoron=0, variable=self.grid_limits, value=1, command=self.show_limit_boxes).grid(row=5, column=0)
        #self.toggle_limits=Radiobutton(self.myContainer, text="Set Degrade Params", indicatoron=0, variable=self.grid_limits, value=2, command=self.show_limit_boxes).grid(row=5, column=1)
        #self.toggle_limits=Radiobutton(self.myContainer, text="Add DQ Flags", indicatoron=0, variable=self.grid_limits, value=3, command=self.show_limit_boxes).grid(row=5, column=2)

        '''
        self.toggle_illumination = Checkbutton(self.myContainer, text="Summed Exposure", variable=self.show_illumination)
        self.toggle_illumination.grid(row =2, column = 3, sticky = W, padx = 3, pady = 3)

        self.toggle_dark = Checkbutton(self.myContainer, text="Summed Dark", variable=self.show_dark)
        self.toggle_dark.grid(row =1, column = 4, sticky = W, padx = 3, pady = 3)

        #self.toggle_spec = Checkbutton(self.myContainer, text="Sample spectrum", variable=self.show_spectrum)
        #self.toggle_spec.grid(row =3, column = 3, sticky = W, padx = 3, pady = 3)
        '''

        self.toggle_dq = Checkbutton(self.myContainer, text="DQ Flags", variable=self.show_dq, indicatoron=1)
        self.toggle_dq.grid(row =2, column = 3, sticky = W, padx = 3, pady = 3)
        
        self.toggle_mydq = Checkbutton(self.myContainer, text="New DQ Flags", variable=self.show_mydq, indicatoron=1)
        self.toggle_mydq.grid(row =3, column = 3, sticky = W, padx = 3, pady = 3)

        #--------DQ flags to show-------#
        self.dq_show_label=Label(self.myContainer, text = 'DQ to plot (bitwise)')
        self.dq_show_label.config(fg = 'blue')
        self.dq_show_box=Entry(self.myContainer, textvariable=self.dq_show)

        self.dq_show_box.grid(row=4, column=3, sticky=W+E, padx = 3, pady = 3)
        self.dq_show_label.grid(row=0, column=3, sticky=W+E, padx = 3, pady = 3)
        self.dq_show_box.delete(0, END)
        self.dq_show_box.insert(0, 184)

        #self.toggle_dq.grid_remove()
        #self.toggle_mydq.grid_remove()
        #self.dq_show_label.grid_remove()
        #self.dq_show_box.grid_remove()


        #------------Master draw button-----------#
        self.draw_dq = Button(self.myContainer, text="Add DQ", command=self.add_dq)
        self.draw_dq.grid(row=1, column=3, sticky=N)

        #---------------------------------------------------------------------#
        # assign two enter boxes for the low color bin and the high color bin #

        self.xmin_box_label= Label(self.myContainer, text = 'Xmin')
        self.xmin_box_label.config(fg = 'blue')
        self.xmin_box = Entry(self.myContainer, textvariable=self.xmin)

        self.xmax_box_label=Label(self.myContainer, text = 'Xmax')
        self.xmax_box_label.config(fg = 'blue')
        self.xmax_box=Entry(self.myContainer, textvariable=self.xmax)

        self.xmin_box.grid(row =1, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.xmax_box.grid(row =1, column = 7, sticky = W+E, padx = 3, pady = 3)
        self.xmin_box_label.grid(row =0, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.xmax_box_label.grid(row =0, column = 7, sticky = W+E, padx = 3, pady = 3)
        self.xmin_box.delete(0, END)
        self.xmax_box.delete(0, END)
        self.xmin_box.insert(0, 'ALL')
        self.xmax_box.insert(0, 'ALL')

        self.ymin_box_label= Label(self.myContainer, text = 'Ymin')
        self.ymin_box_label.config(fg = 'blue')
        self.ymin_box = Entry(self.myContainer, textvariable=self.ymin)

        self.ymax_box_label=Label(self.myContainer, text = 'Ymax')
        self.ymax_box_label.config(fg = 'blue')
        self.ymax_box=Entry(self.myContainer, textvariable=self.ymax)

        self.ymin_box.grid(row =3, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.ymax_box.grid(row =3, column = 7, sticky = W+E, padx = 3, pady = 3)
        self.ymin_box_label.grid(row =2, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.ymax_box_label.grid(row =2, column = 7, sticky = W+E, padx = 3, pady = 3)
        self.ymin_box.delete(0, END)
        self.ymax_box.delete(0, END)
        self.ymin_box.insert(0, 'ALL')
        self.ymax_box.insert(0, 'ALL')


        #self.ymin_box.grid_remove()
        #self.ymax_box.grid_remove()
        #self.ymin_box_label.grid_remove()
        #self.ymax_box_label.grid_remove()

        #self.xmin_box.grid_remove()
        #self.xmax_box.grid_remove()
        #self.xmin_box_label.grid_remove()
        #self.xmax_box_label.grid_remove()


        #---show extraction region box-------#
        '''
        self.extract_label=Label(self.myContainer, text = 'SEG/OPT_ELEM/CENWAVE/APER')
        self.extract_label.config(fg = 'blue')
        self.extract_box=Entry(self.myContainer, textvariable=self.extract)

        self.extract_box.grid(row=5, column=5, sticky=W+E, padx = 3, pady = 3)
        self.extract_label.grid(row=4, column=5, sticky=W+E, padx = 3, pady = 3)
        self.extract_box.delete(0, END)
        self.extract_box.insert(0, 'None')
        

        #--------offset box-------#
        self.extract_offset_label=Label(self.myContainer, text = 'Offset extraction box (pixels)')
        self.extract_offset_label.config(fg = 'blue')
        self.extract_offset_box=Entry(self.myContainer, textvariable=self.extract_offset)

        self.extract_offset_box.grid(row=5, column=6, sticky=W+E, padx = 3, pady = 3)
        self.extract_offset_label.grid(row=4, column=6, sticky=W+E, padx = 3, pady = 3)
        self.extract_offset_box.delete(0, END)
        self.extract_offset_box.insert(0, 0)
        '''
        #--------vmin,vmax boxes---------#

        self.vmin_box_label= Label(self.myContainer, text = 'Contrast min')
        self.vmin_box_label.config(fg = 'blue')
        self.vmin_box = Entry(self.myContainer, textvariable=self.vmin)

        self.vmax_box_label=Label(self.myContainer, text = 'Contrast max')
        self.vmax_box_label.config(fg = 'blue')
        self.vmax_box=Entry(self.myContainer, textvariable=self.vmax)

        self.vmin_box.grid(row =5, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.vmax_box.grid(row =5, column = 7, sticky = W+E, padx = 3, pady = 3)
        self.vmin_box_label.grid(row =4, column = 6, sticky = W+E, padx = 3, pady = 3)
        self.vmax_box_label.grid(row =4, column = 7, sticky = W+E, padx = 3, pady = 3)


        #self.vmax_box_label.grid_remove()
        #self.vmax_box.grid_remove()
        #self.vmin_box_label.grid_remove()
        #self.vmin_box.grid_remove()
        
        #-------Plots-------#
        init_plots()
        pylab.ioff()
        self.fig=pylab.figure(1, figsize=(20, 8.5))
        self.cax=self.fig.add_axes([.92, .1, .01, .82])
        self.ax=self.fig.add_axes()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.myContainer)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row = 7,  column = 0, rowspan=10, columnspan=15, sticky=W+E+S) 
        #self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, parent)
        self.toolbar.update()
        #self.canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=1)
        self.canvas._tkcanvas.grid(row = 7, column = 0,rowspan=10,columnspan=15,sticky=W+E+S)
        self.canvas.show()

        self.myContainer.pack()

    def add_dq(self):
        seg=self.segment.get()
        colors={2:'b',4:'g',8:'y',16:'r',32:'c',1024:'gold',4096:'m',8192:'orange'}
        if self.show_dq.get():
            dq_to_show=int(self.dq_show.get())
            bpix=pyfits.open('bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg and (line[5]&dq_to_show):
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']
                    dy=line['DY']
                    dq=line['DQ']
                    x=[lx,lx+dx,lx+dx,lx,lx]
                    y=[ly,ly,ly+dy,ly+dy,ly]
                    pylab.plot(x,y,'-',lw=2,color=colors[int(dq)])
                    pylab.annotate(str(dq),(lx,ly),color=colors[int(dq)])
        if self.show_mydq.get():
            dq_to_show=int(self.dq_show.get())
            bpix=pyfits.open('new_bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg and (line[5]&dq_to_show):
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']
                    dy=line['DY']
                    dq=line['DQ']
                    x=[lx,lx+dx,lx+dx,lx,lx]
                    y=[ly,ly,ly+dy,ly+dy,ly]
                    pylab.plot(x,y,'--',lw=2,color=colors[int(dq)])
                    pylab.annotate(str(dq),(lx,ly+dy),color=colors[int(dq)])
        self.canvas.show()

    def clear_axis(self):
        pylab.figure(1)
        pylab.subplot(1,1,1)
        pylab.cla()
        #self.toggle_dq.deselect()
        #self.toggle_mydq.deselect()
        #self.toggle_spec.deselect()
        #self.canvas.delete(all)
        self.canvas.show()
        self.extract.set('None')
        #self.Againmap.close()
        #self.Bgainmap.close()
        self.Againmap=None
        self.Bgainmap=None
	self.N_degraded.set(0)

    def degrade(self,array):
        '''
        Adapted from code from D. Massa
        '''
        date0=(55772.0-54983.0)/365.25
        date1=(56108.0-54983.0)/365.25
        xlim_a=(1200,15099)
        ylim_a=(335,699)
        xlim_b=(950,15049)
        ylim_b=(400,749)
        limits=(300,16300)
        seg=self.segment.get()
        #years=(56108-54983.0)/365.25
        years=float(self.degrade_years.get())
        width=int(self.degrade_width.get())
        y=int(self.degrade_loc.get())
        xshift=int(self.degrade_xshift.get())
        if seg=='A':
            infile = os.path.join( data_dir, '12676-all-A-2.fits' )
            y0=487
        elif seg=='B':
            infile = os.path.join( data_dir, '12676-LV-B-2.fits' )
            y0=546
        gain=pyfits.getdata(infile,2)
        gain0=pyfits.getdata(infile,3)
        delta=(gain0-gain)/date0
	if not self.N_degraded.get():
            degraded=gain0-delta*date1
	    self.N_degraded.set(1)
	else:
	    degraded=array
        delta=numpy.roll(delta,shift=xshift)
        delta=numpy.roll(delta,shift=y-y0,axis=0)
        #degraded[y-(width-1)/2:y+(width-1)/2,:]-=delta[y0-(width-1)/2:y0+(width-1)/2,:]*years
        degraded-=delta*years


	hdu=pyfits.HDUList(pyfits.PrimaryHDU())
        hdu[0].header.update('TELESCOP','HST')
        hdu[0].header.update('INSTRUME','COS')
        hdu[0].header.update('DETECTOR','FUV')
        #hdu[0].header.update('OPT_ELEM','ANY')
        #hdu[0].header.update('FILETYPE','')
        #hdu[0].header.update('EXPTIME',exptime)
        #hdu[0].header.update('PHA_MIN',2)
        #hdu[0].header.update('PHA_MAX',30)
        hdu.append(pyfits.core.ImageHDU(data=degraded))
        hdu.writeto('gainmap.fits')


        return degraded

    def dq_tools(self):
        self.toolswindow=Toplevel()
        self.frame2=Frame(master=self.toolswindow)
        self.toolswindow.title('DQ tools')
        #-----------variables------------------#
        self.height=StringVar()
        self.height.set('500')
        self.width=StringVar()
        self.width.set('20')
        self.draw_DQ=IntVar()
        self.draw_DQ.set(0)
        self.draw_newDQ=IntVar()
        self.draw_newDQ.set(1)
        self.SDQFLAG=IntVar()
        self.SDQFLAG.set(184)
        self.segment2=StringVar()
        self.segment2.set('A')
        self.affected_pixels=StringVar()
        self.affected_pixels.set('0')
        self.affected_columns=StringVar()
        self.affected_columns.set('0')
        #-----------Define Close Button--------#
        self.closebutton2 = Button(self.toolswindow, text="CLOSE", fg="red", command=self.toolswindow.destroy)
        self.closebutton2.grid(row = 35, column = 6)
        self.extractbutton = Button(self.toolswindow, text="EXTRACT", fg="blue", command=self.extract_dqs)
        self.extractbutton.grid(row = 0, column = 0)
        self.savebutton = Button(self.toolswindow, text="Save Fig", fg="blue", command=self.save_fig)
        self.savebutton.grid(row = 35, column = 0)

        #-----------Find Best Postition--------#
        self.findbest = Button(self.toolswindow, text="Find best", fg="blue", command=self.find_best)
        self.findbest.grid(row = 0, column = 1)


        #-----------Entry Box-------------------#
        self.SDQFLAG_label=Label(self.toolswindow, text = 'SDQFLAG')
        self.SDQFLAG_label.config(fg = 'blue')
        self.SDQFLAG_box=Entry(self.toolswindow,textvariable=self.SDQFLAG)

        self.SDQFLAG_box.grid(row=6,column=0,sticky=W+E, padx = 3, pady = 3)
        self.SDQFLAG_label.grid(row=5,column=0,sticky=W+E, padx = 3, pady = 3)
        self.SDQFLAG_box.delete(0,END)
        self.SDQFLAG_box.insert(0,184)

        #-----------Extraction results---------------#
        self.results_label=Label(self.toolswindow, text = '% affected Pixels')
        self.results_label.config(fg = 'blue')
        self.results_box=Entry(self.toolswindow,textvariable=self.affected_pixels)

        self.results_box.grid(row=6,column=1,sticky=W+E, padx = 3, pady = 3)
        self.results_label.grid(row=5,column=1,sticky=W+E, padx = 3, pady = 3)
        self.results_box.delete(0,END)
        self.results_box.insert(0,'0')

        self.Cresults_label=Label(self.toolswindow, text = '% affected Columns')
        self.Cresults_label.config(fg = 'blue')
        self.Cresults_box=Entry(self.toolswindow,textvariable=self.affected_columns)

        self.Cresults_box.grid(row=6,column=2,sticky=W+E, padx = 3, pady = 3)
        self.Cresults_label.grid(row=5,column=2,sticky=W+E, padx = 3, pady = 3)
        self.Cresults_box.delete(0,END)
        self.Cresults_box.insert(0,'0')

        self.height_box_label= Label(self.toolswindow, text = 'Y Loc')
        self.height_box_label.config(fg = 'blue')
        self.height_box = Entry(self.toolswindow,textvariable=self.height)

        self.width_box_label=Label(self.toolswindow, text = 'Height')
        self.width_box_label.config(fg = 'blue')
        self.width_box=Entry(self.toolswindow,textvariable=self.width)

        self.height_box.grid(row =2, column = 0, sticky = W+E, padx = 3, pady = 3)
        self.width_box.grid(row =2, column = 1, sticky = W+E, padx = 3, pady = 3)
        self.height_box_label.grid(row =1, column = 0, sticky = W+E, padx = 3, pady = 3)
        self.width_box_label.grid(row =1, column = 1, sticky = W+E, padx = 3, pady = 3)
        self.height_box.delete(0,END)
        self.width_box.delete(0,END)
        self.height_box.insert(0,500)
        self.width_box.insert(0,'20')


        #----Segment Button---#
        self.segment_selector_label2= Label(self.toolswindow, text = 'Segment Selector')
        self.segment_selector_label2.config(fg = 'blue')
        self.segment_selector_label2.grid(row=7,column=0,padx=3,pady=3)
        self.segment_selector2=Radiobutton(self.toolswindow, text="A", variable=self.segment2, value='A').grid(row=8,column=0,padx=3,pady=3)
        self.segment_selector2=Radiobutton(self.toolswindow, text="B", variable=self.segment2, value='B').grid(row=9,column=0,padx=3,pady=3)

        #-----------DQ selectors---------------#
        self.draw_dq = Checkbutton(self.toolswindow, text="DQ Flags",variable=self.draw_DQ)
        self.draw_dq.grid(row =8, column = 1, sticky = W, padx = 3, pady = 3)
        
        self.draw_newdq = Checkbutton(self.toolswindow, text="New DQ Flags",variable=self.draw_newDQ)
        self.draw_newdq.grid(row =9, column = 1, sticky = W, padx = 3, pady = 3)

        matplotlib.rcParams['figure.subplot.left']=.15
        matplotlib.rcParams['figure.subplot.right']=.85
        matplotlib.rcParams['figure.subplot.top']=.85
        matplotlib.rcParams['figure.subplot.bottom']=.15
        self.fig2=pylab.figure(2,figsize=(4,6))
        self.ax2=self.fig2.add_axes()

        pylab.figure(2)
        self.canvas2 = FigureCanvasTkAgg(self.fig2,master=self.toolswindow)
        self.canvas2.show()
        self.canvas2.get_tk_widget().grid(row = 13, column = 0,rowspan=17,columnspan=17,sticky=N+W+E+S) 

	self.toolswindow.mainloop()

    def draw_all(self):
        #self.clear_axis()
        pylab.figure(1)
        pylab.cla()
        pylab.subplot(1,1,1)
	pylab.xlabel('X (Dispersion)')
        pylab.ylabel('Y (Cross Dispersion)')
        show_which = self.draw.get()
        seg=self.segment.get()        
	if ( show_which == 'Data File' ):
            vmin=int(self.vmin.get())
            vmax=int(self.vmax.get())
            if (vmin==0 and vmax==0):
		    vmin=1
	    vmax=15
            if vmin<0: vmin=0
            if vmax>20: vmax=20
	    self.vmin.set(vmin)
            self.vmax.set(vmax)
            levels=range(vmin,vmax+1)
	    if seg=='A':
                extension = 1
            elif seg=='B':
                extension = 2

            #data = pyfits.getdata(self.data_file,extension)
            #C1=pylab.contourf(gainmap.clip(max=vmax),levels,cmap=pylab.get_cmap(self.cmap.get()))
            C1 = pylab.imshow(ttag_image(data_file),interpolation='nearest',aspect='auto',cmap=pylab.get_cmap(self.cmap.get()))#,vmin=,vmax=vmax)
            pylab.colorbar(C1,cax=self.cax)

        if show_which=='Modal Gain':
            if self.Againmap==None: self.Againmap=pyfits.getdata(os.path.join( data_dir, '12676-all-A-2.fits' ),2)
            if self.Bgainmap==None: self.Bgainmap=pyfits.getdata(os.path.join( data_dir, '12676-LV-B-2.fits' ),2)
            vmin=int(self.vmin.get())
            vmax=int(self.vmax.get())
            if (vmin==0 and vmax==0):
                vmin=1
                vmax=15
            if vmin<0: vmin=0
            if vmax>20: vmax=20
            self.vmin.set(vmin)
            self.vmax.set(vmax)
            levels=range(vmin,vmax+1)
            if seg=='A':
                if self.degrade_array.get():self.Againmap=self.degrade(self.Againmap)
                gainmap=self.Againmap
            if seg=='B':
                if self.degrade_array.get():self.Bgainmap=self.degrade(self.Bgainmap)
                gainmap=self.Bgainmap
            #if self.degrade_array.get():
            #    gainmap=self.degrade()
            levels=range(1,vmax+1)
            C1=pylab.contourf(gainmap.clip(max=vmax),levels,cmap=pylab.get_cmap(self.cmap.get()))
            #C1=pylab.imshow(gainmap.astype('int32'),aspect='auto',cmap=pylab.get_cmap(self.cmap.get()),vmin=vmin,vmax=vmax)
            pylab.colorbar(C1,cax=self.cax)
            if seg=='A':
                self.Againmap=gainmap
            elif seg=='B':
                self.Bgainmap=gainmap
        if show_which=='Summed Exposure':
            if self.Aexp==None: self.Aexp=pyfits.getdata(os.path.join( data_dir, '12676-all-A-2.fits' ), 1)
            if self.Bexp==None: self.Bexp=pyfits.getdata(os.path.join( data_dir, '12676-LV-B-2.fits'), 1)
            vmin=int(self.vmin.get())
            vmax=int(self.vmax.get())
            if (vmin==0 and vmax==0):
                vmin=0
                vmax=400
            if vmin<0: vmin=0
            if vmax<100: vmax=400
            self.vmin.set(vmin)
            self.vmax.set(vmax)
            levels=[25,50,100,200,250,300,350,400]
            if seg=='A':
                exp=self.Aexp
            elif seg=='B':
                exp=self.Bexp
            C1=pylab.imshow(exp,aspect='auto',vmin=vmin,vmax=vmax,cmap=pylab.get_cmap(self.cmap.get()))
            #C1=pylab.contourf(Aexp,levels)
            pylab.colorbar(C1,cax=self.cax)
        if show_which=='Summed Dark':
            if self.Adark==None: self.Adark=pyfits.getdata('all_pha_0.fits',1)
            if self.Bdark==None: self.Bdark=pyfits.getdata('all_pha_0.fits',2)
            vmin=int(self.vmin.get())
            vmax=int(self.vmax.get())
            if (vmin==0 and vmax==0):
                vmin=0
                vmax=4
            if vmin<0: vmin=0
            if vmax>13: vmax=4
            self.vmin.set(vmin)
            self.vmax.set(vmax)
            levels=[25,50,100,200,250,300,350,400]
            if seg=='A':
                dark=self.Adark
            if seg=='B':
                dark=self.Bdark
	    #dark = blur_image(dark,6,10)
            C1=pylab.imshow(dark,aspect='auto',vmin=vmin,vmax=vmax,cmap=pylab.get_cmap(self.cmap.get()))
            #C1=pylab.contourf(Aexp,levels)
            pylab.colorbar(C1,cax=self.cax)
        if self.extract.get()!='None':
            text=self.extract.get()
            offset=int(self.extract_offset.get())
            SEGMENT=text[:4]
            OPT_ELEM=text[5:10]
            CENWAVE=int(text[11:15])
            APERTURE=text[16:]
            XTRACTAB=pyfits.getdata('u8k1433nl_1dx.fits',1)
            index=numpy.where((XTRACTAB['SEGMENT']==SEGMENT) & (XTRACTAB['OPT_ELEM']==OPT_ELEM) & (XTRACTAB['CENWAVE']==CENWAVE) &(XTRACTAB['APERTURE']==APERTURE))
            B_SPEC=XTRACTAB[index]['B_SPEC']+offset
            HEIGHT=XTRACTAB[index]['HEIGHT']
            B_BKG1=XTRACTAB[index]['B_BKG1']+offset
            B_BKG2=XTRACTAB[index]['B_BKG2']+offset
            BHEIGHT=XTRACTAB[index]['BHEIGHT']
            
            pylab.axhspan(B_SPEC-(HEIGHT-1)/2,B_SPEC+(HEIGHT-1)/2,facecolor='0.5',alpha=0.4)
            pylab.axhspan(B_BKG1-(BHEIGHT-1)/2,B_BKG1+(BHEIGHT-1)/2,facecolor='0.2',alpha=0.4)
            pylab.axhspan(B_BKG2-(BHEIGHT-1)/2,B_BKG2+(BHEIGHT-1)/2,facecolor='0.2',alpha=0.4)

        if self.show_spectrum.get():
            sample=pyfits.open('lbp102neq_corrtag_a.fits')
            im=ttag_image(sample[1].data)
            #pylab.imshow(im,cmap=pylab.get_cmap('gist_yarg'),aspect='auto',interpolation='nearest',alpha=.7)

        if self.xmin.get()!='ALL' and self.xmax.get()!='ALL':
            pylab.xlim(int(self.xmin.get()),int(self.xmax.get()))
        if self.ymin.get()!='ALL' and self.ymax.get()!='ALL':
            pylab.ylim(int(self.ymin.get()),int(self.ymax.get()))
        self.canvas.show()

    def exit(self):
        self.myContainer.quit()
        sys.exit(1)

    def extract_dqs(self):
        pylab.figure(1)
        pylab.subplot(1,1,1)
        xlim_a=(1200,15099)
        ylim_a=(335,699)
        xlim_b=(950,15049)
        ylim_b=(400,749)
        #seg=self.segment.get()
        seg=self.segment2.get()
        if seg=='A':
            xlim=xlim_a
            ylim=ylim_a
        elif seg=='B':
            xlim=xlim_b
            ylim=ylim_b
        dq_array=numpy.zeros((1024,16384),int)
        if self.draw_DQ.get():
            seg=self.segment2.get()
            bpix=pyfits.open('bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg:
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']-1  #width, not delta
                    dy=line['DY']-1
                    dq=line['DQ']
                    subarray=dq_array[ly:ly+dy,lx:lx+dx]
                    index=numpy.where(subarray!=dq)
                    dq_array[ly:ly+dy,lx:lx+dx][index]+=dq
        if self.draw_newDQ.get():
            bpix=pyfits.open('new_bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg:
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']-1 #width, not delta
                    dy=line['DY']-1
                    dq=line['DQ']
                    subarray=dq_array[ly:ly+dy,lx:lx+dx]
                    index=numpy.where(subarray!=dq)
                    dq_array[ly:ly+dy,lx:lx+dx][index]+=dq
        #pylab.contourf(dq_array,aspect='auto',levels=[0,2,4,8,16,32,64,128,256])
        sdqflags=int(self.SDQFLAG.get())
        height=int(self.height.get())
        width=int(self.width.get())
        extracted=dq_array[height-width/2:height+width/2,xlim[0]:xlim[1]]
        index=numpy.where(sdqflags&extracted)
        xs=[xlim[0],xlim[1],xlim[1],xlim[0],xlim[0]]
        ys=[height-width/2,height-width/2,height+width/2,height+width/2,height-width/2]
        pylab.plot(xs,ys,'b--',lw=3)
        total=extracted.shape[0]*extracted.shape[1]
        self.affected_columns.set(100*len(set(index[1]))/float(xlim[1]-xlim[0]))
        self.affected_pixels.set(100*len(index[0])/float(total))
        self.canvas.show()
        self.histogram()

    def find_best(self):
        pylab.figure(1)
        pylab.subplot(1,1,1)
        number_affected=[]
        index_affected=[]
        hs=[]
        sdqflags=int(self.SDQFLAG.get())
        height=int(self.height.get())
        width=int(self.width.get())
        xlim_a=(1200,15099)
        ylim_a=(335,699)
        xlim_b=(950,15049)
        ylim_b=(400,749)
        seg=self.segment2.get()
        if seg=='A':
            xlim=xlim_a
            ylim=ylim_a
        elif seg=='B':
            xlim=xlim_b
            ylim=ylim_b
        dq_array=numpy.zeros((1024,16384),int)
        if self.draw_DQ.get():
            seg=self.segment2.get()
            bpix=pyfits.open('bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg:
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']-1   #width, not delta
                    dy=line['DY']-1   #width, not delta
                    dq=line['DQ']
                    subarray=dq_array[ly:ly+dy,lx:lx+dx]
                    index=numpy.where(subarray!=dq)
                    dq_array[ly:ly+dy,lx:lx+dx][index]+=dq
        if self.draw_newDQ.get():
            bpix=pyfits.open('new_bpix.fits')
            for line in bpix[1].data:
                if line[0]=='FUV'+seg:
                    lx=line['LX']
                    ly=line['LY']
                    dx=line['DX']-1
                    dy=line['DY']-1
                    dq=line['DQ']
                    subarray=dq_array[ly:ly+dy,lx:lx+dx]
                    index=numpy.where(subarray!=dq)
                    dq_array[ly:ly+dy,lx:lx+dx][index]+=dq
        sweep=range(height-1*width,height+1*width)
        for h in sweep:
            if (h-1*width>ylim[0]) & (h+1*width<ylim[1]):
                extracted=dq_array[h-width/2:h+width/2,xlim[0]:xlim[1]]
                index=numpy.where(sdqflags&extracted)
                hs.append(h)
                number_affected.append(len(index[1]))
        number_affected=numpy.array(number_affected)
        hs=numpy.array(hs)
        min_aff=number_affected[number_affected.argmin()]
        min_indexes=numpy.where(number_affected==min_aff)
        hs=hs[min_indexes]
        number_affected=number_affected[min_indexes]
        closest_index=(numpy.fabs(hs-height)).argmin()
        total=extracted.shape[0]*extracted.shape[1]
        self.affected_pixels.set(100*(min_aff/float(total)))
        height=hs[closest_index]
        xs=[xlim[0],xlim[1],xlim[1],xlim[0],xlim[0]]
        ys=[height-width/2,height-width/2,height+width/2,height+width/2,height-width/2]
        pylab.plot(xs,ys,'-.',lw=6)
        extracted=dq_array[height-width/2:height+width/2,xlim[0]:xlim[1]]
        index=numpy.where(sdqflags&extracted)
        self.height.set(height)
        self.affected_columns.set(100*len(set(index[1]))/float(xlim[1]-xlim[0]))
        self.canvas.show()
        self.histogram(height,width)
        
    def find_low(self,thresh=.5):
        which='dark'
        seg=self.segment.get()
        xlim_a=(1176,15219)
        ylim_a=(414,558)
        xlim_b=(870,14984)
        ylim_b=(470,608)
        ylen=3
        xlen=5
        if seg=='A':
            infile=os.path.join( data_dir, '12676-all-A-2.fits' )
            lx=xlim_a[0]
            ux=xlim_a[1]
            ly=ylim_a[0]
            uy=ylim_a[1]
        elif seg=='B':
            infile = os.path.join( data, '12676-LV-B-2.fits' )
            lx=xlim_b[0]
            ux=xlim_b[1]
            ly=ylim_b[0]
            uy=ylim_b[1]
        if which=='dark':
            thresh=2.0
            ylen=21
            xlen=13
            if seg=='A':
                dark=pyfits.getdata('all_pha_0.fits',1) 
                lx=1172
                ly=340
                ux=15150
                uy=700
            elif seg=='B':
                dark=pyfits.getdata('all_pha_0.fits',2)
                lx=1040
                ly=402
                ux=14945
                uy=754
            print 'Blurring array'
            array=blur_image(dark,10,6)
            #array=dark
        else:
            #print 'Blurring array'
            #array=blur_image(pyfits.getdata(infile,1),10,6)
            array=pyfits.getdata(infile,1)
        coords=[]
        for ystep in range(ly,uy-ylen):
            print ystep
            for xstep in range(lx,ux-xlen):
                sub_array=array[ystep:ystep+ylen,xstep:xstep+xlen]
                med=numpy.median(sub_array)
                std=sub_array.std()
                if which=='dark':
                    index=numpy.where(sub_array>((med)*thresh))
                else:
                    index=numpy.where(sub_array<((med)*thresh))
                if len(index[0])>0:
                    for i in range(len(index[0])):
                        coords.append((index[0][i]+ystep,index[1][i]+xstep))
        coord_set=set(coords)
	'''
        final_coords=[]
        length=len(coord_set)
        for i,item in enumerate(coord_set):
            print i,'/',length
            if (coords.count(item)>=30) and (self.has_neighbor(item,coord_set)):
                print 'True'
                final_coords.append(item)
	'''
	final_coords = list( coord_set )
        if len(final_coords)>0:
            xs=[line[1] for line in final_coords]
            ys=[line[0] for line in final_coords]
            '''
            grid=numpy.zeros((1024,16384))
            for y in ys:
                for x in xs:
                    grid[y,x]=1
            y_cont=[]
            x_cont=[]
            for y in range(1024):
                for x in range(16384):
                    if len(y_cont)==0 & grid[y,x]==1:
                        y_cont.append(y)
                        x_cont.append(x)
                        for i in range(len(60)):
                            for j in range(len(60)):
                                X=x+i
                                Y=y+j
                                if grid[Y,X]==1 & ((X-1 in x_cont) | (X+1 in x_cont) | (Y-1 in y_cont) | (Y+1 in y_cont)):
                        y_cont.append(y)
                        x_cont.append(x)
            '''
	    
            pylab.scatter(xs,ys)
            self.canvas.show()

	    data_array = numpy.zeros((1024,16384))
	    for x,y in zip(xs,ys):
		    data_array[y,x] = 1

	    regions = find_contiguous(data_array)
	    for region in regions:
		    lx = region[0]
		    dx = region[1]
		    ly = region[2]
		    dy = region[3]

                    x=[lx,lx+dx,lx+dx,lx,lx]
                    y=[ly,ly,ly+dy,ly+dy,ly]
		    print lx,ly,dx,dy
                    pylab.plot(x,y,'-',lw=4,color='r',linestyle='-')
	    print 'Found all regions'
	    self.canvas.show()

    def has_neighbor(self,point,coord_list):
        neighbor=False
        for y in (-1,0,1):
            for x in (-1,0,1):
                if not (x==0 and y==0):
                    if (point[0]+y,point[1]+x) in coord_list:
                        neighbor=True
        return neighbor

    def help_file(self):
        self.helpwindow=Toplevel()
        self.frame3=Frame(master=self.helpwindow)
        self.helpwindow.title('Help')
        self.closebutton3 = Button(self.helpwindow, text="CLOSE", fg="red", command=self.helpwindow.destroy)
        self.closebutton3.grid(row = 12, column = 0)

        self.text=ScrolledText(self.helpwindow)
        self.text.grid(row=0,column=0,sticky=W+E+N+S)
        help_file=open('help.txt')
        help_string=help_file.readlines()
        self.text.importfile('help.txt')
        #for line in help_file.readlines():

    def histogram(self,height=-99,width=-99):
        xlim_a=(1183,15230)
        ylim_a=(301,750)
        xlim_b=(993,15015)
        ylim_b=(371,771)
        caption=False
        seg=self.segment2.get()
        if seg=='A':
            xlim=xlim_a
            ylim=ylim_a
            gainmap=self.Againmap
        if seg=='B':
            xlim=xlim_b
            ylim=ylim_b
            gainmap=self.Bgainmap
        if height==-99: height=int(self.height.get())
        if width==-99: width=int(self.width.get())
        gain_1d=(gainmap[height-(width-1)/2:height+(width-1)/2,xlim[0]:xlim[1]].flatten())
        #index=numpy.where(gain_1d>0)
        pylab.figure(2)
        pylab.clf()
        if caption:
            ax2=self.fig2.add_axes((.15,.37,.8,.5))
        else:
            ax2=self.fig2.add_axes((.15,.15,.8,.7))
        pylab.cla()
        pylab.suptitle('Modal Gain')
        pylab.title('Segment= %s, Y=%s, dy=%s, years=%s, xshift=%s' %(self.segment2.get(),self.height.get(),self.width.get(),self.degrade_years.get(),self.degrade_xshift.get()),fontsize=12)
        #ax2.hist(gain_1d[index]*100,normed=True,bins=range(20),align='mid',color='r',alpha=.3)
        ax2.hist(gain_1d,normed=True,bins=range(20),align='mid',color='b')
        N_lte_3=len(numpy.where(gain_1d<=3)[0])
        N_tot=len(gain_1d)
        ys=pylab.ylim()
        pylab.text(1,ys[1]*.9,'Pixels below gain=3:\n%d/%d' %(N_lte_3,N_tot),fontsize=12)
        try:
            ax2.locator_params(nbins=20,tight=True,axis='x')
        except:
            pass
        pylab.xlabel('Gain')
        pylab.ylabel('Normalized Counts')
        t='''Modal Gain in extracted region of Segment %s.
At Y=%s, dy=%s and SDQFLAGS=%s, %1.3f percent of
pixels and %1.3f percent of colums will be flagged.''' %(self.segment2.get(),self.height.get(),self.width.get(),self.SDQFLAG.get(),float(self.affected_pixels.get()),float(self.affected_columns.get()))
        #t='''Modal Gain in extracted region of Segment %s.
#At Y=%s, dy=%s and SDQFLAGS=%s, %1.3f percent of
#pixels and %1.3f percent of colums will be flagged.''' %(self.segment2.get(),self.height.get(),self.width.get(),self.SDQFLAG.get(),float(self.affected_pixels.get()),float(self.affected_columns.get()))
        if caption:
            self.fig2.text(.1,.1,t,fontsize=12)
        self.canvas2.show()

    def open_file(self):
        filename = tkFileDialog.askopenfilename(filetypes=[("allfiles","*"),("pythonfiles","*.py")])
        print filename
	self.data_file = filename


    def save_current(self):
        tkFileDialog.asksaveasfile()

    def save_fig(self):
        pylab.figure(2)
        out_file=tkFileDialog.asksaveasfilename()
        pylab.savefig(out_file)

    def set_cmap(self,cmap):
        if cmap=='autumn': pylab.autumn()
        elif cmap=='spring': pylab.spring()
        self.canvas.show()

    def show_limit_boxes(self):
        self.vmin_box.grid_remove()
        self.vmax_box.grid_remove()
        self.vmin_box_label.grid_remove()
        self.vmax_box_label.grid_remove()
        self.ymin_box.grid_remove()
        self.ymax_box.grid_remove()
        self.ymin_box_label.grid_remove()
        self.ymax_box_label.grid_remove()
        self.xmin_box.grid_remove()
        self.xmax_box.grid_remove()
        self.xmin_box_label.grid_remove()
        self.xmax_box_label.grid_remove()

        self.degrade_years_entry.grid_remove()
        self.degrade_loc_entry.grid_remove()
        self.toggle_degrade.grid_remove()

        self.toggle_dq.grid_remove()
        self.toggle_mydq.grid_remove()
        self.dq_show_label.grid_remove()
        self.dq_show_box.grid_remove()

        which=int(self.grid_limits.get())
        if which==1:
            self.vmin_box.grid()
            self.vmax_box.grid()
            self.vmin_box_label.grid()
            self.vmax_box_label.grid()
            self.ymin_box.grid()
            self.ymax_box.grid()
            self.ymin_box_label.grid()
            self.ymax_box_label.grid()
            self.xmin_box.grid()
            self.xmax_box.grid()
            self.xmin_box_label.grid()
            self.xmax_box_label.grid()
        elif which==2:
            self.degrade_years_entry.grid()
            self.degrade_loc_entry.grid()
            self.toggle_degrade.grid()
        elif which==3:
            self.toggle_dq.grid_remove()
            self.toggle_mydq.grid_remove()
            self.dq_show_label.grid_remove()
            self.dq_show_box.grid_remove()

    def update_plot(self,event):
        pylab.figure(1)
        if self.xmin.get()!='ALL' and self.xmax.get()!='ALL':
            pylab.xlim(int(self.xmin.get()),int(self.xmax.get()))
        else:
            pylab.xlim(0,16384)
        if self.ymin.get()!='ALL' and self.ymax.get()!='ALL':
            pylab.ylim(int(self.ymin.get()),int(self.ymax.get()))
        else:
            pylab.ylim(0,1024)
        self.canvas.show()

def main():
	root = Tk()
	app = App(root)
	root.mainloop()
