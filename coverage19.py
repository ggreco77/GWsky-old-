# -*- coding: utf-8 -*-

from __future__ import print_function

try:
   import cPickle as pickle
except:
   import pickle

import fileinput
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.dates as mdates

import dateutil

# Python 3 support
try:
    from Tkinter import *
    import tkMessageBox
    import tkFont
except ImportError:
    from tkinter import *
    from tkinter import font, messagebox

from math import cos, sin, acos, asin, atan, degrees, radians, pi

import astropy
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table

import numpy as np

from mocpy import MOC

from aladinSAMP import AladinScriptCommands 
aladin = AladinScriptCommands()

aladin.setconf_icrsd() # setting ICRSd

from config_values import UserValues
from lvc_skymap import LVCskymap
from query import Query
from airmass import Airmass
from moon import Moon
from moc_region import MOC_confidence_region

# creating folders in Aladin planes
from utils import Utils
Utils.create_folders(folders=["Queries", "Coords", "FoV"])
Utils.load_user_fov("GWsky_fov.vot")

# init. votable
table = Table()
table['a'] = [0, 0]
table.write("GWsky_query_items", format =
            'votable', overwrite = True)



class SkyCoverage(object): 
    """Moving the FoV-footprint coverage by choosing a starting pointing."""
    
    SHIFT_CORRECTION = 0.00001  # A shift correction of 0.00001 is added
                                 # --> to escape math error during the FoV sequence

    def __init__(self, infile_coords='GWsky_coords'):
        """Creating a class in which the instance attributes are based on the dictionary
       "GWsky_coords" by default. GWsky_coords is created and pickled by the "config_values"
       module and will be deleted when the "SkyCoverageGUI" will be closed.
       It contains the keys: "ra", "dec". ra and dec represents the central location of a FoV. 
        
        Starting sky coordinates:
             self.input_ra: right ascension [deg]
             self.input_dec: declination [deg]
         """
       
        self.infile_coords = infile_coords
        self.entries_GWsky_new =[] # new entries during the FoV sequence
        
        self.user = UserValues() # compositions
        self.lvc = LVCskymap()
        self.query = Query()
        self.airmass = Airmass()
        self.moon = Moon()
        
        with open(infile_coords, 'rb') as data:  
            coords_GWsky = pickle.load(data)
            
        for k, v in coords_GWsky.items():          
            setattr(self, k, v)
            
            
    def ra0ra1(self, A, dec0, dec1):
        """From the angular distance:
           cos(A) = sin(Dec1)sin(Dec2)+cos(Dec1)cos(Dec2)cos(ra1-ra2) --> 
           cos(ra1-ra2) = [cos(A)-sin(dec0)sin(dec1)]/[cos(dec0)cos(dec1)]."""

        dec0, dec1, A = radians(dec0),  radians(dec1), radians(A)
        cos_ra0_ra1 = ( cos(A)-sin(dec0)*sin(dec1) )/( cos(dec0)*cos(dec1) )
        ra0ra1 = degrees( acos(cos_ra0_ra1) )

        return  round(ra0ra1, 5)         
    
    def __updating_center_coords(self, ra, dec):
        """Getting/Updating FoV-center (ra, dec) in the dict named by default "GWsky_coords".
.          For each tile across the sky the file is updated."""""

        with open('GWsky_coords', 'rb') as data:
            coords_GWsky = pickle.load(data)
            
        coords_GWsky['input_ra'], coords_GWsky ['input_dec'] = ra, dec

        with open('GWsky_coords', 'wb') as data:
            pickle.dump(coords_GWsky, data)

    def __fov_stats(self, ra, dec, moon_illumination, sep_fov_moon, table,
                    r, dp_dr, airmass_values, datestrings, prob_fov):    
        """Managing the descriptive statistic window."""
            
        self.__updating_center_coords(ra,dec) # updating file to read value
        
        fov_statistics = FoVstatistics()
        fov_statistics.plot_stats(ra,dec,moon_illumination, sep_fov_moon, table,
                                  r, dp_dr,airmass_values, datestrings, prob_fov)
        
    def update_pointings_file(self, infile, ra, dec, prob_fov, skymap):
         """The central location (ra[deg], dec[deg]) and the integrated probability of
             a selected FoV are saved locally in an external file.
             By default the file is named "GWsky_pointings.txt"."""
           
         with open(infile, 'a') as pointing:
             pointing.write(str(ra) + ' ' + str(dec)+ ' ' + str(prob_fov) + ' ' + skymap +'\n')

    def __query_shape(self, ra, dec, fov_shape):
        """Return the catalog query according with the defined-user FoV shape:
                   (1) box and (2) circle. """
        
        if self.user.get_fov_shape() != 2:  # box FoV
                    query_result = self.query.query_box(
                       ra, dec, self.user.get_fov_width(), self.user.get_fov_height(), self.user.get_catalog(),
                       self.user.get_column_1(), self.user.get_filter_1())              
        else: # circle FoV
            query_result = self.query.query_circle(
               ra, dec, self.user.get_fov_radius(), self.user.get_catalog(),
               self.user.get_column_1(), self.user.get_filter_1())
            
        return query_result

    def __prob_shape(self, ra, dec, fov_shape):
        """Return the integrated probability according with the defined-user FoV shape:
                   (1) box and (2) circle."""
        
        if self.user.get_fov_shape() !=2:  # box FoV
            prob_fov = self.lvc.prob_in_box(
               ra, dec, self.user.get_fov_width(), self.user.get_fov_height())                
        else: # circle FoV
            prob_fov = self.lvc.prob_in_circle(
               ra, dec, self.user.get_fov_radius())
            
        return  prob_fov
                  
    def pick_coverage(self, ra, dec):
        """Setting GWsky: with statistic window (A); without statistic window (D)."""                                      

        if self.user.get_GWsky_basic() == "A":  # full version 
            
            prob_fov = self.__prob_shape(ra, dec,
                                         self.user.get_fov_shape()) # integrated prob
            
            moon_illumination =  self.moon.illumination()           # moon_illumination
            
            sep_fov_moon = self.moon.from_fov(ra, dec)*u.deg        # moon_dist
            sep_fov_moon = sep_fov_moon.round(1)

            query_result = self.__query_shape(ra, dec,
                                              self.user.get_fov_shape())  # query

            airmass_values, datestrings = self.airmass.airmass_step(ra, dec)
            
            # TEST------------------------------------------------------#    
            #fov_sep = Utils.separation(self.input_ra, self.input_dec,   #
            #                            ra, dec)                        #
            #print ('The distance between 2 consecutive FoV centers is', #
            #       fov_sep.round(6))                                    #
            #-----------------------------------------------------------#
            
            r, dp_dr = self.lvc.conditional_distance_fov_center(ra, dec)
            
            self.__fov_stats(ra, dec, moon_illumination, sep_fov_moon, query_result,
                             r, dp_dr, airmass_values, datestrings, prob_fov)          #  Stats win

            #print airmass_values, datestrings
            
        elif self.user.get_GWsky_basic() == "D":  # basic version-> no Stats win
            prob_fov = self.__prob_shape(ra, dec,
                                         self.user.get_fov_shape()) # query

            # TEST------------------------------------------------------#  
            #fov_sep = Utils.separation(self.input_ra, self.input_dec,  #
            #                            ra, dec)                       #
            #print ('The distance between 2 consecutive FoV centers is',#
            #       fov_sep.round(6))                                   #
            #-----------------------------------------------------------#
            
            self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov, skymap=self.user.get_skymap())
            
    def intercardinal_distance(self, ra, dec, shift_up_down, shift_right_left):
        """Moving from the fixed cardinal direction using the bi-directional windows;
           shift_up_down ↕ and/or shift_right_left ↔."""

        if shift_right_left > 0:
           shift_east_west = self.ra0ra1((shift_right_left-self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra + shift_east_west 
         
        elif shift_right_left < 0 :
           shift_east_west = self.ra0ra1((shift_right_left + self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra - shift_east_west
         
        else:
           dist = ra

        return dist 

    def load_entries(self, infile_entries):
        """Opening the file in which the input entries are stored: ra_1 dec_1 ra_2 dec_2 ... ra_n dec_n
           By default the file is named "GWsky_entries". "GWsky_entries" is created from the
           "show_starting_fov" method in "StartingFoV" class.
           A error message invites users to press the "Start FoV" button before carrying out any"""
       
        try:
            with open(infile_entries, 'rb') as data:
                entries_GWsky = pickle.load(data)
                return entries_GWsky
        except IOError as io_error:
            message = "Press the button 'Start FoV' and choose a starting position; \
                        by default the sky position of the max probability pixel is shown"
            tkMessageBox.showerror ('Error', message)
                         
    def north(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in North direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔ """

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            north_pointing = [(dist),
                               (float(dec_start) + self.user.get_fov_height() + shift_up_down)]
             
            ra, dec = north_pointing[0], north_pointing[1]
            
            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Latitude angle(s) must be within -90 deg <= angle <=90 deg, got' + ' ' + str(dec) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def south(self, shift_up_down=0, shift_right_left=0):    
        """Moving the FoV tiles in South direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔"""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            south_pointing = [(dist), (float(dec_start) - self.user.get_fov_height() - shift_up_down)]
                    
            ra, dec = south_pointing[0], south_pointing[1]

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Latitude angle(s) must be within -90 deg <= angle <=90 deg, got' + ' ' + str(dec) + ' ' + 'deg'
            
            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
                                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)    
          
    def east(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in East direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
           A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):              
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                        float(dec_start), float(dec_start))
                
            east_pointing = [(float(ra_start) + ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = east_pointing[0], east_pointing[1]

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Longitude angle(s) must be within 0 deg <= angle <=360 deg, got' + ' ' + str(ra) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)           

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
                   
    def west(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in West direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
            A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  
      
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
               
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                      float(dec_start), float(dec_start))

            west_pointing = [(float(ra_start) - ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = west_pointing[0], west_pointing[1]

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Longitude angle(s) must be within 0 deg <= angle <=360 deg, got' + ' ' + str(ra) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
        

class SkyCoverageGUI(Toplevel):
    """
        The GUI consists of 9 buttons; 4  directional buttons to move the FoVs
        in cardinal directions (North, South, East, West) 4  buttons to shift the FoVs
        from a consecutive cardinal direction (↕, ↕, ↔, ↔) and 1 button to get a
        new FoV position (Start FoV). ***Input values in deg***.
    """
    
    def __init__(self, tkMainWin):
        frame = Frame(tkMainWin, border=9, bg="dim grey")
        frame.pack()
        
        self.B02 = Button(frame,text="N", )   
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.grid(row=0, column=2)
        
        self.B12 = Button(frame,text="↕↔",
                          fg="grey")  
        self.B12.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B12.grid(row=1, column=2)
  
        self.B30 = Button(frame, text="E", )  
        self.B30.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B30.grid(row=3,column=0)

        self.B31 = Button(frame,text="↕↔",
                          fg="grey")   
        self.B31.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B31.grid(row=3, column=1)

        self.B32 = Button(frame, text="Start FoV",
                          fg="red4",)
        self.B32.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B32.grid(row=3, column=2)

        self.B33 = Button(frame,text="↕↔",
                          fg="grey")  
        self.B33.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B33.grid(row=3, column=3)

        self.B34 = Button(frame, text="W",) 
        self.B34.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B34.grid(row=3, column=4)

        self.B42 = Button(frame,text="↕↔",
                          fg ="grey") 
        self.B42.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B42.grid(row=4, column=2)
        
        self.B52 = Button(frame,text="S", ) 
        self.B52.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B52.grid(row=5, column=2)

        ## Adjustments Btns
        self.B60 = Button(frame,text="↞",
                          fg ="grey", pady=3) 
        self.B60.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B60.grid(row=6, column=0)

        self.B61 = Button(frame,text="↠",
                          fg ="grey", pady=3,) 
        self.B61.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B61.grid(row=6, column=1)

        self.B62 = Button(frame,text="✓ Accept",fg ="green4",
                          pady=1, padx=11) 
        self.B62.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B62.grid(row=6, column=2)

        self.B63 = Button(frame,text="↟",
                          fg ="grey",pady=3) 
        self.B63.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B63.grid(row=6, column=3)

        self.B64 = Button(frame,text="↡",
                          fg ="grey",pady=3) 
        self.B64.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B64.grid(row=6, column=4)

        # Folder
        self.B72 = Button(frame,text="⏩▶ Folder",
                          fg ="gold4",pady=1,padx=10) 
        self.B72.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B72.grid(row=7, column=2)

        # MOC plot default
        self.B82 = Button(frame,text="MOC plot",
                          fg ="steel blue",pady=1, padx=12) 
        self.B82.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B82.grid(row=8, column=2)

        # ObsMOC
        self.B80 = Button(frame,text="⧗",
                          fg ="black",pady=3) 
        self.B80.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B80.grid(row=8, column=0)

        # ObsSkymap
        self.B81 = Button(frame,text="☰",
                          fg ="black",pady=3) 
        self.B81.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B81.grid(row=8, column=1)

        # localize sources
        self.B83 = Button(frame,text="◉",
                          fg ="black",pady=3) 
        self.B83.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B83.grid(row=8, column=3)

        # load new skymap
        self.B84 = Button(frame,text="⬊",
                          fg ="black",pady=3) 
        self.B84.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B84.grid(row=8, column=4)
    
    # Actions    
    def clicked(self, event):
        """Moving the user-defined FoV footprint."""
        
        run_sequence = SkyCoverage()
        
        if event.widget == self.B02:
            run_sequence.north()         # north

        if event.widget == self.B12:
            move_fov = ShiftFoV()
            move_fov.north_shift()       # ↕↔
                      
        if event.widget == self.B30:
            run_sequence.east()          # east

        if event.widget == self.B31:
            move_fov = ShiftFoV()
            move_fov.east_shift()        # ↕↔
            
        if event.widget == self.B32: 
            new_starting_fov = StartingFoV()  # start FoV
            new_starting_fov

        if event.widget == self.B33:
            move_fov = ShiftFoV()
            move_fov.west_shift()        # ↕↔
            
        if event.widget == self.B34:   
            run_sequence.west()          # west
            
        if event.widget == self.B42:
            move_fov = ShiftFoV()
            move_fov.south_shift()       # ↕↔
            
        if event.widget == self.B52:    
            run_sequence.south()         # south

        if event.widget == self.B60:     # ↞        
            adj = Adjustments()
            adj.adj_east()
            
        if event.widget == self.B61:     # ↠       
            adj = Adjustments()
            adj.adj_west()            
                
        if event.widget == self.B63:     # ↟          
            adj = Adjustments()
            adj.adj_north()

        if event.widget == self.B64:     # ↡           
            adj = Adjustments()
            adj.adj_south()         
               
        if event.widget == self.B62:     # ✓ Accept
            adj = Adjustments()
            adj.adj_accept()
            
        if event.widget == self.B72:     # ▶ Folder
            Utils.move_to_folder(planes=['Q:*', 'P:*', 'C_*'],
                                 folders=['Queries', 'FoV', 'Coords'])

        if event.widget == self.B80: # ⧗  ObsMOC
            obs_in_MOC = ObsInMOC()

        if event.widget == self.B81: # ☰  ObsSkymap
            obs_skymap = ObsSkymap()

        if event.widget == self.B82:     # MOC plot
            self.moc = MOC_confidence_region()
            self.user = UserValues()
            self.moc.contour_default(_from=10, _to=100, _step=10, skymap=self.user.get_skymap())

        if event.widget == self.B83:   # ◉ localize sources
            localize = LocalizeSources()
            
        if event.widget == self.B84:   # ⬊ loading a new skymap
            self.load = LoadSkymap()


class LocalizeSources(Toplevel):
    """Localize an astrophysical source(s) within a skymap."""
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()
        self.moc = MOC_confidence_region()
        self.lvc = LVCskymap()
        
        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()

        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Localize Sources in probability skymap")
        self.attributes("-topmost", True)

        # label 1
        self.label_1 = Label(self, text="Is/Are the [ID Source(s) RA (°) DEC (°)]", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: sky coords of maximum probability pixel
        #fov_coords = str(self.user.get_ra_max_pixel()), str(self.user.get_dec_max_pixel()) 
        help_string = ("source_1 079.91935 +43.34139; source_2 063.87703 +33.36612 ")
        
        help_string_default = StringVar(self, value=help_string) 
        self.entry_in = Entry(self, width=30, justify=CENTER,
                             textvariable=help_string_default)

        self.entry_in.grid(row=0, padx=15, column=1)

        # label 2
        self.label_2 = Label(self, text="whithin the", bg="slate grey")
        self.label_2.grid(row=0, column=3, sticky=E, pady=0)

        moc_value = 90  # default     
        moc_default = StringVar(self, value=moc_value)
        
        self.entry_percentage = Entry(self, width=5, justify=CENTER,
                             textvariable=moc_default)
        self.entry_percentage.grid(row=0, padx=2, column=5)

        # label 3
        self.label_2 = Label(self, text="% MOC?  ",
                             bg="slate grey")
        self.label_2.grid(row=0, column=6, sticky=E, pady=0)

        # label 4
        folder = "transients"  # default     
        folder_default = StringVar(self, value=folder)
        
        self.entry_folder = Entry(self, width=15, justify=CENTER,
                             textvariable=folder_default)
        self.entry_folder.grid(row=1, padx=2, column=0)

        # label 3.1
        self.label_3 = Label(self, text="  Folder: ",justify=LEFT,
                             bg="slate grey")
        self.label_3.grid(row=1, column=0, sticky=W, pady=0)

        #from gui_utils import scrollHandler
        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_in['xscrollcommand'] = self.entryScroll.set

        #Btns
        self.show = Button(self, text='Ask',
                           command=self.in_skymap)
        self.show.grid(column=7, row=0, sticky=W, padx=2, pady=5)
        
        self.checkbox = Button(self, text="Pinpoint", fg='black',     
                               command=self.pinpoint)
        self.checkbox.grid(column=8,row=0, sticky=E, padx=2, pady=5)

        self.checkbox = Button(self, text="Dist", fg='black',     
                               command=self.cond_distance_source_for)
        self.checkbox.grid(column=15,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=25,row=0, sticky=W, padx=2, pady=5)


    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_in.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_in.xview_moveto(howMany)
            
    def split_entries_3(self):
        """Splitting the entries in ra and dec; # odd: ra and # even: dec."""
        
        entry_sources = self.entry_in.get().replace(';',' ').replace(',',' ').split()
        
        label = entry_sources[::3]
        source_ra = entry_sources[1::3]
        source_dec = entry_sources[2::3]

        print (label[0::], source_ra[0::], source_dec[0::])
        return label[0::], source_ra[0::], source_dec[0::]
        
    def in_skymap(self):
        """Checking if an object falls in a given probability level."""

        aladin.md(self.entry_folder.get())
        aladin.remove(self.entry_folder.get() + '~1') # removing multiple copy of the folder
        
        percentage = float(self.entry_percentage.get())/100.0
        
        skymap=self.user.get_skymap()
        
        import healpy as hp
               
        labels, ra_transients, dec_transients = self.split_entries_3()

        print (labels, ra_transients, dec_transients)
        
        prob = self.moc.read_prob(self.user.get_skymap())
      
        for ra_transient, dec_transient, label in zip(
           ra_transients, dec_transients, labels):

            aladin.draw_newtool(label)
            aladin.draw_source(ra_transient, dec_transient, label)
            aladin.mv(label, self.entry_folder.get())
           
            ipixs = self.moc.ipixs_in_percentage(prob, percentage )
            self.nside = int(self.user.get_nside())
              
            theta = 0.5 * np.pi - np.deg2rad(float(dec_transient))
            phi = np.deg2rad(float(ra_transient))
            ipix = hp.ang2pix(self.nside, theta, phi)
        
            is_there = ipix in ipixs

            if is_there == True:
                res_true = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + " " + "dec="+str(dec_transient)+"°" + " " + \
                            "(labels: " + label+")" + " " +  "lies within the" + " " + str(percentage*100)+'%' + " " + "c.l.\n" +"["+skymap+"]")
                
                print (res_true)
                tkMessageBox.showinfo("Localize Result: " + label, res_true)
            else:
                res_false = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + " " +"dec="+str(dec_transient)+"°" + " " + \
                             "(labels: " + label+")" + " " + "is outside the" + " " + str(percentage*100)+'%' + " " + "c.l.\n" + "["+skymap+"]")

                print (res_false)
                tkMessageBox.showinfo("Localize Results: " + label, res_false)


    def pinpoint(self):
        pinpoint_localize = Pinpoint ()


    def cond_distance_source_for(self):
        """"""
        import numpy
        import time
        from matplotlib import pyplot as plt
        
        labels, ra_transients, dec_transients = self.split_entries_3()
        
        plt.ion()
        for label, ra_transient, dec_transient in zip(labels, ra_transients, dec_transients):

           r, dp_dr = self.lvc.conditional_distance_fov_center(float(ra_transient), float(dec_transient))
           #print (r, dp_dr)
           #plt.close()
           #fig = plt.figure()
           #plt.plot(r,dp_dr)
           #ax.set_xlabel('distance (Mpc)')
           #ax.set_ylabel('prob Mpc$^{-1}$')
           #plt.show()
           #_ = raw_input("Press [enter] to continue.") # wait for input from the user
           #time.sleep(2)
           #plt.close(fig)
           #plt.gcf().clear()
           #plt.draw()
           #plt.close(plot1)
           #time.sleep(2)
           #_ = raw_input("Press [enter] to continue.") # wait for input from the user
           #plt.close()
           
           self.cond_distance_source(label, r, dp_dr)


    def cond_distance_source(self, label, r, dp_dr):
        """"""

        #labels, ra_transients, dec_transients = self.split_entries_3()

        #for ra_transient, dec_transient, label in zip(
        #   ra_transients, dec_transients, labels):
        
        #r, dp_dr = self.lvc.conditional_distance_fov_center(float(ra_transient), float(dec_transient))
            #print (label, r,dp_dr)

         

        import matplotlib.pyplot as plt
        #plt.rcdefaults()
        
        

        #plt.figure()
        #plt.plot(r,dp_dr)
        
            #conditional_distance_line_sight = f.add_subplot(221) 
            #conditional_distance_line_sight.plot(r, dp_dr)
            
            #title_string = 'Conditional distance distribution \n along the FoV center'
            #conditional_distance_line_sight.set_title(title_string,fontsize=10)
            #conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            #conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        skymap=self.user.get_skymap()
        #plt.rcdefaults()
        fig, ax = plt.subplots()

        ax.plot(r, dp_dr)
        title_string = label + ':' + ' ' + ' \n conditional distance distribution  along the line of sight \n' + '['+skymap+']'
        ax.set_title(title_string,fontsize=10)
        ax.set_xlabel('distance (Mpc)')
        ax.set_ylabel('prob Mpc$^{-1}$')
            #ax.set_title('How fast do you want to go today?')

        plt.show()

    def close_window(self):
        self.destroy()


class Pinpoint(Toplevel):
    """The class is designed to determine in which level of probability a source is localized."""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")

        #self.localize = LocalizeSources()
        
        
        self.user = UserValues()

        self.moc = MOC_confidence_region()
        
        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()

        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Pinpoint Localization")
        self.attributes("-topmost", True)

        # label 1
        self.label_1 = Label(self, text=" In which level of probability the source(s) falls/fall", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        help_string = ("source_1 079.91935 +43.34139; source_2 063.87703 +33.36612 ")
        
        help_string_default = StringVar(self, value=help_string) 
        self.entry_pin = Entry(self, width=30, justify=CENTER,
                             textvariable=help_string_default)

        self.entry_pin.grid(row=0, padx=15, column=1)

        # label 2
        self.label_2 = Label(self, text=" from  ", bg="slate grey")
        self.label_2.grid(row=0, column=3, sticky=E, pady=0) 

        from_default = StringVar(self, value="10") # default

        self.entry_from = Entry(self, width=5, justify=CENTER,
                                textvariable=from_default)
        self.entry_from.grid(row=0, column=4, columnspan=8, pady=2, sticky='WE')

        # label 3
        self.label_3 = Label(self, text=" to  ", bg="slate grey")
        self.label_3.grid(row=0, column=20, sticky=E, pady=0) 

        to_default = StringVar(self, value="90") # default

        self.entry_to = Entry(self, width=5, justify=CENTER,
                                textvariable=to_default)
        self.entry_to.grid(row=0, column=25, columnspan=8, pady=2, sticky='WE')

        # label 4
        self.label_4 = Label(self, text=" grid  ", bg="slate grey")
        self.label_4.grid(row=0, column=40, sticky=E, pady=0) 

        grid_default = StringVar(self, value="10") # default

        self.entry_grid = Entry(self, width=5, justify=CENTER,
                                textvariable=grid_default)
        self.entry_grid.grid(row=0, column=45, columnspan=8, pady=2, sticky='WE')


        # label 4
        folder = "transients"  # default     
        folder_default = StringVar(self, value=folder)
        
        self.entry_folder = Entry(self, width=20, justify=CENTER,
                             textvariable=folder_default)
        self.entry_folder.grid(row=1, padx=2, column=0)

        # label 3.1
        self.label_3 = Label(self, text="  Folder:",justify=LEFT,
                             bg="slate grey")
        self.label_3.grid(row=1, column=0, sticky=W, pady=0)
        
        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_pin['xscrollcommand'] = self.entryScroll.set

        #Btns
        self.show = Button(self, text='Do',
                           command=self.pinpoint_for)
        self.show.grid(column=55, row=0, sticky=W, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=70,row=0, sticky=W, padx=2, pady=5)

        
    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_pin.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_pin.xview_moveto(howMany)
            
    def split_entries_3(self):
        """Splitting the entries in ra and dec; # odd: ra and # even: dec."""
        
        entry_sources = self.entry_pin.get().replace(';',' ').replace(',',' ').split()
        
        label = entry_sources[::3]
        source_ra = entry_sources[1::3]
        source_dec = entry_sources[2::3]

        print (label[0::], source_ra[0::], source_dec[0::])
        return label[0::], source_ra[0::], source_dec[0::]

    def pinpoint_for(self):
        """"""

        aladin.md(self.entry_folder.get())
        aladin.remove(self.entry_folder.get() + '~1') # removing multiple copy of the folder
        
        labels, ra_transients, dec_transients = self.split_entries_3()
        for ra_transient, dec_transient, label in zip(ra_transients, dec_transients, labels):

            aladin.draw_newtool(label)
            aladin.draw_source(ra_transient, dec_transient, label)
            aladin.mv(label, self.entry_folder.get()) 
           
            self.pinpoint(ra_transient, dec_transient, label)
           
        
            
    def pinpoint(self, ra_transient, dec_transient, label):
            
        """Find in which confidence level the source falls.
        
        Input parameters
        ---------------
        infile : string
            LVC probability sky localization in healpix format
        from_percentage : float
            fractional percentage from 0 to 1
        to_percentage : float
            fractional percentage from 0 to 1
        resolution_percentage : float
            fractional percentage from 0 to 1
        ra, dec : float
            sky coordinates in degrees
        label : string
            name of the object (optional, by default = '')
        show = True
            show the MOC confidence region and the object in the Aladin planes;
        otherwise no (optional; by default = True )
        """
            
        #self.read_skymap(infile)
        #percentage = float(self.entry_percentage.get())/100.0

        prob = self.moc.read_prob(self.user.get_skymap())
        self.nside = int(self.user.get_nside())
        
        from_percentage = float(self.entry_from.get())/100.0
        to_percentage = float(self.entry_to.get())/100.0
        resolution_percentage = float(self.entry_grid.get())/100.0



        import healpy as hp

        #labels, ra_transients, dec_transients = self.split_entries_3()
        #print (labels, ra_transients, dec_transients)



        #for ra_transient, dec_transient, label in zip(ra_transients, dec_transients, labels):
        #    print (ra_transient, dec_transient, label )
        
            # finding ipix
        theta = 0.5 * np.pi - np.deg2rad(float(dec_transient))
        phi = np.deg2rad(float(ra_transient))
        ipix = hp.ang2pix(self.nside, theta, phi)
            

        find="n"
        #print (find)
         
        while from_percentage < to_percentage or find =="y":
            ipixs = self.moc.ipixs_in_percentage(prob, from_percentage)
            is_there = ipix in ipixs

            if is_there != True:                                
                from_percentage = from_percentage + resolution_percentage #adj
                #from_percentage = to_percentage
                #res_n = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + ' ' + "dec="+str(dec_transient)+"°"+" " + "(label:" + label+")" \
                #           + " " + "is not localized within the" + " " + str(from_percentage*100)+'%' + " " + "c.l.")
                
                #tkMessageBox.showinfo("Localize Result: " + label, res_n)
                
            else:
                find="y"      
                res_y = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + ' ' + "dec="+str(dec_transient)+"°"+" " + "(label:" + label+")" \
                           "lies within the" + " " + str(from_percentage*100)+'%' + " " + "c.l.")
                   
                tkMessageBox.showinfo("Localize Result: " + label, res_y)

                return find

        from_percentage = to_percentage
        res_n = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + ' ' + "dec="+str(dec_transient)+"°"+" " + "(label:" + label+")" \
                 + " " + "is not localized within the" + " " + str(from_percentage*100)+'%' + " " + "c.l.")
                
        tkMessageBox.showinfo("Localize Result: " + label, res_n)

                           
    def close_window(self):
        self.destroy()
        

class LoadSkymap(Toplevel):
    """Loading a new skymap."""
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()

        self.wait_visibility()
        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title("Load a new skymap")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="LVC skymap", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: input skymap
        skymap_input = StringVar(self, value=self.user.get_skymap()) 
        self.entry_new_skymap = Entry(self, width=30, justify=CENTER,
                             textvariable=skymap_input)

        self.entry_new_skymap.grid(row=0, padx=15, column=1)

        #Btns
        self.show = Button(self, text='Load',
                           command=self.new_skymap)
        self.show.grid(column=2, row=0, sticky=W, padx=2, pady=5)
        
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=5,row=0, sticky=E, padx=2, pady=5)
    
    def new_skymap(self):
        """Loading a new LVC skymap."""     

        import healpy as hp
        from aladinSAMP import AladinViaSAMP, AladinScriptCommands 
        samp = AladinViaSAMP()

        # update nside
        prob = hp.read_map(self.entry_new_skymap.get(), verbose = False)
        npix = len(prob)
        nside = hp.npix2nside(npix)

        # update coord maximum prob pixel
        ipix_max = np.argmax(prob)
        prob[ipix_max] 
        theta, phi = hp.pix2ang(nside, ipix_max)
        
        ra_max = round(np.rad2deg(phi),5) 
        dec_max = round(np.rad2deg(0.5 * np.pi - theta),5)
        
        with open('GWsky_config', 'rb') as data:
            config_GWsky = pickle.load(data)
            
        config_GWsky['skymap'], config_GWsky['nside'],config_GWsky['ra_max_pixel'],config_GWsky['dec_max_pixel']=\
        self.entry_new_skymap.get(), nside, ra_max, dec_max

        with open('GWsky_config', 'wb') as data:
            pickle.dump(config_GWsky, data)

        samp.send_file(  self.entry_new_skymap.get() )  # sending to Aladin
        aladin.rename(self.entry_new_skymap.get())

        os.remove("GWsky_entries")

    def close_window(self):
        self.destroy()
        
class Adjustments(SkyCoverage):
    """Adjustments FoV position."""

    def __init__ (self):

       SkyCoverage.__init__(self)
       
       self.shift_up = 0.15     # default adjustments   (up)
       self.shift_down = 0.15   #      "               (down)
       self.shift_left = 0.15   #      "               (left)
       self.shift_right = 0.15  #      "               (right)

       
    def adj_north(self):
        """Adjustments FoV position -> north direction"""
            
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_up, shift_right_left=0)
            north_adj = [(dist),
                         (dec_start + 0 + self.shift_up)]
             
            ra, dec = north_adj[0], north_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_south(self):
        """Adjustments FoV position -> south direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_down, shift_right_left=0)
            south_adj = [(dist),
                         (dec_start + 0 - self.shift_down)]
             
            ra, dec = south_adj[0], south_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_east(self):
        """Adjustments FoV position -> east direction"""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
            
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_left),
                                        float(dec_start), float(dec_start))
                          
            east_adj = [(float(ra_start) + ra_distance), (float(dec_start) + 0)]
            ra, dec = round(east_adj[0],5), round(east_adj[1],5)

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def adj_west(self):
        """Adjustments FoV position -> west direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
             
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_right),
                                        float(dec_start), float(dec_start))
            
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            west_adj = [(float(ra_start) - ra_distance), (float(dec_start) + 0)]
            ra, dec = west_adj[0], west_adj[1]

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))

            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_accept(self):
        """Confirming the adjustments FoV position -> open statistic win."""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
                   
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
                    
            self.pick_coverage(float(ra_start), float(dec_start))

            new_sky_pos = [ra_start,dec_start]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    
class StartingFoV(Toplevel):
    """Starting a sequence from a list of FoV(s). The window contains 1 keyboard entries and 3 Buttons.

        entry:
             ra_1 dec_1 ra_2 dec_2 ra_3 dec_3 ... ra_n dec_n [deg]
             By default: sky coords of maximum probability pixel

        Btns:
           Show : draw user-defined FoV footprint(s) in Aladin Plane(s)
           No show : no draw user-defined FoV footprint(s) in Aladin Plane(s)
           Close : close the widget
        """
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()

        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()

        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Starting FoV")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="RA (°) DEC (°)", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: sky coords of maximum probability pixel
        fov_coords = str(self.user.get_ra_max_pixel()), str(self.user.get_dec_max_pixel()) 
        
        max_pixel_default = StringVar(self, value=fov_coords) 
        self.entry_1 = Entry(self, width=30, justify=CENTER,
                             textvariable=max_pixel_default)

        self.entry_1.grid(row=0, padx=15, column=1)

        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_1['xscrollcommand'] = self.entryScroll.set

        #Btns
        self.show = Button(self, text='Show',
                           command=self.show_starting_fov)
        self.show.grid(column=2, row=0, sticky=W, padx=2, pady=5)
        
        self.checkbox = Button(self, text="Not show",      
                               command=self.no_show_starting_fov)
        self.checkbox.grid(column=3,row=0, sticky=E, padx=2, pady=5)

        self.browse_1 = Button(self, text="Browse ...",
                               command=self.browsecsv)
        self.browse_1.grid(row=0, column=4, padx=2, pady=2)

        #self.close = Button(self, text="Obs",  fg='dark green',
        #                    command=self.obs)  
        #self.close.grid(column=4,row=0, sticky=W, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=6,row=0, sticky=E, padx=2, pady=5)


    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry in starting FoV window."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_1.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_1.xview_moveto(howMany)
            
    def __split_entries(self):
        """Splitting the entries in ra and dec; # odd: ra and # even: dec."""
        
        current_fov_coords = self.entry_1.get().replace(';',' ').replace(',',' ').split()
        fov_center_ra = current_fov_coords[0::2]
        fov_center_dec = current_fov_coords[1::2]

        return current_fov_coords, fov_center_ra, fov_center_dec
    
    def show_starting_fov(self):
        """Drawing the FoV footprint(s) in the Aladin plane(s).
         By default: sky coords (ra[deg], dec[deg]) of maximum probability pixel."""   
   
        show_sky_coverage = SkyCoverage()
        
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()
        
        try:
            for ra_starting, dec_starting in zip (fov_center_ra, fov_center_dec):
                aladin.get_FoV(float(ra_starting), float(dec_starting))
                show_sky_coverage.pick_coverage(float(ra_starting), float(dec_starting))
                              
        except ValueError as value_error:
            tkMessageBox.showerror ('Error', value_error)
          
        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky, data)
            
        self.entries_GWsky=[] # re-init.
        
    def no_show_starting_fov(self):
        """No Draw the FoV footprint(s) in the Aladin plane(s);
               useful to re-initialize the sequence."""

        self.entries_GWsky=[] # re-init.     
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()

        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            return pickle.dump(self.entries_GWsky, data)

    def browsecsv(self):
        from tkFileDialog import askopenfilename

        #Tk().withdraw() 
        filename_from_browser = askopenfilename(filetypes=[("fits files","*.fits.gz")])
        #fits = filename.read()
        print (filename_from_browser)

     
    def close_window(self):
        self.destroy()

class ObsSkymap(Toplevel):
    """Initializi"""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()
        self.moc = MOC_confidence_region()

        self.observatory = astropy.coordinates.EarthLocation(
           lat=self.user.get_latitude()*u.deg, lon=self.user.get_longitude()*u.deg,
           height=self.user.get_altitude()*u.m)

        self.obs_time = Time(self.user.get_obs_time())
        self.nside = self.user.get_nside()

        self.wait_visibility()
        
        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title("ObsSkymap --" + "Observability" + " " + "starting from" + \
                   " " + self.user.get_obs_time() + " " + "UT")
        self.attributes("-topmost", True)

        self.bttn_clicks = 0 # counter ">>" Button
        
        # first label
        self.label_1 = Label(self, text="Show the airmass regions over the skymap from 1 to 4 in step of 1",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, padx=0)

        #Btn
        self.show = Button(self, text='Show',
                           command=self.snapshot_airmass)
        self.show.grid(column=4, row=0, sticky=W, padx=2, pady=5)
        
        self.forward = Button(self, text=">>",      
                               command=self.snapshot_airmass_up)
        self.forward.grid(column=5,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=7,row=0, sticky=E, padx=2, pady=5)

    #Actions
    def update_count(self):
        """Return the time in step of 1h when the button ">>" is clicked."""
   
        self.bttn_clicks += 1
              
        dt = TimeDelta(3600.0, format='sec')
        update_time = int(self.bttn_clicks) * dt
        
        obs_time = Time(self.user.get_obs_time())
        time_update = obs_time + update_time
        
        return time_update

    def from_ipixs_to_moc(self, time_input, percentage = 0.9):
        """Return ipix table with the associated airmass"""

        prob = self.moc.read_prob(self.user.get_skymap())
        #percentage = 
        #float(self.entry_percentage.get())/100.0
        
        ipixs = self.moc.ipixs_in_percentage(prob, percentage )
        nside = int(self.user.get_nside())
        
        ra, dec = self.moc.sky_coords(ipixs, nside)
        
        sky_coord = SkyCoord(ra = ra*u.deg,dec=dec*u.deg, frame='icrs')
        altaz = sky_coord.transform_to(AltAz(obstime=time_input, location=self.observatory))
        airmass_values = altaz.secz
        
        contour_ipix = Table([ ra, dec, airmass_values, ipixs ],
                             names = ('RA[deg]', 'DEC[deg]', 'airmass', 'ipix'),
                             meta = {'ipix': 'ipix table'})             # astropy table       

        return contour_ipix


    def snapshot_airmass(self):
        """Showing a snapshot af the airmass at a given time."""
        
        #obs_time = Time(self.user.get_obs_time())
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(self.obs_time)
        
        airmass_start = [1, 2, 3, 4]
        airmass_end = [2, 3, 4, 5.8]
        #colors = ["red", "blue", "ciao","ciao"]

        aladin.md( "airmass@"+str(self.obs_time))
        
        for i,j in zip(airmass_start, airmass_end):
           
            #print i,j, color
            snap_1 = contour_ipix[(contour_ipix['airmass'] >= i) & (contour_ipix['airmass'] < j) ]
            moc_order = self.moc.moc_order(nside)
            moc = MOC.from_table( snap_1, 'RA[deg]', 'DEC[deg]',
                                  moc_order )                # moc creation
            moc.write( 'snap_airmass_'+'i', format = 'fits' )     # fits file

            if len(snap_1)!=0:
                aladin.send_file('snap_airmass_'+'i')
                #aladin.rename("cc")
                #airmass_plane = str(i) + "< airmass <" + str(j)
                aladin.rename(str(i)+"=<airmass<"+str(j))
                #aladin.set_planeID(str(i)+"=<airmass<"+str(j))
                aladin.mv(str(i)+'=<airmass<'+str(j), '"'+'airmass@'+str(self.obs_time)+'"' )
                

    def snapshot_airmass_up(self):
        """Showing a snapshot af the airmass at a given time using 5 colors."""

        time_update = self.update_count()
        #ipixs_update, moc_update = self.from_ipixs_to_moc(time_update)

        
        #obs_time = Time(self.user.get_obs_time())
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(time_update)
        print (time_update)
############################################################################################        
        snap_1 = contour_ipix[(contour_ipix['airmass'] >= 1) & (contour_ipix['airmass'] < 2) ]
        #print snap_1     
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_1, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_1', format = 'fits' )     # fits file

        if len(snap_1)!=0:
            aladin.send_file('snap_airmass_1')
            aladin.draw_airmass(140, 54, 1-2)
#####################################################################################################
        snap_2 = contour_ipix[(contour_ipix['airmass'] >= 2) & (contour_ipix['airmass'] < 3) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_2, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_2', format = 'fits' )     # fits file

        if len(snap_2)!=0:
            aladin.send_file('snap_airmass_2')
############################################################################################
        snap_3 = contour_ipix[(contour_ipix['airmass'] >= 3) & (contour_ipix['airmass'] < 4) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_3, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_3', format = 'fits' )     # fits file

        if len(snap_3)!=0:
            aladin.send_file('snap_airmass_3')
########################################################################
        snap_4 = contour_ipix[(contour_ipix['airmass'] >= 4) & (contour_ipix['airmass'] < 5.8) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_4, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_4', format = 'fits' )     # fits file

        if len(snap_4)!=0:
            aladin.send_file('snap_airmass_4')

    def close_window(self):
        self.destroy()



class ObsInMOC(Toplevel):
    """Initializi"""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()
        self.moc = MOC_confidence_region()

        self.observatory = astropy.coordinates.EarthLocation(
           lat=self.user.get_latitude()*u.deg, lon=self.user.get_longitude()*u.deg,
           height=self.user.get_altitude()*u.m)

        self.wait_visibility()
        
        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title("Observability" + " " + "starting from" + " " + self.user.get_obs_time() + " " + "UT")
        self.attributes("-topmost", True)

        self.bttn_clicks = 0 # counter ">>" Button
        
        # first label
        self.label_1 = Label(self, text="Show the region in the",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, padx=0)

        moc_value = 90  # default     
        moc_default = StringVar(self, value=moc_value)
        
        self.entry_percentage = Entry(self, width=5, justify=CENTER,
                             textvariable=moc_default)
        self.entry_percentage.grid(row=0, padx=2, column=1)

        # second label
        self.label_2 = Label(self, text="% MOC in which the airmass is ≤",
                             bg="slate grey")
        self.label_2.grid(row=0, column=2, sticky=E, pady=0)

        airmass_value = "2.5" # default
        airmass_default = StringVar(self, value=airmass_value)
        
        self.entry_airmass = Entry(self, width=5, justify=CENTER,
                             textvariable=airmass_default)
        self.entry_airmass.grid(row=0, padx=2, column=3)

        #Btn
        self.show = Button(self, text='Show',
                           command=self.snapshot_airmass)
        self.show.grid(column=4, row=0, sticky=W, padx=2, pady=5)
        
        self.moon = Button(self, text="Sun/Moon",
                            command=self.get_moon_position)  
        self.moon.grid(column=6,row=0, sticky=W, padx=2, pady=5) 
        
        self.forward = Button(self, text=">>",      
                               command=self.snapshot_airmass_up)
        self.forward.grid(column=5,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=7,row=0, sticky=E, padx=2, pady=5)


    #Actions
    def update_count(self):
        """Return the time in step of 1h when the button ">>" is clicked."""
   
        self.bttn_clicks += 1
              
        dt = TimeDelta(3600.0, format='sec')
        update_time = int(self.bttn_clicks) * dt
        
        obs_time = Time(self.user.get_obs_time())
        time_update = obs_time + update_time
        
        return time_update

    def from_ipixs_to_moc(self, time_input):
        """Return ipix table with the associated airmass"""

        prob = self.moc.read_prob(self.user.get_skymap())
        percentage = float(self.entry_percentage.get())/100.0
        
        ipixs = self.moc.ipixs_in_percentage(prob, percentage )
        nside = int(self.user.get_nside())
        
        ra, dec = self.moc.sky_coords(ipixs, nside)
        
        sky_coord = SkyCoord(ra = ra*u.deg,dec=dec*u.deg, frame='icrs')
        altaz = sky_coord.transform_to(AltAz(obstime=time_input, location=self.observatory))
        airmass_values = altaz.secz
        
        contour_ipix = Table([ ra, dec, airmass_values, ipixs ],
                             names = ('RA[deg]', 'DEC[deg]', 'airmass', 'ipix'),
                             meta = {'ipix': 'ipix table'})             # astropy table       

        return contour_ipix


    def snapshot_airmass(self):
        """Showing a snapshot af the airmass at a given time using 5 colors."""

        
        obs_time = Time(self.user.get_obs_time())
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(obs_time)
############################################################################################        
        snap_1 = contour_ipix[(contour_ipix['airmass'] >= 1) & (contour_ipix['airmass'] < 2) ]
        #print snap_1     
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_1, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_1', format = 'fits' )     # fits file

        if len(snap_1)!=0:
            aladin.send_file('snap_airmass_1')
#####################################################################################################
        snap_2 = contour_ipix[(contour_ipix['airmass'] >= 2) & (contour_ipix['airmass'] < 3) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_2, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_2', format = 'fits' )     # fits file

        if len(snap_2)!=0:
            aladin.send_file('snap_airmass_2')
############################################################################################
        snap_3 = contour_ipix[(contour_ipix['airmass'] >= 3) & (contour_ipix['airmass'] < 4) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_3, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_3', format = 'fits' )     # fits file

        if len(snap_3)!=0:
            aladin.send_file('snap_airmass_3')
########################################################################
        snap_4 = contour_ipix[(contour_ipix['airmass'] >= 4) & (contour_ipix['airmass'] < 5.8) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_4, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_4', format = 'fits' )     # fits file

        if len(snap_4)!=0:
            aladin.send_file('snap_airmass_4')


    def snapshot_airmass_up(self):
        """Showing a snapshot af the airmass at a given time using 5 colors."""

        time_update = self.update_count()
        #ipixs_update, moc_update = self.from_ipixs_to_moc(time_update)

        
        #obs_time = Time(self.user.get_obs_time())
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(time_update)
        print (time_update)
############################################################################################        
        snap_1 = contour_ipix[(contour_ipix['airmass'] >= 1) & (contour_ipix['airmass'] < 2) ]
        #print snap_1     
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_1, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_1', format = 'fits' )     # fits file

        if len(snap_1)!=0:
            aladin.send_file('snap_airmass_1')
            aladin.draw_airmass(140, 54, 1-2)
#####################################################################################################
        snap_2 = contour_ipix[(contour_ipix['airmass'] >= 2) & (contour_ipix['airmass'] < 3) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_2, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_2', format = 'fits' )     # fits file

        if len(snap_2)!=0:
            aladin.send_file('snap_airmass_2')
############################################################################################
        snap_3 = contour_ipix[(contour_ipix['airmass'] >= 3) & (contour_ipix['airmass'] < 4) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_3, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_3', format = 'fits' )     # fits file

        if len(snap_3)!=0:
            aladin.send_file('snap_airmass_3')
########################################################################
        snap_4 = contour_ipix[(contour_ipix['airmass'] >= 4) & (contour_ipix['airmass'] < 5.8) ]
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_table( snap_4, 'RA[deg]', 'DEC[deg]',
                                      moc_order )                # moc creation
        moc.write( 'snap_airmass_4', format = 'fits' )     # fits file

        if len(snap_4)!=0:
            aladin.send_file('snap_airmass_4')
                      
##################################################################################################
##### codice funzionante
##
##        mask = (contour_ipix['airmass']) >= 1 # clearing
##        obs1 = contour_ipix[mask]
##        
##        mask2 = (obs1['airmass']) <= float(self.entry_airmass.get())  # airmass user values
##        obs = obs1[mask2]
##
##        # TEST
##        #print obs, "sono qui"
##        
##        nside = self.user.get_nside()
##
##        #if len(obs)!=0:
##        moc_order = self.moc.moc_order(nside)
##        moc = MOC.from_table( obs, 'RA[deg]', 'DEC[deg]',
##                               moc_order )                # moc creation
##            
##        moc.write( 'obs_airmass_', format = 'fits' )     # fits file
##
##        if len(obs)!=0:
##           aladin.send_file('obs_airmass_')
##         
##        return obs, moc
#####################################################################################################
      
    def moc_obs(self):
        """Return the MOC region in which the airmass is <= the airmass value defined by the user."""

        time_start = Time(self.user.get_obs_time())
        ipixs_start, moc_start = self.from_ipixs_to_moc(time_start)

        percentage = float(self.entry_percentage.get())/100.0

        if len(ipixs_start) !=0:
            aladin.rename('obs_airmass_'+self.entry_airmass.get()\
                      +'MOC_'+str(percentage)+ '@' + str(time_start.isot))

            # printing area
            square_degrees_sphere = (360.0**2)/pi           
            area_sq2 = round( ( moc_start.sky_fraction * square_degrees_sphere ), 1 )
            
            tkMessageBox.showinfo('MOC visibility'+ '@' + str(time_start.isot),
                                 '   sky coverage: ' + str(area_sq2)+ ' ' + 'sq. deg')          
        else:
            tkMessageBox.showinfo('MOC visibility'+ '@' + str(time_start.isot),
                                  'No region for the selected  airmass')

        #print time_start
        #print type (time_start)

        # TEST
        #print (time_start)

    def moc_obs_update(self):
        """Return the MOC region in which the airmass is <= the airmass value defined by the user."""

        time_update = self.update_count()
        ipixs_update, moc_update = self.from_ipixs_to_moc(time_update)

        percentage = float(self.entry_percentage.get())/100.0

        if len(ipixs_update) != 0:
            aladin.rename('obs_airmass_'+self.entry_airmass.get()+\
                          'MOC_'+str(percentage)+ '@' + str(time_update.isot))

            # printing area
            square_degrees_sphere = (360.0**2)/pi
            area_sq2 = round( ( moc_update.sky_fraction * square_degrees_sphere ), 1 )
            
            tkMessageBox.showinfo('MOC visibility'+ '@' + str(time_update.isot),
                                 '   sky coverage: ' + str(area_sq2)+ ' ' + 'sq. deg')
            
        else:
            tkMessageBox.showinfo('MOC visibility'+ '@' + str(time_update.isot),
                                  'No region for the selected  airmass')

        # TEST
        #print (time_update)

    def get_moon_position(self):
        moon = Moon()
        moon.sky_position()

        #moon.steps()

    def close_window(self):
        self.destroy()

        
class ShiftFoV(Toplevel):
    """Shifting the FoV footprint(s) from a consecutive cardinal direction (↕, ↕, ↔, ↔);
       The widget contains 2 entries and 2 Buttons."""
    
    def __init__(self):
        Toplevel.__init__(self, border=7, bg="slate grey")
        self.attributes("-topmost", True)
        self.wait_visibility()

        self.user = UserValues() # get trasparency windows
        trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)     
        
        self.label_3 = Label(self, text="↕ (°)",bg="slate grey")  
        self.entry_3 = Entry(self, width=6, justify=CENTER)
        self.label_3.grid(row=0, sticky=E) 
        self.entry_3.grid(row=0, column=1)

        self.label_4 = Label(self, text="↔ (°)",bg="slate grey")  
        self.entry_4 = Entry(self, width=6, justify=CENTER)
        self.label_4.grid(row=0,column=3) 
        self.entry_4.grid(row=0, column=4)

        self.close = Button(self, text="Close",
                            command = self.close_window)  
        self.close.grid(column=4,row=2)
    
    def north_shift(self):
        self.title(" Shifting - North")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_north)
        self.checkbox.grid(column=3,row=2)

    def south_shift(self):
        self.title(" Shifting - South")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_south)
        self.checkbox.grid(column=3,row=2)

    def east_shift(self):
        self.title(" Shifting - East")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_east)
        self.checkbox.grid(column=3,row=2)

    def west_shift(self):
        self.title(" Shifting - West")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_west)
        self.checkbox.grid(column=3,row=2)
        
    # Actions
    def shift_north(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.north(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_south(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.south(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_east(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.east(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_west(self):              
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.west(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def close_window(self):
        self.destroy()

class FoVstatistics(Toplevel, SkyCoverage):
    """ FoV statistics window consists of 3 buttons:
         (1) Confirm the FoV-footprint
         (2) Delete FoV-footprint
         (3) Zoom in FoV. """
    
    def __init__(self):
        Toplevel.__init__(self)
        SkyCoverage.__init__(self)

        self.title("FoV statistics")
        self.attributes("-topmost", True)
              
        self.B00 = Button(self,text="Confirm the pointing in the GWsky_pointings txt file")  
        self.B00.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B00.pack(side=TOP,fill=BOTH)

        self.B01 = Button(self,text="Delete the FoV")  
        self.B01.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B01.pack(side=TOP,fill=BOTH)

        self.B02 = Button(self,text="Zoom in the FoV")  
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.pack(side=TOP,fill=BOTH)

        # default
        fov_center_ra_dec = str(self.input_ra), str(self.input_dec), 'on fly notes, NO DELETE OR CHANGE THE COORDS!'       
        current_fov = StringVar(self, value=fov_center_ra_dec)  
        self.entry_current_fov = Entry(self, width=30, justify=LEFT,
                                       textvariable=current_fov)
        
        #self.entry_current_fov.bind("<Key>", lambda e: "break") # saving entries
        self.entry_current_fov.pack(side=TOP,fill=BOTH)

    def __rm_from_stack(self, ra, dec):
        """Removing from Aladin stack the associated planes."""
        
        aladin.remove_FoV(ra, dec)       
        aladin.remove("Q:"+ ra +"/"+ dec)           
        aladin.remove("C_" + ra+ "/" + dec)

    def clicked(self, event):
        """Retain or delate the FoV-footprint(s).
            The FoV-center positions are saved in "GWsky_pointings.txt" """      
        
        if event.widget == self.B00: # Retain and Close the FoV.
                                     # creating FoV moc
            current_fov_coords = self.entry_current_fov.get().split() # getting entries
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]
            #print current_fov_ra, current_fov_dec
            #print self.user.get_fov_shape()
            #prob_fov = self.lvc.prob_in_box(
            #   ra, dec, self.user.get_fov_width(), self.user.get_fov_height()) 
            #if self.user.get_fov_shape() != 2:  # box FoV
            #print self.user.get_fov_width(), self.user.get_fov_height()

            self.moc = MOC_confidence_region()
            
            ipixs, NSIDE = self.moc.ipix_in_box(float(current_fov_ra), float(current_fov_dec),
                                 float(self.user.get_fov_width()), float(self.user.get_fov_height()))
            #import healpy as hp
            #nside=512
            #theta, phi = hp.pix2ang(nside, list(ipixs))
            #print type(ipixs)
            ra,dec = self.moc.sky_coords(ipixs, NSIDE)
            #print coords
            # creating an astropy.table with RA[deg] and DEC[deg]
            box_ipix = Table([ra, dec], names = ('RA[deg]', 'DEC[deg]'),
                                 meta = {'ipix': 'ipix table'})

            moc_order = self.moc.moc_order(NSIDE)
            moc = MOC.from_table( box_ipix, 'RA[deg]', 'DEC[deg]',
                                  moc_order )                # moc creation
            
            moc.write( 'z_fov', format = 'fits' )     # fits file

            # progressive union.
            # cobn gli entry  calcoli i vertici
            #fai una query
            # la trasformi in ra dec
            # moc table
            
            self.destroy()
            
        if event.widget == self.B01: # Delete FoV
            current_fov_coords = self.entry_current_fov.get().split() # getting entries
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]
            
            self.__rm_from_stack(current_fov_ra, current_fov_dec)
            
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                 ra=str(current_fov_ra), dec=str(current_fov_dec))           
            self.destroy()
       
        if event.widget == self.B02:  # Zoom in the  FoV
            aladin.location(str(self.input_ra), str(self.input_dec))             
            aladin.zoom('1x')

    def __are_all_same(self, items):
        """Check if all elements of a list are the same."""
        
        return all(x == items[0] for x in items)
      
    # Plots        
    def plot_stats(self, ra,dec, moon_illumination, sep_fov_moon, table_q,
                   r, dp_dr, airmass_values, datestrings, prob_fov):
        """Showing the plots in the FoV statistic window."""
        
        f = Figure(figsize=(9, 5.2), facecolor='white')
        f.subplots_adjust(left=.13, bottom=.16, right=.93, top=.84, wspace=.26, hspace=.3)
        

        def airmass_subplot():
            """SubPlot Airmass."""
            
            airmass_plt = f.add_subplot(223)
            airmass_plt.set_ylabel('airmass')
            
            same = self.__are_all_same(airmass_values)

            if same !=True:
                time_step = [dateutil.parser.parse(s) for s in datestrings]
            
                ax=plt.gca()
                ax.xaxis.set_major_formatter(
                mdates.DateFormatter('%H:%M:%S'))

                airmass_plt.set_ylabel('airmass')
                airmass_plt.set_xlabel('Universal Time')
                airmass_plt.invert_yaxis()
                airmass_plt.grid(True)
           
                airmass_plt.plot(time_step, airmass_values, "o-", linestyle="--")
                f.autofmt_xdate()
            else:
                tkMessageBox.showinfo('Warning',"airmass outside the range of 1 - 5.8")
                pass

        airmass_subplot = airmass_subplot()

        def subplot_cat_distrib():
            """SubPlot catalog distributions."""

            # setting column 1
            query_catalog_1 = f.add_subplot(222)
            query_catalog_1.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                       'col: ' + self.user.get_column_1())
            query_catalog_1.set_ylabel('Count')    

            # setting column 2
            query_catalog_2 = f.add_subplot(224)
            query_catalog_2.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                       'col: ' + self.user.get_column_2())
            query_catalog_2.set_ylabel('Count')         
            
            for table_name in table_q.keys():
                table = table_q[table_name]
           
                try:
                    query_catalog_1.hist(table[self.user.get_column_1()])
                    query_catalog_2.hist(table[self.user.get_column_2()])
         
                except KeyError as key_error:
                    tkMessageBox.showerror(' Error: no key found', key_error)
                except ValueError:
                    c1, c2  = np.array(table[self.user.get_column_1()]), \
                              np.array(table[self.user.get_column_2()])
                    newc1, newc2 = c1[~np.isnan(c1)], \
                                   c2[~np.isnan(c2)]
                    query_catalog_1.hist(newc1)
                    query_catalog_2.hist(newc2)

                table.write("GWsky_query_items", format =
                            'votable', overwrite = True)
                #aladin.send_file("GWsky_query_items")
                #aladin.rename("Q:"+str(ra)+"/"+str(dec))
            
        subplot_cat_distrib = subplot_cat_distrib()       

        def subplot_cond_dist():
            """Conditional Distance Distribution Along a Line of Sight (FoV center position)."""
            
            conditional_distance_line_sight = f.add_subplot(221) 
            conditional_distance_line_sight.plot(r, dp_dr)
            
            title_string = 'Conditional distance distribution \n along the FoV center'
            conditional_distance_line_sight.set_title(title_string,fontsize=10)
            conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        subplot_cond_dist = subplot_cond_dist()

        def draw_area():
            """Drawing Area of the FoV Statistic window."""

            fov_information_title = "FoV center (ra "+str(ra) + "  " +"dec "+ str(dec)+")" + "; " + "prob: " + str(prob_fov)+ \
                                    ";" + " " + "Moon" + " " + "(illumi.:" + " " + str(moon_illumination) + " " + \
                                    "dist.:" + " " + str(sep_fov_moon) + ")"
            
            f.suptitle(fov_information_title, fontsize=10)
           
            canvas = FigureCanvasTkAgg(f, self) 
            canvas.show()
            canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

            toolbar = NavigationToolbar2TkAgg(canvas, self)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
            
        draw_area = draw_area()

        self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov, skymap=self.user.get_skymap())

        aladin.send_file("GWsky_query_items")
        aladin.rename("Q:"+str(ra)+"/"+str(dec))
        
        aladin.draw_newtool("C_" + str( ra ) + "/" + str( dec ))
        aladin.draw_string(ra, dec, str( ra ) + "/" + str( dec ))
        print (str(sep_fov_moon))
        print (sep_fov_moon.round(1))

def on_closing():
    """Asking the closure of the coverage window. If "Quit" the files in the list "temp_files" are deleted.
          ***Improving with tempfile module***"""
    
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
        try:
            temp_files=["GWsky_entries", "GWsky_query_items", "GWsky_fov.vot",
                        "GWsky_config", "GWsky_coords", "obs_airmass_"]
            for temp_file in temp_files:
               os.remove(temp_file)
        except OSError:
            pass
        mainWin.destroy()


# running
# ***TO BE COMPLETED***

mainWin = Tk()
sscGUI = SkyCoverageGUI(mainWin)

mainWin.title('GWsky')
mainWin.attributes("-topmost", True)

mainWin.wait_visibility(mainWin)

user = UserValues() # getting transparency
trasparency = user.get_win_trasparency()
mainWin.wm_attributes('-alpha', trasparency) 

mainWin.protocol("WM_DELETE_WINDOW", on_closing)
mainWin.mainloop()
