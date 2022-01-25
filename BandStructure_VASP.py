from tkinter import OptionMenu, Button, Checkbutton, Radiobutton, Label, Entry, Text, Menu, Frame
from tkinter import Tk, Toplevel, colorchooser
from PIL import ImageTk, Image
from tkinter import INSERT, END, RIDGE, NORMAL, DISABLED, SUNKEN
from tkinter import messagebox, filedialog
from tkinter import StringVar, IntVar, DoubleVar, BooleanVar
from tkinter import font as tkFont

import json
import numpy as np
import os
import datetime
import re 
from Kpoints_new import K_points

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib.collections import LineCollection
from matplotlib.collections import PathCollection
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


class FullScreenApp(object):
    """
    Class to increase the window to the maximum screen size
    """

    def __init__(self, screen, **kwargs):
        """
        Change screen to fullscreen
        Input:
        ----------------
        screen: The window to maximize
        """
        self.screen = screen
        edge = 3
        self._small = '400x200+0+0'
        screen.geometry("{0}x{1}+0+0".format(
            screen.winfo_screenwidth() - edge, screen.winfo_screenheight() - edge))
        screen.bind('<Escape>', self.toggle_screen)


    def toggle_screen(self, event):
        small = self.screen.winfo_geometry()
        self.screen.geometry(self._small)
        self._small = small


class Help:
    """
    Help Menu to link to other websites/documentations
    """

    def __init__(self):
        """
        Help menu window
        """
        self.screen_help = Toplevel()
        self.screen_help.configure(bg = MainApplication._from_rgb(self, (11, 165, 193)))
        self.screen_help.geometry('700x300')
        self.screen_help.iconbitmap('icon_band.ico')

    def welcome(self):
        """
        Welcome window
        """
        text = Text(self.screen_help, height = 15)
        text.insert(INSERT, 'Welcome to the VASP Electronic Band Structure Plotter! \n')
        text.insert(INSERT, '\n')
        text.insert(INSERT, 'You can create KPOINT files and plot electronic band structures, density of     states (DOS), and partial DOS which are computed using VASP. ')
        text.insert(INSERT, ' \n')
        text.insert(INSERT, '\n')
        text.insert(END, ' Please check out the documentaries for more information!  \n \n Thank you for choosing the VASP Electronic Structure Plotter')
        text.grid(row = 0, column = 0, padx = 10, pady = (30, 10))

    def documentary(self, *args):
        """
        Documentary window
        """
        text = Text(self.screen_help, height = 15)
        text.insert(INSERT, 'Welcome to the VASP Electronic Band Structure Plotter! \n')
        text.insert(INSERT, '\n')
        text.insert(END, 'I will update the documentaries at a later point!')
        text.grid(row = 0, column = 0, padx = 10, pady = (30, 10))

    def about(self):
        """
        About window
        """
        text = Text(self.screen_help)
        text.insert(INSERT, 'This is a program to plot electronic band structures, density of states (DOS),  and partial DOS, computed using the commercial software VASP. \n')
        text.insert(INSERT, '\n')
        text.insert(INSERT, 'The software was written in Python and Tkinter!  This is version v1.0 and I am  looking for any suggestions and reports of errors. \n')
        text.insert(INSERT, '\n')
        text.insert(INSERT, 'This is a free software and should not be used for commercial reasons. \n')
        text.insert(INSERT, '\n')
        text.insert(INSERT, 'If you have suggestions, concerns, or find errors, please send me an email: \n Jan.Poehls@dal.ca \n ')
        text.insert(INSERT, '\n')
        text.insert(INSERT, 'I would like to acknowledge the FRQNT PBEEE postdoctoral fellowship! \n \n')
        text.insert(INSERT, 'Thank you for choosing the VASP Electronic Band Structure Plotter. \n \n  --Jan-- \n \n')
        text.insert(INSERT,  '\xa9 Jan-Hendrik Poehls, PhD, MSc, BSc, 2021')
        text.grid(row = 0, column = 0, padx = 10, pady = (30, 10))


class Energy:
    """
    Get information and parameters from PROCAR files
    """

    def __init__(self, procar, kpoints=list(), energy=list(), occupation=list(), total_DOS=list(), coordinates=list(),
        weight=list(), DOS_elements=list(), DOS_orbitals=list(), DOS_element_new=list()):
        """
        Get energy, DOS, and other parameters from PROCAR file
        Input:
        -----------------------
        procar: name of procar file, str
        kpts: ndarray, shape (N), dtype=int
            array of kpoints
        energy: ndarray, shape (M, N), dtype=float
            array of energy where M is number of bands and N is number of kpoints
        occ: ndarray, shape (M, N), dtype=float
            array of occupations
        totDOS: ndarray, shape (M, N), dtype=float
            array of totalDOS
        coord: ndarray, shape (N, 3), dtype=float
            array of coordinates
        weight: ndarray, shape (N)
            array of weight of each kpoint (degeneracy)
        DOS_elements: ndarray, shape (M, N, Ions), dtype=float
            array of elemental DOS where Ions is the number of ions
        DOS_orbitals: ndarray, shape (M, N, orb), dtype=float
            array of orbital DOS where orb is the number of different orbitals (s, p, d, f)
        DOS_element_new: ndarray, shape (M, N, Cmp), dtype=float
            array of elemental DOS summing up the same element where Cmp is the number of elements
        """
        self.procar = procar
        self.kpts = kpoints
        self.energy = energy
        self.occ = occupation
        self.totDOS = total_DOS
        self.coord = coordinates
        self.weight = weight
        self.DOS_elements = DOS_elements
        self.DOS_orbitals = DOS_orbitals
        self.DOS_element_new = DOS_element_new


    def get_energies(self):
        """
        Get information from PROCAR file
        """
        Input_pro = self.procar[1].split()

        if Input_pro[3] != "#":
            Nmb_kpts = int(Input_pro[3])
            Nmb_bands = int(Input_pro[7])
            Nmb_ions = int(Input_pro[11])

        else:
            Nmb_kpts = int(re.split('(\d+)', Input_pro[2])[1])
            Nmb_bands = int(Input_pro[6])
            Nmb_ions = int(Input_pro[10])

        self.kpts = np.zeros(Nmb_kpts, dtype=int); self.coord = np.zeros((Nmb_kpts, 3), dtype=float)
        self.weight = np.zeros(Nmb_kpts, dtype=float); self.energy = np.zeros((Nmb_bands, Nmb_kpts), dtype=float)
        self.occ = np.zeros((Nmb_bands, Nmb_kpts), dtype=float); self.totDOS = np.zeros((Nmb_bands, Nmb_kpts), dtype=float)
        self.DOS_elements = np.zeros((Nmb_bands, Nmb_kpts, Nmb_ions), dtype=float)
        self.DOS_orbitals = np.zeros((Nmb_bands, Nmb_kpts, len(self.procar[7].split()) - 2), dtype=float)

        for i in range(int(Nmb_kpts)):

            line = ((int(Nmb_ions) + 5) * int(Nmb_bands) + 5) * i + 3 - i * 2
            self.kpts[i] = self.procar[line].split()[1]

            for xx in range(3):
                self.coord[i][xx] = self.procar[line].split()[3 + xx]

            self.weight[i] = self.procar[line].split()[8]


            for j in range(int(Nmb_bands)):
                lines = (int(Nmb_ions) + 5) * j + 2 + line
                self.energy[j][i] = self.procar[lines].split()[4]
                self.occ[j][i] = self.procar[lines].split()[7]
                self.totDOS[j][i] = self.procar[lines + 3 + int(Nmb_ions)].split()[len(self.procar[7].split()) - 1]

                for k in range(int(Nmb_ions)):
                    self.DOS_elements[j][i][k] = self.procar[lines + 3 + k].split()[-1]
                for l in range(len(self.procar[7].split()) - 2):
                    self.DOS_orbitals[j][i][l] = self.procar[lines + 3 + int(Nmb_ions)].split()[l + 1]


    def get_band_gap(self):
        """
        Get the band gap in eV
        """
        VBM = -1000; CBM = 1000
        for i in range(len(self.energy)):

            for k in range(len(self.energy[0])):

                if self.occ[i][k] > 0.01:
                    if VBM < self.energy[i][k]:
                        VBM = self.energy[i][k]

                else:
                    if CBM > self.energy[i][k]:
                        CBM = self.energy[i][k]

        return CBM - VBM, VBM


    def tick_label(self, kpoint):
        """
        Write elements in LaTeX format
        """
        if len(kpoint.split()[-1]) > 1:
            return "${}$".format(kpoint.split()[-1])
        else:
            return '{}'.format(kpoint.split()[-1])


    def get_distance(self, Kpoint_mesh, kpoints):
        """
        Write the list of ticks and position of the ticks in the plot
        Input:
        --------------------------
        Kpoint_mesh: list
            List of number of kpoints between high-symmetry points
        kpoints: Lines from KPOINTS file

        Output:
        -------------------------
        ticks: List
            List of labels of high-symmetry points
        distance: List
            List of distances between 0 and 1
        """

        ticks = []; distance_new = [0]
        for i in range(len(Kpoint_mesh)):

            distance_new.append(Kpoint_mesh[i] + distance_new[-1])

        ticks.append(self.tick_label(kpoints[3]))

        for k in range(4, len(kpoints)):

            if len(kpoints[k].split()) < 1 and k + 1 != len(kpoints):

                if kpoints[k - 1].split()[-1] == kpoints[k + 1].split()[-1]:
                    ticks.append(self.tick_label(kpoints[k - 1]))

                else:
                    ticks.append('{}$\\mid$ {}'.format(self.tick_label(kpoints[k -1]), self.tick_label(kpoints[k + 1])))

            elif k + 1 == len(kpoints):
                ticks.append(self.tick_label(kpoints[k - 1]))

        distance = [x / distance_new[-1] for x in distance_new]
        self.kpts = [x / distance_new[-1] for x in self.kpts]

        return ticks, distance


    def element_DOS(self, contcar):
        """
        Sum up all elemental DOS of the same element
        Input:
        ----------------------------
        contcar: CONTCAR file includes a list and order of the ions
        """
        Nmb_Cmp = np.array(contcar[6].split(), dtype=int)
        self.DOS_element_new = np.zeros((len(self.DOS_elements), len(self.DOS_elements[0]), len(Nmb_Cmp)), dtype=float)

        for i in range(len(self.DOS_elements)):
            for j in range(len(self.DOS_elements[i])):
                counter = 0
                for c in range(len(Nmb_Cmp)):
                    for k in range(Nmb_Cmp[c]):
                        self.DOS_element_new[i][j][c] += self.DOS_elements[i][j][counter]
                        counter += 1


    def sum_partial_DOS(self, totalDOS, minE, maxE, Eres):
        """
        Get partial DOS over the entire Brillouin zone
        Input:
        --------------------------
        totalDOS: ndarray, shape (M, N), dtype=float
            total DOS (elemental, orbital) for M bands and N kpoints
        minE: float
            minimum energy in eV
        maxE: float
            maximum energy in eV
        Eres: float
            stepsize to solve
        """
        steps = int((maxE - minE) / Eres)
        Energy_DOS = np.zeros(steps, dtype=float); TOTAL_DOS = np.zeros_like(Energy_DOS)

        for i in range(steps):
            Energy_DOS[i] = minE + i * Eres + 0.001

            for j in range(len(self.energy)):
                for k in range(len(self.energy[0])):
                    if  Energy_DOS[i] -0.5 * Eres < self.energy[j][k] < Energy_DOS[i] + 0.5 * Eres:
                        TOTAL_DOS[i] += totalDOS[j][k] * self.weight[k]

        return Energy_DOS, TOTAL_DOS


class EntryItem:
    """
    Combine different tkinter items with Label items
    """

    def __init__(self, parent, name, row = 0, column = 1, padx = 10, pady = 6, width = 20, columnspan = 1, state = NORMAL, ipadx = 0, options = ['0']):
        """
        Combine Entry/Menu items with Label items
        Input:
        -------------------------------------
        parent: class
        name: str
            Name of the Label
        row: int
            number of row
        column: int
            number of column
        padx: int
            space in x-direction of grid
        pady: int
            space in y-direction of grid
        width: int
            width of the Entry
        columnspan: int
            column span of Entry and Label
        state: NORMAL or DISABLED
        options: list
            List of options for OptionMenu
        initial_val: StringVar()
            initial value of OptionMenu
        var: StringVar()
            name and type of variable of Entry
        entry: Entry item
        label: Label item
        menu: OptionMenu item
        """
        self.parent = parent
        self.name = name
        self.row = row
        self.column = column
        self.padx = padx
        self.pady = pady
        self.width = width
        self.columnspan = columnspan
        self.state = state
        self.ipadx = ipadx
        self.options = options
        self.initial_val = StringVar()
        self.var = StringVar()
        self.entry = Entry(self.parent, textvariable = self.var, state = self.state, width = self.width)
        self.label = Label(self.parent, text = name, relief = RIDGE, anchor = 'w')
        self.menu = OptionMenu(self.parent, self.initial_val, *self.options)


    def create_EntryItem(self, padx_label = 10, pady_label = 6, ipadx_label = 0):
        """
        Combine Label with Entry item next to each other
        Input:
        ------------------------------
        padx_label: int
            space in x-direction for the Label
        pady_label: int
            space in y-direction for the Label
        ipadx_label: int
            size of the Label
        """
        self.entry.grid(row = self.row, column = self. column, columnspan = self.columnspan, padx = self.padx, pady = self.pady, ipadx = self.ipadx)
        self.label.grid(row = self.row, column = self.column - 1, columnspan = self.columnspan, padx = padx_label, pady = pady_label, ipadx = ipadx_label)


    def set_name(self, new_name = 'NaN'):
        """
        Change the name of the Entry
        Input:
        ------------------
        new_name: str
            Name of the Label
        """
        self.var.set(new_name)


    def get_name(self):
        """
        Get value of the Entry
        """
        return self.var.get()


    def delete(self):
        """
        Delete an Entry
        """
        self.entry.delete(0, END)


    def update(self, name):
        """
        Remove entry and insert a new entry with a new value
        Input:
        ----------------------
        name: str
            New value of Entry
        """
        self.delete()
        self.entry.insert(0, name)
        return


class MainApplication:
    """
    Main window and functions to plot VASP data
    """

    def __init__(self, parent, *args, **kwargs):
        self.parent = parent
        self.parent.configure(bg = self._from_rgb((11, 165, 193)))
        self.title = self.parent.title('VASP Electronic Band Structure App')
        self.icon = self.parent.iconbitmap('icon_band.ico')
        self.font_window = tkFont.Font(family='Helvetica', size= 10, weight='bold')

        # Create Frame
        self.input = Frame(self.parent, height = 308, width = 395, bg = self._from_rgb((11, 112, 141)))
        self.input.grid(row = 0, column = 0, columnspan = 2, rowspan = 9, pady = (10, 5))
        self.label_input = Label(self.parent, text = 'Input parameters')
        self.label_input.grid(row = 0, column = 0, pady = (10, 5))
        self.label_input['font'] = self.font_window
        self.output = Frame(self.parent, height = 302, width = 395, bg = self._from_rgb((175, 188, 205)))
        self.output.grid(row = 9, column = 0, columnspan = 2, rowspan = 9, pady = (0, 5))
        self.label_pDOS = Label(self.parent, text = 'Partial Density of States')
        self.label_pDOS.grid(row = 9, column = 0, pady = (0, 5))
        self.label_pDOS['font'] = self.font_window

        # Create MenuBar
        my_Menu = Menu(self.parent)
        self.parent.config(menu = my_Menu)

        file_menu = Menu(my_Menu)
        my_Menu.add_cascade(label = 'File', menu = file_menu)
        file_menu.add_command(label = 'New', command = self.clear)
        file_menu.add_command(label = 'Open File', command = self.open_file)
        file_menu.add_separator()
        file_menu.add_command(label = 'Exit', command = self.close_program)

        edit_menu = Menu(my_Menu)
        my_Menu.add_cascade(label = 'Edit', menu = edit_menu)
        edit_menu.add_command(label = 'Edit Graph', command = self.Edit_graph)

        kpoint_menu = Menu(my_Menu)
        my_Menu.add_cascade(label = 'KPOINTS', menu = kpoint_menu)
        kpoint_menu.add_command(label='Create KPOINTS', command=self.create_kpoint)

        plot_menu = Menu(my_Menu)
        my_Menu.add_cascade(label = 'Plot', menu = plot_menu)
        plot_menu.add_command(label='Plot VASP', command=self.plot_electronic_structure)

        help_menu = Menu(my_Menu)
        my_Menu.add_cascade(label = 'Help', menu = help_menu)
        help_menu.add_command(label = 'Welcome', command = self.welcome)
        help_menu.add_command(label = 'Documentations', command = self.documentary)
        help_menu.add_separator()
        help_menu.add_command(label = 'About', command = self.about)

        self.app = FullScreenApp(self.parent)

        # Create Plot Data
        self.font_size_band_x = DoubleVar(); self.font_size_band_y = DoubleVar()
        self.font_size_band_ticks = DoubleVar(); self.font_size_band_energy = DoubleVar()
        self.font_size_DOS_x = DoubleVar(); self.font_size_DOS_y = DoubleVar()
        self.font_size_DOS_number = DoubleVar()
        self.size_x = DoubleVar(); self.size_y = DoubleVar()
        self.size_x_space = DoubleVar(); self.size_x_length = DoubleVar()
        self.size_y_space = DoubleVar(); self.size_y_length = DoubleVar()

        self.filename = EntryItem(self.parent, name = 'Filename', row = 1, state=DISABLED)
        self.filename.create_EntryItem(ipadx_label = 45)
        self.minE = EntryItem(self.parent, name = 'min. Energy / eV', row = 2)
        self.minE.create_EntryItem(ipadx_label = 26); self.minE.initial_val = DoubleVar()
        self.maxE = EntryItem(self.parent, name = 'max. Energy / eV', row = 3)
        self.maxE.create_EntryItem(ipadx_label = 26); self.maxE.initial_val = DoubleVar()
        self.Eres = EntryItem(self.parent, name = 'Energy resolution / eV', row = 4)
        self.Eres.create_EntryItem(ipadx_label = 12); self.Eres.initial_val = DoubleVar()
        self.ymax = EntryItem(self.parent, name = 'maximum DOS', row = 5)
        self.ymax.create_EntryItem(ipadx_label = 29); self.ymax.set_name('2')

        self.pDOS = IntVar(); self.pDOS.set(1); self.pDOS_E_var = BooleanVar(); self.pDOS_O_var = BooleanVar()
        self.label_energy_var = BooleanVar(); self.label_ticks_var = BooleanVar()
        self.label_DOS_var = BooleanVar(); self.label_energy_DOS_var = BooleanVar()
        self.grid_energy_var = BooleanVar(); self.grid_DOS_var = BooleanVar()
        self.ticks_energy_var = BooleanVar(); self.ticks_wavevector_var = BooleanVar()
        self.ticks_energy_DOS_var = BooleanVar(); self.ticks_DOS_var = BooleanVar()
        self.DOS_label = Label(self.parent, text='Choose DOS to plot')
        self.DOS_label.grid(row=10, column=0, ipadx=30)
        self.TotalDOS = Radiobutton(self.parent, text="Total DOS", variable=self.pDOS, value=1)
        self.TotalDOS.grid(row = 10, column = 1, ipadx=25)
        self.partialDOS = Radiobutton(self.parent, text="Elemental DOS", variable=self.pDOS, value=2)
        self.partialDOS.grid(row = 11, column = 1, ipadx=12)
        self.orbitalDOS = Radiobutton(self.parent, text="Orbital DOS", variable=self.pDOS, value=3)
        self.orbitalDOS.grid(row = 12, column = 1, ipadx=20)
        self.pDOS_E = Checkbutton(self.parent, text='Elemental DOS on Electronic Band Structure', variable=self.pDOS_E_var)
        self.pDOS_E.grid(row=13, column=0, columnspan=2, ipadx=20)
        self.pDOS_O = Checkbutton(self.parent, text='Orbital DOS on Electronic Band Structure', variable=self.pDOS_O_var)
        self.pDOS_O.grid(row=14, column=0, columnspan=2, ipadx=28)
        self.pDOS_E_var.trace('w', self.check_orbital); self.pDOS_O_var.trace('w', self.check_elemental)

        self.load_button = Button(self.parent, text='Load', command=self.load_electronic_properties, state=DISABLED)
        self.load_button.grid(row=15, column=0, pady=10, ipadx=20)
        self.load_button['font'] = self.font_window
        self.plot_button = Button(self.parent, text='Plot', command=self.plot_electronic_structure, state=DISABLED)
        self.plot_button.grid(row=15, column=1,  pady=10, ipadx=20)
        self.plot_button['font'] = self.font_window

        self.set_dpi = EntryItem(self.parent, name = 'dpi', row = 14, column=8)
        self.set_dpi.create_EntryItem(ipadx_label=40)
        self.save_figure_button = Button(self.parent, text='Save', command=self.save_electronic_structure, state=DISABLED)
        self.save_figure_button.grid(row=14, column=9, columnspan=2, pady=10, ipadx=28)
        self.save_figure_button['font'] = self.font_window

        self.save_fig_csv_button = Button(self.parent, text='Save as .csv', command=self.save_csv_file, state=DISABLED)
        self.save_fig_csv_button.grid(row=15, column=9, columnspan=2, pady=10, ipadx=4)
        self.save_fig_csv_button['font'] = self.font_window

        self.font_options = [
            'Times New Roman',
            'Arial',
            'Gabriola',
            'Courier New',
            'Cambria',
            'Calibri',
        ]
        self.initial_font = StringVar()

        self.color_2plot_options = [
            'red-blue',
            'red-green',
            'green-blue'
        ]
        self.initial_color_2plot = StringVar()
        self.foldername = ''

        self.initial_parameters()
        self.create_empty_plot()


    def check_orbital(self, *args):
        """
        If check elemental, turn off orbital
        """
        if self.pDOS_O_var.get() and self.pDOS_E_var.get():
            self.pDOS_O_var.set(False)


    def check_elemental(self, *args):
        """
        If check orbital, turn off elemental
        """
        if self.pDOS_E_var.get() and self.pDOS_O_var.get():
            self.pDOS_E_var.set(False)


    def initial_parameters(self):
        """
        Initial parameters for the plot, if '~default.json' exists parameters are taking from this file
        """

        if '~default.json' in os.listdir():
            with open('~default.json') as d:
                dic = json.load(d)

            self.font_size_band_x.set(dic['size_band_x']); self.font_size_band_y.set(dic['size_band_y'])
            self.font_size_band_ticks.set(dic['size_band_ticks']); self.font_size_band_energy.set(dic['size_band_energy'])
            self.font_size_DOS_x.set(dic['size_DOS_x']); self.font_size_DOS_y.set(dic['size_DOS_y'])
            self.font_size_DOS_number.set(dic['size_DOS_ticks'])
            self.size_x.set(dic['figure_size_x']); self.size_y.set(dic['figure_size_y'])
            self.size_x_space.set(dic['figure_space_x']); self.size_x_length.set(dic['figure_length_x'])
            self.size_y_space.set(dic['figure_space_y']); self.size_y_length.set(dic['figure_length_y'])
            self.label_energy_var.set(dic['label_energy']); self.label_energy_DOS_var.set(dic['label_energy_DOS'])
            self.label_ticks_var.set(dic['label_ticks']); self.label_DOS_var.set(dic['label_DOS'])
            self.grid_energy_var.set(dic['grid_energy']); self.grid_DOS_var.set(dic['grid_DOS'])
            self.ticks_energy_var.set(dic['ticks_energy']); self.ticks_wavevector_var.set(dic['ticks_wavevector'])
            self.ticks_energy_DOS_var.set(dic['ticks_energy_DOS']); self.ticks_DOS_var.set(dic['ticks_DOS'])
            self.color = dic['color']; self.hx = dic['color_hex']
            for i in range(len(self.font_options)):
                if dic['font'] == self.font_options[i]:
                    self.initial_font.set(self.font_options[i])
            for i in range(len(self.color_2plot_options)):
                if dic['color_2plot'] == self.color_2plot_options[i]:
                    self.initial_color_2plot.set(self.color_2plot_options[i])

        else:
            self.font_size_band_x.set(16); self.font_size_band_y.set(16)
            self.font_size_band_ticks.set(14); self.font_size_band_energy.set(16)
            self.font_size_DOS_x.set(16); self.font_size_DOS_y.set(16)
            self.font_size_DOS_number.set(16)
            self.size_x.set(8.5); self.size_y.set(5.)
            self.size_x_space.set(0.18); self.size_x_length.set(0.78)
            self.size_y_space.set(0.23); self.size_y_length.set(0.68)
            self.initial_font.set(self.font_options[0]); self.initial_color_2plot.set(self.color_2plot_options[0])
            self.label_energy_var.set(True); self.label_DOS_var.set(True)
            self.label_ticks_var.set(True); self.label_energy_DOS_var.set(False)
            self.grid_energy_var.set(True); self.grid_DOS_var.set(True)
            self.ticks_energy_var.set(True); self.ticks_wavevector_var.set(True)
            self.ticks_energy_DOS_var.set(True); self.ticks_DOS_var.set(True)
            self.color= np.array([0, 0, 0]); self.hx = '#000000'

        self.minE.set_name(-5); self.maxE.set_name(5)
        self.Eres.set_name(0.05); self.pDOS_E_var.set(False)
        self.pDOS_E.config(state=DISABLED); self.pDOS_O.config(state=DISABLED)
        self.set_dpi.set_name('100')


    def create_empty_plot(self):
        """
        Create an empty plot using default values
        """

        plt.rcParams["font.family"] = self.initial_font.get()
        plt.rcParams.update({'font.size': self.font_size_band_energy.get()})

        self.fig = Figure(figsize= (self.size_x.get(), self.size_y.get()), dpi = 100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.parent)
        self.canvas.draw()
        self.plot_widget = self.canvas.get_tk_widget()
        self.plot_widget.grid(row = 1, column = 3, columnspan=11, rowspan = 13)

        gs = self.fig.add_gridspec(1, 2, width_ratios=[2, 1,])
        gs.update(left=0.1, right=0.95, wspace=0.15)
        self.ax1 = self.fig.add_subplot(gs[0])
        self.ax2 = self.fig.add_subplot(gs[1])

        gs2 = self.fig.add_gridspec(2, 2, width_ratios=[1, 10,], height_ratios=[1, 8,])
        gs2.update(bottom=0.3, top =0.95)
        self.ax3 = self.fig.add_subplot(gs2[0])
        self.ax3.axis('off')
        if self.label_ticks_var.get():
            self.ax1.set_xlabel('Wavevector $k$', fontsize=self.font_size_band_x.get(), family=self.initial_font.get())
        if self.label_energy_var.get():
            self.ax1.set_ylabel('Energy / eV', fontsize=self.font_size_band_y.get(), family=self.initial_font.get())
        if self.grid_energy_var.get():
            self.ax1.grid()
        self.ax1.set_ylim(float(self.minE.get_name()), float(self.maxE.get_name()))
        self.ax1.set_xlim(0, 1)
        self.ax1.hlines(y=0, xmin=0, xmax=1, color="k", lw=1.5)
        self.ax1.tick_params(axis='x', which='major', labelsize=self.font_size_band_ticks.get())
        self.ax2.set_ylim(float(self.minE.get_name()), float(self.maxE.get_name()))
        self.ax2.set_xlim(0, 1)
        if self.label_DOS_var.get():
            self.ax2.set_xlabel('Density of States', fontsize=self.font_size_DOS_x.get(), family=self.initial_font.get())
        if self.label_energy_DOS_var.get():
            self.ax2.set_ylabel('Energy / eV', fontsize=self.font_size_DOS_number.get(), family=self.initial_font.get())
        if self.grid_DOS_var.get():
            self.ax2.grid()
        self.ax2.hlines(y=0, xmin=0, xmax=1, color='k', lw=1.5)
        self.ax2.tick_params(axis='x', which='major', labelsize=self.font_size_DOS_number.get())
        if not self.ticks_energy_var.get():
            self.ax1.set_yticklabels([])
        if not self.ticks_wavevector_var.get():
            self.ax1.set_xticklabels([])
        if not self.ticks_energy_DOS_var.get():
            self.ax2.set_yticklabels([])
        if not self.ticks_DOS_var.get():
            self.ax2.set_xticklabels([])

        toolbar_frame = Frame(self.parent) 
        toolbar_frame.grid(row=16,column=2,columnspan=4) 
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()


    def _from_rgb(self, rgb):
        """
        translates an rgb tuple of int to a tkinter friendly color code
        """
        return "#%02x%02x%02x" % rgb


    def clear(self):
        """
        Default values to start new project
        """
        self.initial_parameters()
        self.plot_widget.grid_forget()
        self.create_empty_plot()
        self.pDOS_E_var.set(False); self.pDOS_O_var.set(False)
        self.pDOS_E.config(state=DISABLED); self.pDOS_O.config(state=DISABLED)
        self.load_button.config(state=DISABLED); self.plot_button.config(state=DISABLED)


    def create_kpoint(self):
        """
        Create KPOINTS file
        """

        self.KPoints = Toplevel()
        self.KPoints.configure(bg = self._from_rgb((11, 165, 193)))
        self.KPoints.geometry("700x650")
        self.KPoints.iconbitmap('icon_band.ico')
        self.KPoints.grab_set()

        self.crystal_options = [
            'Triclinic, primitive, aP',
            'Monoclinic, primitive, mP',
            'Monoclinic, C-centered, mS',
            'Orthorhombic, primitive, oP',
            'Orthorhombic, face-centered, oF',
            'Orthorhombic, body-centered, oI',
            'Orthorhombic, C-centered, oS',
            'Tetragonal, primitive, tP',
            'Tetragonal, body-centered, tI',
            'Hexagonal, primitive, hP',
            'Rhombohedral, primitive, hR',
            'Cubic, primitive, cP',
            'Cubic, face-centered, cF',
            'Cubic, body-centered, cI'
        ]
        self.initial_crystal = StringVar(); self.initial_crystal.set(self.crystal_options[0])
        self.crystal_menu = OptionMenu(self.KPoints, self.initial_crystal, *self.crystal_options)
        self.crystal_menu.grid(row = 0, column = 1, padx = 10, pady = 10)
        self.crystal_label = Label(self.KPoints, text = 'Crystal System', relief = RIDGE, anchor = 'w')
        self.crystal_label.grid(row = 0, column = 0, padx = 10, pady = 10, ipadx = 32)
        self.lattice_a = EntryItem(self.KPoints, name = 'Lattice parameter a / Ang', row = 1)
        self.lattice_a.create_EntryItem(ipadx_label = 5); self.lattice_a.initial_val = DoubleVar(); self.lattice_a.set_name(1)
        self.lattice_b = EntryItem(self.KPoints, name = 'Lattice parameter b / Ang', row = 2)
        self.lattice_b.create_EntryItem(ipadx_label = 5); self.lattice_b.initial_val = DoubleVar(); self.lattice_b.set_name(1)
        self.lattice_c = EntryItem(self.KPoints, name = 'Lattice parameter c / Ang', row = 3)
        self.lattice_c.create_EntryItem(ipadx_label = 5); self.lattice_c.initial_val = DoubleVar(); self.lattice_c.set_name(1)
        self.lattice_alpha = EntryItem(self.KPoints, name = 'Alpha / deg', row = 4)
        self.lattice_alpha.create_EntryItem(ipadx_label = 41); self.lattice_alpha.initial_val = DoubleVar(); self.lattice_alpha.set_name(90)
        self.lattice_beta = EntryItem(self.KPoints, name = 'Beta / deg', row = 5)
        self.lattice_beta.create_EntryItem(ipadx_label = 45); self.lattice_beta.initial_val = DoubleVar(); self.lattice_beta.set_name(90)
        self.lattice_gamma = EntryItem(self.KPoints, name = 'Gamma / deg', row = 6)
        self.lattice_gamma.create_EntryItem(ipadx_label = 36); self.lattice_gamma.initial_val = DoubleVar(); self.lattice_gamma.set_name(90)

        self.minimum_d = EntryItem(self.KPoints, name = 'Minimum # of points', row = 7)
        self.minimum_d.create_EntryItem(ipadx_label = 17); self.minimum_d.initial_val = DoubleVar(); self.minimum_d.set_name(10)

        self.create_list_button = Button(self.KPoints, text='Create List of High-Symmetry Points', command=self.create_list)
        self.create_list_button.grid(row=8, column=0, columnspan=2, pady=10, ipadx=4)
        self.create_list_button['font'] = self.font_window

        self.create_kpoints_button = Button(self.KPoints, text='Create KPOINTS file', command=self.create_kpoints, state=DISABLED)
        self.create_kpoints_button.grid(row=9, column=0, columnspan=2, pady=10, ipadx=4)
        self.create_kpoints_button['font'] = self.font_window

        self.close_kpoints_button = Button(self.KPoints, text='Close', command=self.close_kpoints)
        self.close_kpoints_button.grid(row=10, column=0, columnspan=2, pady=10, ipadx=4)
        self.close_kpoints_button['font'] = self.font_window

        self.Check_path_button = list(); self.path_var = list()
        self.kpt_nmb = EntryItem(self.KPoints, 'Number of Kpts', column=4, state=DISABLED)


    def create_list(self):
        """
        Create paths between high-symmetry points to choose from
        """

        for c in self.Check_path_button:
            c.destroy()
        self.path_var = list()

        self.create_kpoints_button.config(state=NORMAL)

        self.kpts = K_points(
            ibrav=self.initial_crystal.get(),
            lattice_a=float(self.lattice_a.get_name()),
            lattice_b=float(self.lattice_b.get_name()),
            lattice_c=float(self.lattice_c.get_name()),
            alpha=float(self.lattice_alpha.get_name()),
            beta=float(self.lattice_beta.get_name()),
            gamma=float(self.lattice_gamma.get_name()),
            factor=int(self.minimum_d.get_name())
            )

        self.label_Path = Label(self.KPoints, text = 'KPoint_Path')
        self.label_Path.grid(row = 0, column = 3, columnspan=2, pady = (10, 5))
        self.label_Path['font'] = self.font_window

        Points, Path = self.kpts.Kpoint_path()

        nmb = 1
        for p in Path:
            self.path_var.append(BooleanVar())
            self.Check_path_button.append(Checkbutton(self.KPoints, text='{} --> {}'.format(p[0], p[1]), variable=self.path_var[-1]))
            self.Check_path_button[-1].grid(row=nmb, column=3, columnspan=2, ipadx = 20, pady=10); self.path_var[-1].set(True)
            nmb += 1

        self.List_data = self.kpts.min_distance(Path, Points)

        self.kpt_nmb.row = nmb
        self.kpt_nmb.create_EntryItem(ipadx_label=10); self.kpt_nmb.set_name(sum(self.List_data))

        for var in self.path_var:
            var.trace('w', self.click_path)


    def get_new_path(self):
        """
        Create a new path from chosen high-symmetry points
        """
        Points, Path = self.kpts.Kpoint_path()
        New_path = []
        for i in range(len(self.path_var)):
            if self.path_var[i].get() == True:
                New_path.append(Path[i])

        return Points, New_path


    def click_path(self, *args):
        """
        Update the number of kpoints depending on the chosen path
        """
        Points, New_path = self.get_new_path()

        self.kpt_nmb.set_name(sum(self.kpts.min_distance(New_path, Points)))


    def create_list_points(self, New_Path, Points, data):
        """
        Write list of kpoints for KPOINTS file
        Input:
        --------------------------
        New_path: list, shape (N, 2)
            List of chosen path between high-symmetry points, N is the number of paths
        Points: dictionary; keys: 'point' (Point of high-symmetry point), 'Sym' (Symbol of the high-symmetry point)
            Dictionary of high-symmetries point in the Brillouin zone
        data: List, shape (N)
            list of length between high-symmetry path
        """

        list_kpts = list()
        for p in range(len(data)):
            pt1 = np.array(Points[New_Path[p][0]]['point'])
            pt2 = np.array(Points[New_Path[p][1]]['point'])

            list_kpts.append('{:.4f} {:.4f} {:.4f}     1   {} \n'.format(*tuple(pt1), Points[New_Path[p][0]]['Sym']))
            for i in range(1, data[p] - 1):
                kpt_point = (pt2 - pt1) / (data[p] - 1) * i
                list_kpts.append('{:.4f} {:.4f} {:.4f}     1 \n'.format(*tuple(kpt_point)))

            list_kpts.append('{:.4f} {:.4f} {:.4f}     1   {} \n'.format(*tuple(pt2), Points[New_Path[p][1]]['Sym']))
            list_kpts.append('\n')

        return list_kpts


    def create_kpoints(self):
        """
        Create KPOINTS and POINTS.json file which is used to plot the electronic band structure
        """
        foldername_kpts = filedialog.askdirectory()
        Points, New_path = self.get_new_path()
        List_data = self.kpts.min_distance(New_path, Points)

        with open(foldername_kpts + '//POINTS.json', 'w') as fil:
            json.dump(List_data, fil)

        list_kpts = self.create_list_points(New_path, Points, List_data)
        with open(foldername_kpts + '//KPOINTS', 'w') as fil:
            fil.write('Electronic band structure of {} \n'.format(self.initial_crystal.get()))
            fil.write('{} \n'.format(sum(List_data)))
            fil.write('Reciprocal \n')
            for x in range(len(list_kpts)):
                fil.write(list_kpts[x])


    def close_kpoints(self):
        """
        Close the kpoint window
        """
        answer = messagebox.askokcancel('Exit?', 'Are you sure you want to exit?', icon='warning')

        if answer:
            self.KPoints.destroy()
        else:
            return


    def open_file(self):
        """
        Open folder which needs to include CONTCAR, KPOINTS, PROCAR_band, PROCAR_DOS, and POINTS.json files
        """
        self.clear()
        self.foldername = filedialog.askdirectory()

        if self.foldername == '':
            return

        if set(['CONTCAR', 'KPOINTS', 'PROCAR_band', 'PROCAR_DOS', 'POINTS.json']).issubset(set(os.listdir(self.foldername))):
            self.filename.set_name(self.foldername.split('/')[-1])
            self.create_empty_plot()

        else:
            messagebox.showerror(message = 'Folder needs to include CONTCAR, KPOINTS, PROCAR_band, PROCAR_DOS from VASP and POINTS.json calculated from this app! ' +
            'Please label them as stated.')
            return

        with open(self.foldername + '/CONTCAR') as con:
            contcar = con.readlines()

        self.cmp = contcar[5].split()

        self.load_button.config(state=NORMAL)


    def close_program(self):
        """
        Close the VASP program
        """
        self.parent.quit()
        self.parent.destroy()


    def plot_electronic_structure(self):
        """
        Plot the electronic band structure, connected to the Plot button
        """
        if self.foldername == '':
            return

        self.save_figure_button.config(state=NORMAL)
        self.save_fig_csv_button.config(state=NORMAL)
        self.plot_widget.grid_forget()
        self.plot()


    def get_energies(self):
        """
        Get information from PROCAR_DOS and PROCAR_band files; remove the valence band maximum from the energies
        """

        with open(self.foldername +"/PROCAR_DOS") as pro_DOS:
            procar_DOS = pro_DOS.readlines()

        self.DOS = Energy(procar_DOS)
        self.DOS.get_energies()

        with open(self.foldername +"/PROCAR_band") as pro_band:
            procar_band = pro_band.readlines()

        self.Band = Energy(procar_band)
        self.Band.get_energies()

        Eg, VBM = self.DOS.get_band_gap()

        self.Band.energy -= VBM
        self.DOS.energy -= VBM


    def get_contribution(self, energy, DOS_elements_new):
        """
        Get the contributions for each band and kpoint
        Input:
        ----------------------
        energy: ndarray, shape (M, N), dtype=float
            Array of energies for M bands and N kpoints
        DOS_element_new: ndarray, shape (M, N, X), dtype=float
            Array of DOS for each band and kpoint as well as element (orbital)
        """

        contrib = np.zeros((len(energy), len(energy[0]), len(DOS_elements_new[0][0])))
        for i in range(len(energy)):

            for j in range(len(energy[0])):
                sume = 0
                for c in range(len(DOS_elements_new[i][j])):
                    sume += DOS_elements_new[i][j][c]**2
                total_DOS = np.sqrt(sume)
                for k in range(len(DOS_elements_new[0][0])):
                    contrib[i, j, k] = DOS_elements_new[i][j][k] / total_DOS

        return contrib


    def get_kpoints(self):
        """
        Get ticks and number of kpoints between high-symmetry points
        """

        with open(self.foldername + "/Points.json") as json_file:
            Kpoint_mesh = json.load(json_file)

        with open(self.foldername + "/KPOINTS") as k:
            kpoints = k.readlines()

        self.ticks, self.distance = self.Band.get_distance(Kpoint_mesh, kpoints)


    def Edit_graph(self):
        """
        Edit plot by changing fonts, font size, labels, or grids
        """

        self.Top = Toplevel()
        self.Top.configure(bg = self._from_rgb((11, 165, 193)))
        self.Top.geometry("800x650")
        self.Top.iconbitmap('icon_band.ico')
        self.Top.grab_set()

        self.label_input = Label(self.Top, text = 'Labels/Ticks')
        self.label_input.grid(row = 0, column = 0, pady = (10, 5))
        self.label_input['font'] = self.font_window

        self.label_input = Label(self.Top, text = 'Figure')
        self.label_input.grid(row = 0, column = 2, pady = (10, 5))
        self.label_input['font'] = self.font_window

        self.font_size_band_x_entry = Entry(self.Top, textvariable = self.font_size_band_x, width = 24)
        self.font_size_band_x_entry.grid(row = 1, column = 1, padx = 10, pady=10)
        self.font_size_band_x_label = Label(self.Top, text = 'Font Size Bandstructure x', relief = RIDGE, anchor = 'w')
        self.font_size_band_x_label.grid(row = 1, column = 0, padx = 10, pady=10, ipadx = 14)

        self.font_size_band_y_entry = Entry(self.Top, textvariable = self.font_size_band_y, width = 24)
        self.font_size_band_y_entry.grid(row = 2, column = 1, padx = 10)
        self.font_size_band_y_label = Label(self.Top, text = 'Font Size Bandstructure y', relief = RIDGE, anchor = 'w')
        self.font_size_band_y_label.grid(row = 2, column = 0, padx = 10, ipadx = 14)

        self.font_size_band_ticks_entry = Entry(self.Top, textvariable = self.font_size_band_ticks, width = 24)
        self.font_size_band_ticks_entry.grid(row = 3, column = 1, padx = 10, pady=10)
        self.font_size_band_ticks_label = Label(self.Top, text = 'Font Size Bandstructure ticks x', relief = RIDGE, anchor = 'w')
        self.font_size_band_ticks_label.grid(row = 3, column = 0, padx = 10, pady=10, ipadx = 0)

        self.font_size_band_energy_entry = Entry(self.Top, textvariable = self.font_size_band_energy, width = 24)
        self.font_size_band_energy_entry.grid(row = 4, column = 1, padx = 10, pady=10)
        self.font_size_band_energy_label = Label(self.Top, text = 'Font Size Bandstructure ticks y', relief = RIDGE, anchor = 'w')
        self.font_size_band_energy_label.grid(row = 4, column = 0, padx = 10, pady=10, ipadx = 0)

        self.font_size_DOS_x_entry = Entry(self.Top, textvariable = self.font_size_DOS_x, width = 24)
        self.font_size_DOS_x_entry.grid(row = 5, column = 1, padx = 10, pady=10)
        self.font_size_DOS_x_label = Label(self.Top, text = 'Font Size DOS x', relief = RIDGE, anchor = 'w')
        self.font_size_DOS_x_label.grid(row = 5, column = 0, padx = 10, pady=10, ipadx = 39)

        self.font_size_DOS_y_entry = Entry(self.Top, textvariable = self.font_size_DOS_y, width = 24)
        self.font_size_DOS_y_entry.grid(row = 6, column = 1, padx = 10, pady=10)
        self.font_size_DOS_y_label = Label(self.Top, text = 'Font Size DOS y', relief = RIDGE, anchor = 'w')
        self.font_size_DOS_y_label.grid(row = 6, column = 0, padx = 10, pady=10, ipadx = 39)

        self.font_size_DOS_ticks_entry = Entry(self.Top, textvariable = self.font_size_DOS_number, width = 24)
        self.font_size_DOS_ticks_entry.grid(row = 7, column = 1, padx = 10, pady=10)
        self.font_size_DOS_ticks_label = Label(self.Top, text = 'Font Size DOS ticks', relief = RIDGE, anchor = 'w')
        self.font_size_DOS_ticks_label.grid(row = 7, column = 0, padx = 10, pady=10, ipadx = 31)

        self.size_x_entry = Entry(self.Top, textvariable = self.size_x, width = 24)
        self.size_x_entry.grid(row = 1, column = 3, padx = 10, pady=10)
        self.size_x_label = Label(self.Top, text = 'Figure Width', relief = RIDGE, anchor = 'w')
        self.size_x_label.grid(row = 1, column = 2, padx = 10, pady=10, ipadx = 37)

        self.size_y_entry = Entry(self.Top, textvariable = self.size_y, width = 24)
        self.size_y_entry.grid(row = 2, column = 3, padx = 10, pady=10)
        self.size_y_label = Label(self.Top, text = 'Figure Height', relief = RIDGE, anchor = 'w')
        self.size_y_label.grid(row = 2, column = 2, padx = 10, pady=10, ipadx = 35)

        self.size_x_space_entry = Entry(self.Top, textvariable = self.size_x_space, width = 24)
        self.size_x_space_entry.grid(row = 3, column = 3, padx = 10, pady=10)
        self.size_x_space_label = Label(self.Top, text = 'Figure Move Horizontal', relief = RIDGE, anchor = 'w')
        self.size_x_space_label.grid(row = 3, column = 2, padx = 10, pady=10, ipadx=8)

        self.size_y_space_entry = Entry(self.Top, textvariable = self.size_y_space, width = 24)
        self.size_y_space_entry.grid(row = 4, column = 3, padx = 10, pady=10)
        self.size_y_space_label = Label(self.Top, text = 'Figure Move Vertical', relief = RIDGE, anchor = 'w')
        self.size_y_space_label.grid(row = 4, column = 2, padx = 10, pady=10, ipadx=17)

        self.size_x_length_entry = Entry(self.Top, textvariable = self.size_x_length, width = 24)
        self.size_x_length_entry.grid(row = 5, column = 3, padx = 10, pady=10)
        self.size_x_length_label = Label(self.Top, text = 'Figure Length Horizontal', relief = RIDGE, anchor = 'w')
        self.size_x_length_label.grid(row = 5, column = 2, padx = 10, pady=10, ipadx=5)

        self.size_y_length_entry = Entry(self.Top, textvariable = self.size_y_length, width = 24)
        self.size_y_length_entry.grid(row = 6, column = 3, padx = 10, pady=10)
        self.size_y_length_label = Label(self.Top, text = 'Figure Length Vertical', relief = RIDGE, anchor = 'w')
        self.size_y_length_label.grid(row = 6, column = 2, padx = 10, pady=10, ipadx=14)

        self.font_menu = OptionMenu(self.Top, self.initial_font, *self.font_options)
        self.font_menu.grid(row = 8, column = 1, padx = 10, pady = 10)
        self.font_label = Label(self.Top, text = 'Font', relief = RIDGE, anchor = 'w')
        self.font_label.grid(row = 8, column = 0, padx = 10, pady = 10, ipadx = 32)

        self.color_2plot_menu = OptionMenu(self.Top, self.initial_color_2plot, *self.color_2plot_options)
        self.color_2plot_menu.grid(row = 8, column = 3, padx = 10, pady = 10)
        self.color_2plot_label = Label(self.Top, text = 'Color for 2 Elements', relief = RIDGE, anchor = 'w')
        self.color_2plot_label.grid(row = 8, column = 2, padx = 10, pady = 10, ipadx = 2)

        self.label_energy = Checkbutton(self.Top, text='Energy Label Band Structure', variable=self.label_energy_var)
        self.label_energy.grid(row=9, column=0, columnspan=2, pady=10, ipadx=5)

        self.label_ticks = Checkbutton(self.Top, text='Wavevector Label', variable=self.label_ticks_var)
        self.label_ticks.grid(row=9, column=2, pady=10, ipadx=18)

        self.label_energy_DOS = Checkbutton(self.Top, text='Energy Label DOS', variable=self.label_energy_DOS_var)
        self.label_energy_DOS.grid(row=10, column=0, columnspan=2, pady=10, ipadx=33)

        self.label_DOS = Checkbutton(self.Top, text='DOS Label', variable=self.label_DOS_var)
        self.label_DOS.grid(row=10, column=2, pady=10, ipadx=37)

        self.grid_energy = Checkbutton(self.Top, text='Grid Band Structure', variable=self.grid_energy_var)
        self.grid_energy.grid(row=11, column=0, columnspan=2, pady=10, ipadx=29)

        self.grid_DOS = Checkbutton(self.Top, text='Grid DOS', variable=self.grid_DOS_var)
        self.grid_DOS.grid(row=11, column=2, pady=10, ipadx=40)

        self.ticks_energy = Checkbutton(self.Top, text='Ticks Energy Band', variable=self.ticks_energy_var)
        self.ticks_energy.grid(row=9, column=3, pady=10, ipadx=20)

        self.ticks_wavevector = Checkbutton(self.Top, text='Ticks Wavevector', variable=self.ticks_wavevector_var)
        self.ticks_wavevector.grid(row=10, column=3, pady=10, ipadx=22)

        self.ticks_energy_DOS = Checkbutton(self.Top, text='Ticks Energy DOS', variable=self.ticks_energy_DOS_var)
        self.ticks_energy_DOS.grid(row=11, column=3, pady=10, ipadx=22)

        self.ticks_DOS = Checkbutton(self.Top, text='Ticks DOS', variable=self.ticks_DOS_var)
        self.ticks_DOS.grid(row=12, column=3, pady=10, ipadx=42)

        btn_color = Button(self.Top, text = 'Color Band Structure', command = self.choose_color)
        btn_color.grid(row = 7, column = 2, padx =10, pady = 10, ipadx = 13)

        self.frame = Frame(self.Top, border=1, relief=SUNKEN, width=60, height=30)
        self.frame.grid(row=7, column=3, padx=10, pady=10, ipadx=10)
        self.frame.config(bg=self.hx)

        btn_close = Button(self.Top, text = 'Save/Close', command = self.close_update_graph)
        btn_close.grid(row = 13, column = 3, padx =10, pady = 10, ipadx = 35)

        btn_default = Button(self.Top, text = 'Save as Default', command = self.default)
        btn_default.grid(row = 12, column = 0, padx =10, pady = 10, ipadx = 35)


    def choose_color(self):
        """
        Choose a color for the electronic band structure
        """
        (color_new, self.hx) = colorchooser.askcolor()
        if color_new != None:
            self.color = np.array(color_new) / 256.
        self.frame.config(bg=self.hx)


    def close_update_graph(self):
        """
        Close Edit window and update the plot
        """

        self.plot_widget.grid_forget()

        if self.plot_button['state'] == NORMAL:
            self.plot()

        else:
            self.create_empty_plot()
        self.Top.destroy()


    def default(self):
        """
        Produce a default file which is used for the program
        """

        dic = {
            'font' : self.initial_font.get(),
            'color_2plot' : self.initial_color_2plot.get(),
            'size_band_x' : self.font_size_band_x.get(),
            'size_band_y' : self.font_size_band_y.get(),
            'size_band_ticks' : self.font_size_band_ticks.get(),
            'size_band_energy' : self.font_size_band_energy.get(),
            'size_DOS_x' : self.font_size_DOS_x.get(),
            'size_DOS_y' : self.font_size_DOS_y.get(),
            'size_DOS_ticks' : self.font_size_DOS_number.get(),
            'figure_size_x' : self.size_x.get(),
            'figure_size_y' : self.size_y.get(),
            'figure_space_x' : self.size_x_space.get(),
            'figure_space_y' : self.size_y_space.get(),
            'figure_length_x' : self.size_x_length.get(),
            'figure_length_y' : self.size_y_length.get(),
            'label_energy' : self.label_energy_var.get(),
            'label_energy_DOS' : self.label_energy_DOS_var.get(),
            'label_ticks' : self.label_ticks_var.get(),
            'label_DOS' : self.label_DOS_var.get(),
            'grid_energy' : self.grid_energy_var.get(),
            'grid_DOS' : self.grid_DOS_var.get(),
            'ticks_energy' : self.ticks_energy_var.get(),
            'ticks_wavevector' : self.ticks_wavevector_var.get(),
            'ticks_energy_DOS' : self.ticks_energy_DOS_var.get(),
            'ticks_DOS' : self.ticks_DOS_var.get(),
            'color' : self.color.tolist(),
            'color_hex' : self.hx,
        }

        with open('~default.json', 'w') as d:
            json.dump(dic, d)


    def sum_DOS_elements(self):
        """
        Sum over the ions to get DOS for one element
        """

        with open(self.foldername + '/CONTCAR') as con:
            contcar = con.readlines()

        self.Band.element_DOS(contcar)
        self.DOS.element_DOS(contcar)


    def rgbline(self, ax, k, e, red, green, blue, alpha=1.):
        """
        Produce segments for rgb values
        Input:
        -------------------------
        ax: Matplotlib subplot
        k: ndarray, shape (N), dtype=float
            Array of kpoints
        e: ndarray, shape (M, N), dtype=float
            Array of energies for M bands and N kpoints
        red: ndarray, shape (N)
            Array of red contribution
        green: ndarray, shape (N)
            Array of green contribution
        blue: ndarray, shape (N)
            Array of blue contribution
        alpha: ndarray, shape (N)
            Array of shading
        http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
        """

        pts = np.array([k, e]).T.reshape(-1, 1, 2)

        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
        nseg = len(k) - 1
        if len(red) == 0:
            r = [0 for i in range(nseg)]
        else:
            r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]

        if len(green) == 0:
            g = [0 for i in range(nseg)]
        else:
            g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]

        if len(blue) == 0:
            b = [0 for i in range(nseg)]
        else:
            b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]

        a = np.ones(nseg, np.float64) * alpha
        lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
        ax.add_collection(lc)


    def colors(self):
        """
        List of potential colors [in this order]
        """
        return  ['red',
            'green',
            'blue',
            'orange',
            'cyan',
            'yellow',
            'lawngreen',
            'pink',
            'magenta',
            'navy',
            'springgreen'
        ]


    def plot_pDOS(self, ax):
        """
        Plot the elemental projected DOS
        Input:
        ---------------------------------
        ax: subplot of Matplotlib
        """
        colormap = self.colors()
        if self.pDOS_E_var.get():

            if len(self.cmp) == 2:

                if self.initial_color_2plot.get() == 'red-blue':
                    colormap = ['red', 'blue']

                elif self.initial_color_2plot.get() == 'green-blue':
                    colormap = ['green', 'blue']

        for cmpd in range(len(self.partial_DOS)):
            ax.plot(self.Energy_DOS, -self.partial_DOS[cmpd],

                    color=colormap[cmpd], label=self.cmp[cmpd], lw=2)
        if self.label_DOS_var.get():
            ax.set_xlabel('Projected DOS', fontsize=self.font_size_DOS_x.get(), family=self.initial_font.get())


    def plot_oDOS(self, ax):
        """
        Plot the orbital projected DOS
        Input:
        -------------------------------
        ax: subplot of Matplotlib
        """

        colormap = self.colors()
        orbital_map = ['s', 'p', 'd', 'f', 'g']

        for k in range(len(self.orbital_DOS)):
            ax.plot(self.Energy_DOS, -self.orbital_DOS[k],

                    color=colormap[k], label=orbital_map[k], lw=2)

        if self.label_DOS_var.get():
            ax.set_xlabel('Projected DOS', fontsize=self.font_size_DOS_x.get(), family=self.initial_font.get())


    def load_electronic_properties(self):
        """
        Load the electronic properties from PROCAR_band and PROCAR_DOS and compute the DOS over the entire Brillouin zone
        Connected to the Load button
        """
        self.list_energy = [self.minE.get_name(), self.maxE.get_name(), self.Eres.get_name()]

        self.get_energies()
        self.get_kpoints()
        self.sum_DOS_elements()

        self.DOS.energy_DOS, self.DOS.totDOS_DOS = self.DOS.sum_partial_DOS(self.DOS.totDOS, float(self.minE.get_name()), float(self.maxE.get_name()), float(self.Eres.get_name()))

        self.contrib = self.get_contribution(self.Band.energy, self.Band.DOS_element_new)
        self.contrib_orbital = self.get_contribution(self.Band.energy, self.Band.DOS_orbitals)

        self.orbital_DOS = []

        for cmpd in range(len(self.DOS.DOS_orbitals[0][0])):
            self.Energy_DOS, oDOS = self.DOS.sum_partial_DOS(np.moveaxis(self.DOS.DOS_orbitals, 2, 0)[cmpd], float(self.minE.get_name()), float(self.maxE.get_name()), float(self.Eres.get_name()))
            self.orbital_DOS.append(oDOS)


        self.partial_DOS = []

        for cmpd in range(len(self.DOS.DOS_element_new[0][0])):
            self.Energy_DOS, pDOS = self.DOS.sum_partial_DOS(np.moveaxis(self.DOS.DOS_element_new, 2, 0)[cmpd], float(self.minE.get_name()), float(self.maxE.get_name()), float(self.Eres.get_name()))
            self.partial_DOS.append(pDOS)

        self.plot_button.config(state=NORMAL)

        if len(self.cmp) in [2, 3]:
            self.pDOS_E.config(state=NORMAL)

        if len(self.partial_DOS) in [2, 3]:
            self.pDOS_O.config(state=NORMAL)


    def plot(self, save=False, filename=''):
        """
        Plot electronic band structure and DOS
        Input:
        --------------------------------
        save: Boolean
            If save is True, the figure will be saved (default is False)
        filename: str
            Filename from Save button
        """

        if self.list_energy != [self.minE.get_name(), self.maxE.get_name(), self.Eres.get_name()]:
            load_new = messagebox.askyesno('New energy range!',
             'You have changed minimum energy, maximum energy, or the resolution. You need to load first the new data. Do you want to load the data?')

            if load_new:
                self.load_electronic_properties()

        plt.rcParams["font.family"] = self.initial_font.get()
        plt.rcParams.update({'font.size': self.font_size_band_energy.get()})

        self.fig = Figure(figsize= (self.size_x.get(), self.size_y.get()), dpi = 100)

        gs = self.fig.add_gridspec(1, 2, width_ratios=[2, 1,])
        gs.update(left=0.1, right=0.95, wspace=0.15)
        self.ax1 = self.fig.add_subplot(gs[0])
        self.ax2 = self.fig.add_subplot(gs[1])

        self.ax1.set_xlim(0, 1.)
        self.ax1.set_ylim(float(self.minE.get_name()), float(self.maxE.get_name()))
        self.ax2.set_xlim(float(self.minE.get_name()), float(self.maxE.get_name()))
        self.ax2.set_ylim(-0.0005, float(self.ymax.var.get()))
        if self.grid_energy_var.get():
            self.ax1.grid()
        if self.label_ticks_var.get():
            self.ax1.set_xlabel('Wavevector $k$', fontsize=self.font_size_band_x.get(), family=self.initial_font.get())
        if self.label_energy_var.get():
            self.ax1.set_ylabel('$E-E_F$ / eV', fontsize=self.font_size_band_y.get(), family=self.initial_font.get())

        for p in self.distance:
            self.ax1.plot([p, p], [float(self.minE.get_name()), float(self.maxE.get_name())], 'k-', color='grey')
        self.ax1.set_xticks(self.distance)
        self.ax1.set_xticklabels(self.ticks)
        self.ax1.tick_params(axis='x', which='major', labelsize=self.font_size_band_ticks.get())
        self.ax2.tick_params(axis='x', which='major', labelsize=self.font_size_DOS_number.get())

        self.ax2.fill_between(self.DOS.energy_DOS, -self.DOS.totDOS_DOS, 0, color=(0.7, 0.7, 0.7), facecolor=(0.7, 0.7, 0.7))
        self.ax2.plot(self.DOS.energy_DOS, -self.DOS.totDOS_DOS, color =(0.6, 0.6, 0.6), label='Total DOS')
        if self.label_DOS_var.get():
            self.ax2.set_xlabel('Density of States', fontsize=self.font_size_DOS_y.get(), family=self.initial_font.get())
        if self.label_energy_DOS_var.get():
            self.ax2.set_ylabel('$E-E_F$ / eV', fontsize=self.font_size_DOS_y.get(), family=self.initial_font.get())

        if self.pDOS_E_var.get() or self.pDOS_O_var.get():
            gs2 = self.fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 2,], height_ratios=[1, 4,])
            gs2.update(bottom=0.6, top=0.95)
            self.ax3 = self.fig.add_subplot(gs2[1]); self.ax3.axis('off')
            self.ax4 = self.fig.add_subplot(gs2[0]); self.ax4.axis('off')
            self.ax5 = self.fig.add_subplot(gs2[2]); self.ax5.axis('off')

            if len(self.cmp) == 3 and self.pDOS_E_var.get():
                rgb_triangle = plt.imread('rgb_triangle.png')
                self.ax3.imshow(rgb_triangle)
                self.ax3.text(200, 0, self.cmp[0], color='red')
                self.ax4.text(33, 0, self.cmp[2], color='blue'); self.ax4.set_xlim(-100, 2)
                self.ax5.text(-0.45, 0, self.cmp[1], color='green'); self.ax5.set_xlim(0, 1)
                for b in range(len(self.Band.energy)):
                    self.rgbline(self.ax1,
                        self.Band.kpts,
                        self.Band.energy[b],
                        self.contrib[b, :, 0],
                        self.contrib[b, :, 1],
                        self.contrib[b, :, 2])

            if len(self.Band.DOS_orbitals[0][0]) == 3 and self.pDOS_O_var.get():
                rgb_triangle = plt.imread('rgb_triangle.png')
                self.ax3.imshow(rgb_triangle)
                self.ax3.text(290, 0, 's', color='red')
                self.ax4.text(39, 0, 'd', color='blue'); self.ax4.set_xlim(-100, 2)
                self.ax5.text(-0.45, 0, 'p', color='green'); self.ax5.set_xlim(0, 1)
                for b in range(len(self.Band.energy)):
                    self.rgbline(self.ax1,
                        self.Band.kpts,
                        self.Band.energy[b],
                        self.contrib_orbital[b, :, 0],
                        self.contrib_orbital[b, :, 1],
                        self.contrib_orbital[b, :, 2])

            elif len(self.Band.DOS_orbitals[0][0]) == 2 and self.pDOS_O_var.get():

                if self.initial_color_2plot.get() == 'red-green':
                    rg_line = plt.imread('rg_line.png')
                    self.ax3.imshow(rg_line)
                    self.ax4.text(0, 0, 's', color='red')
                    self.ax5.text(0, 0, 's', color='green')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            self.contrib_orbital[b, :, 0],
                            self.contrib_orbital[b, :, 1],
                            [])

                elif self.initial_color_2plot.get() == 'red-blue':
                    rb_line = plt.imread('rb_line.png')
                    self.ax3.imshow(rb_line)
                    self.ax4.text(0, 0, 's', color='red')
                    self.ax5.text(0, 0, 'p', color='blue')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            self.contrib_orbital[b, :, 0],
                            [],
                            self.contrib_orbital[b, :, 1],
                            )

                elif self.initial_color_2plot.get() == 'green-blue':
                    gb_line = plt.imread('gb_line.png')
                    self.ax3.imshow(gb_line)
                    self.ax4.text(0, 0, 's', color='green')
                    self.ax5.text(0, 0, 'p', color='blue')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            [],
                            self.contrib_orbital[b, :, 0],
                            self.contrib_orbital[b, :, 1],
                            )


            elif len(self.cmp) == 2 and self.pDOS_E_var.get():

                if self.initial_color_2plot.get() == 'red-green':
                    rg_line = plt.imread('rg_line.png')
                    self.ax3.imshow(rg_line)
                    self.ax4.text(0, 0, self.cmp[0], color='red')
                    self.ax5.text(0, 0, self.cmp[1], color='green')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            self.contrib[b, :, 0],
                            self.contrib[b, :, 1],
                            [])

                elif self.initial_color_2plot.get() == 'red-blue':
                    rb_line = plt.imread('rb_line.png')
                    self.ax3.imshow(rb_line)
                    self.ax4.text(0, 0, self.cmp[0], color='red')
                    self.ax5.text(0, 0, self.cmp[1], color='blue')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            self.contrib[b, :, 0],
                            [],
                            self.contrib[b, :, 1],
                            )

                elif self.initial_color_2plot.get() == 'green-blue':
                    gb_line = plt.imread('gb_line.png')
                    self.ax3.imshow(gb_line)
                    self.ax4.text(0, 0, self.cmp[0], color='green')
                    self.ax5.text(0, 0, self.cmp[1], color='blue')

                    for b in range(len(self.Band.energy)):
                        self.rgbline(self.ax1,
                            self.Band.kpts,
                            self.Band.energy[b],
                            [],
                            self.contrib[b, :, 0],
                            self.contrib[b, :, 1],
                            )

                self.ax3.set_xlim(0, 500); self.ax3.set_ylim(200, 240)
                self.ax4.set_ylim(-1, 1); self.ax5.set_ylim(-1, 1); self.ax4.set_xlim(-100, 2)

        else:
            for b in range(len(self.Band.energy)):
                self.ax1.plot(
                    self.Band.kpts,
                    self.Band.energy[b],
                    color=self.color)

        if self.pDOS.get() == 2:
            self.plot_pDOS(self.ax2)

        elif self.pDOS.get() == 3:
            self.plot_oDOS(self.ax2)

        self.ax1.hlines(y=0, xmin=0, xmax=1, color="k", lw=2)
        r = Affine2D().rotate_deg(90)

        for x in self.ax2.images + self.ax2.lines + self.ax2.collections:
            trans = x.get_transform()
            x.set_transform(r + trans)
            if isinstance(x, PathCollection):
                transoff = x.get_offset_transform()
                x._transOffset = r + transoff

        old = self.ax2.axis()
        self.ax2.axis(old[2:4] + old[0:2])

        if self.grid_DOS_var.get():
            self.ax2.grid()
        if not self.ticks_energy_var.get():
            self.ax1.set_yticklabels([])
        if not self.ticks_wavevector_var.get():
            self.ax1.set_xticklabels([])
        if not self.ticks_energy_DOS_var.get():
            self.ax2.set_yticklabels([])
        if not self.ticks_DOS_var.get():
            self.ax2.set_xticklabels([])
        self.ax2.hlines(y=0, xmin=-0.5, xmax=float(self.ymax.get_name()) + 0.5, color='k', lw =2)
        self.ax2.legend(fancybox=True, shadow=True, prop={'size': 18})

        if save:
            self.fig.set_size_inches(12, 8)
            self.fig.savefig(filename, dpi=int(self.set_dpi.get_name()))

        else:
            self.canvas = FigureCanvasTkAgg(self.fig, master = self.parent)
            self.canvas.draw()
            self.plot_widget = self.canvas.get_tk_widget()
            self.plot_widget.grid(row = 1, column = 3, columnspan = 11, rowspan = 13)

            toolbar_frame = Frame(self.parent) 
            toolbar_frame.grid(row=16,column=2,columnspan=4) 
            toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
            toolbar.update()


    def save_electronic_structure(self):
        """
        Save the Figure as shown in the Main window
        """

        filename = filedialog.asksaveasfilename(title = 'Save file', defaultextension="*.*", filetypes = [('Portable Network Graphics', '*.png'), ('PDF' ,'*.pdf'), ('Encapsulated Postscript', '*.eps')])
        if filename == '':
            return

        self.plot(save=True, filename=filename)


    def save_csv_file(self):
        """
        Save the DOS and bands as a .csv file
        """

        filename = filedialog.asksaveasfilename(title='Save file', filetypes=[('CSV (Comma delimited)', '*.csv')])
        if filename == '':
            return

        with open(filename + '_DOS.csv', 'w') as fil:
            fil.write('Energy / eV, Total DOS')

            for p in range(len(self.cmp)):
                fil.write(" , {} ".format(self.cmp[p]))

            orbitals = ['s', 'p', 'd', 'f', 'g']
            for o in range(len(self.orbital_DOS)):
                fil.write(" , {}".format(orbitals[o]))
            fil.write('\n')

            for i in range(len(self.DOS.energy_DOS)):
                fil.write("{}, {}".format(self.DOS.energy_DOS[i], self.DOS.totDOS_DOS[i]))

                for p in range(len(self.cmp)):
                    fil.write(" , {}".format(self.partial_DOS[p][i]))

                for o in range(len(self.orbital_DOS)):
                    fil.write(" , {} ".format(self.orbital_DOS[o][i]))
                fil.write("\n")

        with open(filename + '_Band.csv', 'w') as fil:
            fil.write('ticks, k-points, Energy / eV \n')

            counter = 0
            for k in range(len(self.Band.kpts)):
                if k + 1 < len(self.Band.kpts):
                    if self.distance[counter] < self.Band.kpts[k + 1]:
                        fil.write(self.ticks[counter])
                        counter += 1
                    else:
                        fil.write('')
                else:
                    fil.write(self.ticks[counter])

                fil.write(" , {}".format(self.Band.kpts[k]))
                for b in range(len(self.Band.energy)):
                    fil.write(" , {} ".format(self.Band.energy[b][k]))
                fil.write("\n")


    def welcome(self):
        """
        Create a welcome window in the Help menu
        """

        welcome = Help()
        welcome.welcome()


    def documentary(self):
        """
        Create a documentary in the Help menu
        """

        documentary = Help()
        documentary.documentary()


    def about(self):
        """
        Create an about window in the Help menu
        """

        about = Help()
        about.screen_help.geometry('700x450')
        about.about()


if __name__ == "__main__":
    root = Tk()
    MainApplication(root)
    root.mainloop()
