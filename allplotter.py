"""
TODO
----
freq. response calculations
refresh lines without refreshing whole plot
export df's
"""
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import make_interp_spline, BSpline #interpolation
from scipy.interpolate import griddata as gd
from tkinter import *
from tkinter import ttk
from PIL import Image, ImageTk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

DEBUG = 0
DEFAULT_VALVE = 'E88CC'
DEFAULT_VSUPPLY = 265
DEFAULT_Ra = 33000
DEFAULT_Rk = 470
DEFAULT_INPUTSIGNAL = 2000
DEFAULT_XMAX = 300
DEFAULT_YMAX = 10
DEFAULT_RL = 1000000
DEFAULT_CO = 22
#DEFAULT_CK = 1
DEFAULT_RG = 68000
#DEFAULT_CF = 0

class mclass:
    def __init__(self,  window):
        self.window = window
        self.fig = Figure(figsize=(11,8))
        self.ax = self.fig.add_subplot(111)

        # read valve specs and grid lines from csv
        self.read_valvedata()
        self.read_valvespecs()

        # build UI
        window.columnconfigure(0, weight=1)
        window.columnconfigure(1, weight=1)
        window.columnconfigure(2, weight=3)
        window.columnconfigure(3, weight=1)
        window.columnconfigure(4, weight=1)

        # title
        self.title = Label(window, text='andmarti Load Line Plotter', fg='#1C5AAC', font=('Courier New', 24, 'bold'))
        self.title.grid(row=0, column=0, columnspan=5, padx=10, sticky=N)

        # circuit
        io = Image.open('images/circuit.jpg')
        self.img = ImageTk.PhotoImage(io)
        self.circuit = Label(window, image=self.img)
        self.circuit.grid(row=6, column=2, rowspan=10, sticky=N)

        # plot
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().grid(row=1, column=4, rowspan=44, columnspan=3, sticky=W, padx=0, pady=0)
        self.canvas.mpl_connect('motion_notify_event', self.motion_hover)
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7, alpha=0.5)
        self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center', labelpad=20)
        self.ax.set_xlabel('Va, V', fontsize=16, loc='center',  labelpad=20)

        # valve
        self.lbl_valve = Label(window, text="VALVE", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_valve.grid(row=1, column=0, sticky=W, padx=50, pady=5)
        self.str_valve = StringVar()
        self.str_valve.set(DEFAULT_VALVE)
        self.cmb_valve = ttk.Combobox(window, values=self.df['valve'].drop_duplicates().values.tolist(), textvariable=self.str_valve, font=('Courier New', 12), width=13)
        self.cmb_valve['state'] = 'readonly'
        self.cmb_valve.grid(row=1, column=1, sticky=W, padx=5, pady=5)
        self.cmb_valve.bind('<<ComboboxSelected>>', self.valve_changed)

        # Vs
        self.lbl_supply = Label(window, text="Vsupply, V", font=('Courier New', 12), background=self.window['bg'], width=20, anchor='w')
        self.lbl_supply.grid(row=2, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_supply = StringVar()
        self.str_supply.set(DEFAULT_VSUPPLY)
        self.etr_supply = Entry(window, textvariable=self.str_supply, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_supply.grid(row=2, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Ra = Label(window, text="Ra, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Ra.grid(row=3, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Ra = StringVar()
        self.str_Ra.set(DEFAULT_Ra)
        self.etr_Ra = Entry(window, textvariable=self.str_Ra, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_Ra.grid(row=3, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Rk = Label(window, text="Rk, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Rk.grid(row=4, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Rk = StringVar()
        self.str_Rk.set(DEFAULT_Rk)
        self.etr_Rk = Entry(window, textvariable=self.str_Rk, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_Rk.grid(row=4, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Vq = Label(window, text="Vq, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Vq.grid(row=5, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Vq = StringVar()
        self.lbl_Vqval = Label(window, textvariable=self.str_Vq, font=('Courier New', 12), width=15, anchor="w")
        self.lbl_Vqval.grid(row=5, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Iq = Label(window, text="Iq, mA", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Iq.grid(row=6, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Iq = StringVar()
        self.lbl_Iqval = Label(window, textvariable=self.str_Iq, font=('Courier New', 12), width=15, anchor="w")
        self.lbl_Iqval.grid(row=6, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Vgk = Label(window, text="Vgk, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Vgk.grid(row=7, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Vgk = StringVar()
        self.lbl_Vgkval = Label(window, textvariable=self.str_Vgk, font=('Courier New', 12), width=15, anchor="w")
        self.lbl_Vgkval.grid(row=7, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Xmax = Label(window, text="Xmax, V", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Xmax.grid(row=8, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Xmax = StringVar()
        self.str_Xmax.set(DEFAULT_XMAX)
        self.etr_Xmax = Entry(window, textvariable=self.str_Xmax, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_Xmax.grid(row=8, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_Ymax = Label(window, text="Ymax, mA", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_Ymax.grid(row=9, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_Ymax = StringVar()
        self.str_Ymax.set(DEFAULT_YMAX)
        self.etr_Ymax = Entry(window, textvariable=self.str_Ymax, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_Ymax.grid(row=9, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_ck = Label(window, text="Ck, uF", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_ck.grid(row=10, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_ck = StringVar()
        #self.str_ck.set(DEFAULT_CK)
        self.etr_ck = Entry(window, textvariable=self.str_ck, font=('Courier New', 12), width=15)
        self.etr_ck.grid(row=10, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        self.lbl_rl = Label(window, text="Rl, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_rl.grid(row=11, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_rl = StringVar()
        self.str_rl.set(DEFAULT_RL)
        self.etr_rl = Entry(window, textvariable=self.str_rl, font=('Courier New', 12), width=15)
        self.etr_rl.grid(row=11, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        # input signal
        self.lbl_inputsignal = Label(window, text="input signal, mVpp", font=('Courier New', 12), background=self.window['bg'], width=18)
        self.lbl_inputsignal.grid(row=12, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_inputsignal = StringVar()
        self.str_inputsignal.set(DEFAULT_INPUTSIGNAL)
        self.etr_inputsignal = Entry(window, textvariable=self.str_inputsignal, font=('Courier New', 12), width=15, background='misty rose')
        self.etr_inputsignal.grid(row=12, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        # Rg
        self.lbl_rg = Label(window, text="Rg, ohms", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_rg.grid(row=13, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_rg = StringVar()
        self.str_rg.set(DEFAULT_RG)
        self.etr_rg = Entry(window, textvariable=self.str_rg, font=('Courier New', 12), width=15)
        self.etr_rg.grid(row=13, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        # Co
        self.lbl_co = Label(window, text="Co, nF", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_co.grid(row=14, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_co = StringVar()
        self.str_co.set(DEFAULT_CO)
        self.etr_co = Entry(window, textvariable=self.str_co, font=('Courier New', 12), width=15)
        self.etr_co.grid(row=14, column=1, rowspan=1, sticky=W, padx=2, pady=5)

        # Cf
        self.lbl_cf = Label(window, text="Cf, pF", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_cf.grid(row=15, column=0, rowspan=1, sticky=W, padx=50, pady=5)
        self.str_cf = StringVar()
        #self.str_cf.set(DEFAULT_CF)
        self.etr_cf = Entry(window, textvariable=self.str_cf, font=('Courier New', 12), width=15)
        self.etr_cf.grid(row=15, column=1, rowspan=1, sticky=W, padx=2, pady=5)


        # bottom frame
        fm = Frame(window)

        # cathode loadline check
        self.lbl_cathodeloadline = Label(fm, text="Cathode loadline", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_cathodeloadline.grid(row = 0, column = 0, pady=0, sticky='w')
        self.chk_cathodeloadline_var = IntVar()
        self.chk_cathodeloadline = Checkbutton(fm, variable=self.chk_cathodeloadline_var, onvalue = 1, offvalue = 0, height=1, font=('Courier New', 12),
                                               command=self.chk_cathodeloadline_click, background=self.window['bg'], width=1, anchor="w")
        self.chk_cathodeloadline.select()
        self.chk_cathodeloadline.grid(row = 0, column = 1, columnspan=1, sticky='W')

        # input signal swing check
        self.lbl_input_signal_swing = Label(fm, text="input signal swing", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_input_signal_swing.grid(row = 1, column = 0, pady=0, sticky='w')
        self.chk_input_signal_swing_var = IntVar()
        self.chk_input_signal_swing = Checkbutton(fm, variable=self.chk_input_signal_swing_var, onvalue = 1, offvalue = 0, height=1, font=('Courier New', 12),
                                               command=self.chk_input_signal_swing_click, background=self.window['bg'], width=1, anchor="w")
        self.chk_input_signal_swing.select()
        self.chk_input_signal_swing.grid(row = 1, column = 1, columnspan=1, sticky='W')

        # interpolation check
        self.lbl_showinterpolation = Label(fm, text="Show interpolation", font=('Courier New', 12), background=self.window['bg'])
        self.lbl_showinterpolation.grid(row = 2, column = 0, pady=0, sticky='w')
        self.chk_showinterpolation_var = IntVar()
        self.chk_showinterpolation = Checkbutton(fm, variable=self.chk_showinterpolation_var, onvalue = 1, offvalue = 0, height=1, font=('Courier New', 12),
                                               command=self.chk_showinterpolation_click, background=self.window['bg'], width=1, anchor="w")
        self.chk_showinterpolation.grid(row = 2, column = 1, columnspan=1, sticky='W')

        # valve specs label
        self.str_specs = StringVar()
        self.str_specs.set("Specs: ")
        self.lbl_specs = Label(fm, textvariable=self.str_specs, font=('Courier New', 12), background=self.window['bg'])
        self.lbl_specs.grid(row = 3, column = 0, columnspan=62, sticky='W')

        # calculations label
        self.str_calculations = StringVar()
        self.str_calculations.set("Calculations: ")
        self.lbl_calculations = Label(fm, textvariable=self.str_calculations, font=('Courier New', 12),  fg='green', background=self.window['bg'])
        self.lbl_calculations.grid(row = 4, column = 0, columnspan=62, sticky='W')

        # cursor position
        self.txt_coordinates = Text(fm, bd=0, bg=window['bg'], fg='red', height=1, wrap="none", state="normal", font=('Courier New', 12), background=self.window['bg'])
        self.txt_coordinates.grid(row = 5, column = 0, columnspan=10, sticky='W')
        self.txt_coordinates.config(highlightthickness = 0, borderwidth=0)
        self.txt_coordinates.config(state=DISABLED)

        fm.grid(row=45, column=0, padx=2, pady=0, columnspan=10, sticky='W')

        #BUTTONS
        fmbut = Frame(window)
        self.button_quit = Button(fmbut, text="QUIT", command=self.quit, font=('Courier New', 12))
        self.button_quit.grid(row=0, column=0, padx=10)
        self.button_start = Button(fmbut, text="PLOT", command=self.change_state, font=('Courier New', 12))
        self.button_start.grid(row=0, column=1, padx=10)
        self.button_clear = Button(fmbut, text="CLEAR", command=self.clear_chart, font=('Courier New', 12), state='normal')
        self.button_clear.grid(row=0, column=2, padx=10)
        self.but_export = Button(fmbut, text="EXPORT", command=self.export, font=('Courier New', 12))
        self.but_export.grid(row=0, column=3, padx=10)
        fmbut.grid(row=20, column=0, padx=50, pady=0, columnspan=5, sticky='W')

        # bind focus out events
        self.etr_Xmax.bind("<FocusOut>", self.parameters_changed)
        self.etr_Ra.bind("<FocusOut>", self.parameters_changed)
        self.etr_supply.bind("<FocusOut>", self.parameters_changed)
        self.etr_Rk.bind("<FocusOut>", self.parameters_changed)
        self.etr_Ymax.bind("<FocusOut>", self.parameters_changed)
        self.etr_rl.bind("<FocusOut>", self.parameters_changed)
        self.etr_ck.bind("<FocusOut>", self.parameters_changed)
        self.etr_inputsignal.bind("<FocusOut>", self.parameters_changed)
        self.etr_Xmax.bind("<Return>", self.parameters_changed)
        self.etr_Ra.bind("<Return>", self.parameters_changed)
        self.etr_supply.bind("<Return>", self.parameters_changed)
        self.etr_Rk.bind("<Return>", self.parameters_changed)
        self.etr_Ymax.bind("<Return>", self.parameters_changed)
        self.etr_rl.bind("<Return>", self.parameters_changed)
        self.etr_ck.bind("<Return>", self.parameters_changed)
        self.etr_inputsignal.bind("<Return>", self.parameters_changed)
        self.etr_co.bind("<FocusOut>", self.parameters_changed)
        self.etr_cf.bind("<FocusOut>", self.parameters_changed)
        self.etr_rg.bind("<FocusOut>", self.parameters_changed)
        self.etr_co.bind("<Return>", self.parameters_changed)
        self.etr_cf.bind("<Return>", self.parameters_changed)
        self.etr_rg.bind("<Return>", self.parameters_changed)
        #end of ui

        self.sechd = NONE
        self.swingl = NONE
        self.swingr = NONE

        # Read XMAX and YMAX from specs
        #self.updateMaxXY()

        self.valve_changed(None)
        #self.cmb_valve.focus_set()
        self.etr_supply.icursor(len(self.str_supply.get()))
        self.etr_supply.focus_set()

    def parameters_changed(self, event):
        self.change_state()

    def valve_changed(self, event):
        self.updateMaxXY()
        self.change_state()
        valve = self.str_valve.get()
        d=self.specs.loc[self.specs['valve']  == valve ]
        self.str_specs.set("Specs: mu: %s, ra: %s ohms, VaMax: %s V, Cga: %s pF, Cg to all but anode: %s pF, Pmax: %s W" % (
            d['mu'].iloc[0],
            d['ra'].iloc[0],
            d['VaMax'].iloc[0],
            d['Cga'].iloc[0],
            d['CgAEA'].iloc[0],
            d['Pmax'].iloc[0]
            ))

    # upodate xmax and ymax from valve specs
    def updateMaxXY(self):
        valve = self.str_valve.get()
        d=self.specs.loc[self.specs['valve']  == valve ]
        DEFAULT_XMAX = d['DefaultXmax'].iloc[0]
        DEFAULT_YMAX = d['DefaultYmax'].iloc[0]
        self.ax.set_ylim(0, d['DefaultYmax'].iloc[0])
        self.ax.set_xlim(0, d['DefaultXmax'].iloc[0])
        if len(self.etr_Ymax.get()) != 0 and self.can_convert_to_float(self.etr_Ymax.get()) == False: return
        if len(self.etr_Xmax.get()) != 0 and self.can_convert_to_float(self.etr_Xmax.get()) == False: return
        self.ax.set_yticks(range(0, int(self.etr_Ymax.get())+1, 5))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        #  self.ax.xaxis.set_ticks(np.arange(0, float(self.etr_Xmax.get()), float(self.etr_Xmax.get())/5), fontsize=20, visible=True)
        self.str_Ymax.set(DEFAULT_YMAX)
        self.str_Xmax.set(DEFAULT_XMAX)

    def clear_chart(self):
        self.ax.clear() # clear previous plot !!!!
        self.ax.set_yticks(range(0, int(self.etr_Ymax.get())+1, 5))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center', labelpad=20)
        self.ax.set_xlabel('Va, V', fontsize=16, loc='center',  labelpad=20)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def export(self):
        pass

    def chk_showinterpolation_click(self):
        self.change_state()

    def chk_cathodeloadline_click(self):
        self.change_state()

    def chk_input_signal_swing_click(self):
        if self.chk_input_signal_swing_var.get() == 0:
           self.etr_inputsignal.config({"background": "white"})
        else:
           self.etr_inputsignal.config({"background": "misty rose"})
        self.change_state()

    def quit(self):
        Tk().quit()

    def read_valvespecs(self):
        df = pd.read_csv('valves_specs.csv', comment='#', sep=',',  engine='python', skipinitialspace=True)
        #if DEBUG: print(df)
        self.specs = df

    def read_valvedata(self):
        pd.set_option('display.float_format', '{:20,.20f}'.format)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        #pd.set_option('display.width', None)
        #pd.set_option('display.max_colwidth', None)

        df = pd.read_csv('valves_data.csv', comment='#', sep='\r\n', names=['entire'], engine='python')
        df[['valve', 'curve', 'values']] = df['entire'].str.split(',', n=2, expand=True)
        df = df.drop(['entire'], axis=1)

        df['values'] = df['values'].apply(np.fromstring, sep=',')

        df = (df
          .assign(x=df['values'].str[::2], y=df['values'].str[1::2])
          .drop(columns='values')
          .explode(['x', 'y'])
        )
        df["x"] = pd.to_numeric(df["x"])
        df["y"] = pd.to_numeric(df["y"])

        df = df.set_index(['valve', 'curve' ])
        res = df.groupby(['valve', 'curve']).apply(lambda x: pd.Series(np.polyfit(x.x, x.y, 2), index=['a', 'b', 'c']))
        df = df.reset_index()
        self.df = pd.merge(df, res, how="inner", on=['valve', 'curve'])
        #if DEBUG: print(df)

    def motion_hover(self, event):
        if event.inaxes is not None:
            self.txt_coordinates.config(state=NORMAL)
            x = format(event.xdata, '.2f')
            y = format(event.ydata, '.2f')
            self.txt_coordinates.delete('1.0', END)
            self.txt_coordinates.insert(END, "(%s V, %s mA)" % (x, y))
            self.txt_coordinates.config(state=DISABLED)

    # function to plot loadline or refresh plot
    def change_state(self):
        # check values are valid
        if len(self.etr_Ra.get()) == 0: return
        if len(self.etr_supply.get()) == 0: return
        if len(self.etr_Rk.get()) == 0: return
        if len(self.etr_Xmax.get()) == 0: return
        if len(self.etr_Ymax.get()) == 0: return

        if len(self.etr_Xmax.get()) != 0 and self.can_convert_to_float(self.etr_Xmax.get()) == False: return
        if len(self.etr_Ra.get()) != 0 and self.can_convert_to_float(self.etr_Ra.get()) == False: return
        if len(self.etr_supply.get()) != 0 and self.can_convert_to_float(self.etr_supply.get()) == False: return
        if len(self.etr_Rk.get()) != 0 and self.can_convert_to_float(self.etr_Rk.get()) == False: return
        if len(self.etr_Ymax.get()) != 0 and self.can_convert_to_float(self.etr_Ymax.get()) == False: return
        if len(self.etr_rl.get()) != 0 and self.can_convert_to_float(self.etr_rl.get()) == False: return
        if len(self.etr_ck.get()) != 0 and self.can_convert_to_float(self.etr_ck.get()) == False: return
        if len(self.etr_inputsignal.get()) != 0 and self.can_convert_to_float(self.etr_inputsignal.get()) == False: return

        self.clear_chart()

        ymax = int(float(self.str_Ymax.get()))
        xmax = int(float(self.str_Xmax.get()))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        # select valve
        df = self.df
        valve = self.str_valve.get()
        de=df.loc[df['valve']  == valve ]

        # setup plot
        self.ax.set_title(valve, fontsize=24, pad=30, weight='bold')
        self.ax.set_yticks(range(0, int(self.etr_Ymax.get())+1, 5))
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        self.ax.grid(which="both", axis='both', color='slategray', linestyle='--', linewidth=0.7)
        #self.ax.set_ylabel('Ia, mA', fontsize=16, loc='center')
        #self.ax.set_xlabel('Va, V', fontsize=16, loc='center')

        """
        #por regresion lineal
        num_pts = 100
        da=de[['valve','curve','a','b','c']]
        res = da.drop_duplicates()
        for i in res.index:
            lower_limit = 0    #depende de curva
            x = np.linspace(lower_limit, xmax, num_pts)
            poly_coefs = res['a'][i], res['b'][i], res['c'][i]
            y = np.polyval(np.asarray(poly_coefs), x)
            #plt.plot(x, y, '-0')
            ind = np.argmin(y)
            x_min = x[ind]
            ax.set_ylim(0, ymax)
            ax.set_xlim(0, xmax)
            y_positive = y
            x_positive = x
            if res['curve'][i].strip() != 'Pmax':
                y_positive = y[x >= x_min]
                x_positive = x[x >= x_min]
            ax.plot(x_positive, y_positive, linewidth=2, label=res['curve'][i].strip())
        """

        # plot cathode line
        df = self.df
        valve = self.str_valve.get()
        de = df.loc[df['valve']  == valve ]
        de = de.drop(de[de['curve'] == ' Pmax'].index)
        data = de[['curve','a','b','c'] ].drop_duplicates()
        data[['curve', 'a', 'b', 'c']] = data[['curve', 'a', 'b', 'c']].apply(pd.to_numeric)
        data['y'] = - data['curve'] / float(self.str_Rk.get()) * 1000
        data['x'] = (-data['b']+(data['b']**2-4*data['a']*(data['c']-data['y']))**0.5)/(2*data['a'])
        x_positive = data['x']
        y_positive = data['y']

        # extrapolation
        interval = int(float(self.etr_supply.get()) / 0.01)
        xnew = np.linspace(x_positive.min(), x_positive.max(), interval)
        spl = make_interp_spline(x_positive, y_positive, k=1)  # type: BSpline
        ynew = spl(xnew)

        if self.chk_cathodeloadline_var.get() == 1:
            self.ax.plot(xnew, ynew, '-', color='green', linewidth=1)
            #self.ax.plot(data['x'], data['y'], '-', color='green', linewidth=1)

        # quiscient
        b = float(self.str_supply.get()) * 1000 / float(self.str_Ra.get())
        a = -b / float(self.str_supply.get())
        yloadline = xnew * a + b
        idx = np.argwhere(np.diff(np.sign(ynew - yloadline))).flatten()
        self.vq = xnew[idx].astype("float")[0]
        self.str_Vq.set(format(self.vq, ".2f"))
        self.iq = yloadline[idx].astype("float")[0]
        self.str_Iq.set(format(self.iq, ".2f"))

        # Vgk
        self.vgk = - self.iq / 1000 * float(self.str_Rk.get())
        self.str_Vgk.set(format(self.vgk, ".2f"))

        # get input signal swing values
        if self.chk_input_signal_swing_var.get() == 1:
            self.input_swing()

        # plot grid lines using existent points and interpolation
        self.ax.set_ylim(0, ymax)
        self.ax.set_xlim(0, xmax)
        valve = self.str_valve.get()
        da = df.loc[df['valve']  == valve ]
        da=da.loc[(da['curve'] != ' Pmax') & (da['x'] <= xmax) & (da['y'] <= ymax), ['valve','curve','x','y'] ]
        for name, g in da.groupby('curve'):
            x_positive = g['x']
            y_positive = g['y']

            # try interpolation
            try:
                interval = int(float(self.etr_supply.get()) / 0.01)
                xnew = np.linspace(x_positive.min(), x_positive.max(), interval)
                spl = make_interp_spline(x_positive, y_positive, k=2)  # type: BSpline
                ynew = spl(xnew)
            except:
                # data points not enought ?
                xnew = x_positive
                ynew = y_positive
            self.ax.plot(xnew, ynew, '-', color='black', linewidth=1)

            # annotate grid curves
            if name == ' 0': name = "Vg = 0"
            self.ax.annotate(name + 'V', xy=(x_positive.max(), y_positive.max()), rotation=0, fontweight='bold')

        # plot Pmax using points
        da=df.loc[(df['curve'] == ' Pmax') & (df['valve']  == valve), ['valve','curve','x','y'] ]
        x_positive = da['x']
        y_positive = da['y']
        # try interpolation
        xnew = np.linspace(x_positive.min(), x_positive.max(), 300)
        spl = make_interp_spline(x_positive, y_positive, k=2)  # type: BSpline
        ynew = spl(xnew)
        # add Pmax legend
        valve = self.str_valve.get()
        pmax=self.specs.loc[self.specs['valve']  == valve ]['Pmax'].iloc[0]
        self.ax.plot(xnew, ynew, 'r--', linewidth=1, label=da['curve'].iloc[0].strip() + ' = %sW' % str(pmax) )

        # plot loadline
        x_values = [float(self.str_supply.get()), 0]
        y_values = [0, float(self.str_supply.get()) * 1000 / float(self.str_Ra.get())]
        self.ax.plot(x_values, y_values, '-', color='blue', linewidth=2)

        # plot input signal suing
        if self.swingl != NONE and self.chk_input_signal_swing_var.get() == 1:
            self.input_swing()
            Vs = float(self.str_supply.get())
            Ra = float(self.str_Ra.get())
            b = Vs * 1000 / Ra
            a = -b / Vs
            self.ax.plot([ (self.swingl[1] - b)/a, (self.swingr[1] - b)/a ], [ self.swingl[1], self.swingr[1] ], '-', color='orange', linewidth=2.7)

        # plot quiscient at the end
        self.ax.plot(self.vq, self.iq, 'ro', markersize=6)

        self.ax.legend(loc='upper left')
        self.canvas.draw()

        # calculate_gain_impedance
        self.calculate_gain_impedance()

        # calculate and show 2HD
        self.calculate_2hd()
        if self.sechd != NONE:
            self.str_calculations.set(self.str_calculations.get() + ', 2HD: %.2f %%' % self.sechd)

    def input_swing(self):
        self.swingl = NONE
        self.swingr = NONE
        inputVpp = self.str_inputsignal.get()
        if len(inputVpp) == 0: return
        valve = self.str_valve.get()
        inputVpp = float(inputVpp)
        Vs = float(self.str_supply.get())
        Ra = float(self.str_Ra.get())
        lft_point_y = self.vgk+inputVpp/2000
        rht_point_y = self.vgk-inputVpp/2000

        if lft_point_y > 0: return

        da = self.df.loc[(self.df['valve'] == valve) & (self.df['curve'] != ' Pmax'), ['valve','curve','x','y'] ]
        da[['curve', 'x', 'y']] = da[['curve', 'x', 'y']].apply(pd.to_numeric, errors='coerce', axis=1)
        b = Vs * 1000 / Ra
        a = -b / Vs

        # left point of swing
        de = da.loc[(da['curve'] < lft_point_y), ['curve'] ]
        lft_prev_curve = de['curve'].max()
        de = da.loc[(da['curve'] > lft_point_y), ['curve'] ]
        lft_next_curve = de['curve'].min()

        de = da.loc[(da['valve'] == valve) & (da['curve'] == lft_prev_curve), ['x'] ]
        xmin = de['x'].min()

        de = da.loc[(da['curve'] < rht_point_y), ['curve'] ]
        rht_prev_curve = de['curve'].max()
        de = da.loc[(da['curve'] > rht_point_y), ['curve'] ]
        rht_next_curve = de['curve'].min()
        de = da.loc[(da['valve'] == valve) & (da['curve'] == rht_next_curve), ['x'] ]
        xmax = de['x'].max()

        if lft_point_y >= float(self.str_Xmax.get()): return
        #if DEBUG: print('swing - min: %, max: %' , xmin, xmax)

        interp_A = lft_point_y
        interp_B = np.arange(xmin, xmax, 0.01) #step=0.1
        ynew = gd((da['curve'], da['x']), da['y'], (interp_A, interp_B), method='cubic')
        yloadline = interp_B * a + b

        idx = np.argwhere(np.diff(np.sign(yloadline - ynew))).flatten()

        # left point of swing
        if self.chk_showinterpolation_var.get() == 1:
        #if DEBUG:
        #    print('left swing: Vg=%s Ia=%s Va=%s' % (lft_point_y, ynew, interp_B))
            self.ax.plot(interp_B, ynew, '-', color='green', linewidth=1)
        self.swingl = [lft_point_y, yloadline[idx][0], (yloadline[idx][0]-b)/a   ]
        #print(lft_point_y)
        #print(yloadline[idx][0])
        #print((yloadline[idx][0]-b)/a)

        interp_A = rht_point_y
        interp_B = np.arange(xmin, xmax, 0.01) #step=0.1
        ynew = gd((da['curve'], da['x']), da['y'], (interp_A, interp_B), method='cubic')

        idx = np.argwhere(np.diff(np.sign(ynew - yloadline))).flatten()

        # right point of swing
        if self.chk_showinterpolation_var.get() == 1:
        #if DEBUG:
        #    print('right swing: Vg=%s Ia=%s Va=%s' % (rht_point_y, ynew, interp_B))
             self.ax.plot(interp_B, ynew, '-', color='green', linewidth=1)
        self.swingr = [rht_point_y, yloadline[idx][-1], (yloadline[idx][-1]-b)/a  ]
        #print(rht_point_y)
        #print(yloadline[idx][-1])
        #print((yloadline[idx][-1]-b)/a)


    def can_convert_to_float(self, string):
        try:
            result = float(string)
            return True
        except ValueError:
            return False

    def calculate_2hd(self):
        self.sechd = NONE
        if len(self.str_Vq.get()) == 0: return
        if self.swingl != NONE or self.swingr != NONE:
            A = self.swingl[2]
            B = float(self.str_Vq.get())
            C = self.swingr[2]
            AB = B - A
            BC = C - B
            if DEBUG:
                print('\n\n\nAB: %f' % AB)
                print('BC: %f' % BC)
                print('A: %f' % A)
                print('B: %f' % B)
                print('C: %f' % C)
            self.sechd = abs((AB - BC) / (2 * (AB+BC)) * 100)
            #print('2: %f' % self.sechd)
            #self.sechd = (( B - ((C+A)/2) ) / (C - A) ) * 100
            #print('2: %f' % self.sechd)

    # calculate gain, power diss, output impedance, etc.
    def calculate_gain_impedance(self):
        valve = self.str_valve.get()
        mu = self.specs.loc[self.specs['valve']  == valve ]['mu'].iloc[0]
        ra = self.specs.loc[self.specs['valve']  == valve ]['ra'].iloc[0]
        if len(self.etr_Ra.get()) == 0: return
        Ra = float(self.etr_Ra.get())
        if len(self.etr_rl.get()) == 0: return
        if len(self.etr_rl.get()) != 0 and self.can_convert_to_float(self.etr_rl.get()) == False: return
        Rl = float(self.etr_rl.get())

        #sin cap
        if len(self.etr_ck.get()) == 0 or self.can_convert_to_float(self.etr_ck.get()) == False or self.etr_ck.get() == '0':
            if len(self.etr_Rk.get()) == 0: return
            if len(self.etr_Rk.get()) != 0 and self.can_convert_to_float(self.etr_Rk.get()) == False: return
            Rk = float(self.etr_Rk.get())
            g=(mu*(1/((1/Ra)+(1/Rl))))/((1/((1/Ra)+(1/Rl)))+ra+(Rk*(mu+1)))
            gdb=(math.log10(g))*20
            aoi =(Ra*(ra+(Rk*(mu+1))))/(Ra+ra+(Rk*(mu+1)))
            coi = 1/((1/((Ra+ra)/(mu+1)))+(1/Rk))
            self.str_calculations.set('A: %.1f %.1f dB, P: %.2f W, anode output impedance: %d ohms, cathode output impedance: %d ohms' % (g, gdb, self.vq * self.iq / 1000, aoi, coi))
            Cga = self.specs.loc[self.specs['valve']  == valve ]['Cga'].iloc[0]
            CgAEA = self.specs.loc[self.specs['valve']  == valve ]['CgAEA'].iloc[0]
            Cf = 0
            if len(self.str_cf.get()) != 0 and self.can_convert_to_float(self.str_cf.get()) == TRUE:
                Cf = float(self.str_cf.get())
            total_input_capacitance = CgAEA +((Cf+Cga)*g)
            self.str_calculations.set(self.str_calculations.get() + ', Total_Input_Capacitance: %d pf' % total_input_capacitance)
        elif len(self.etr_ck.get()) != 0 and self.can_convert_to_float(self.etr_ck.get()) != False:
            gcc=(mu*(1/((1/Ra)+(1/Rl))))/((1/((1/Ra)+(1/Rl)))+ra)
            gdb=(math.log10(gcc))*20
            aoi = (Ra*ra)/(Ra+ra)
            self.str_calculations.set('A: %.1f %.1f dB, P: %.2f W, anode output impedance: %d ohms' % (gcc, gdb, self.vq * self.iq / 1000, aoi))
            Cga = self.specs.loc[self.specs['valve']  == valve ]['Cga'].iloc[0]
            CgAEA = self.specs.loc[self.specs['valve']  == valve ]['CgAEA'].iloc[0]
            Cf = 0
            if len(self.str_cf.get()) != 0 and self.can_convert_to_float(self.str_cf.get()) == TRUE:
                Cf = float(self.str_cf.get())
            total_input_capacitance = CgAEA +((Cf+Cga)*gcc)
            self.str_calculations.set(self.str_calculations.get() + ', Total_Input_Capacitance: %d pf' % total_input_capacitance)


        """
        #sin cap
        g=(mu*(1/((1/Ra)+(1/Rl))))/((1/((1/Ra)+(1/Rl)))+ra+(Rk*(mu+1)))
        gdb=(LOG10(g))*20
        anode output impedance =(Ra*(ra+(Rk*(mu+1))))/(Ra+ra+(Rk*(mu+1)))
        cathode output impedance = 1/((1/((Ra+ra)/(mu+1)))+(1/Rk))
        total input capacitance = Cgaea +((Cf+Cga)*g)
        #####################################
        #con cap
        gcc=(mu*(1/((1/Ra)+(1/Rl))))/((1/((1/Ra)+(1/Rl)))+ra)
        gdb=(LOG10(gcc))*20
        anode output impedance = (Ra*ra)/(Ra+ra)
        total input capacitance = Cgaea +((Cf+Cga)*gcc)
        #####################################
        HF roll-off due to grid stopping (-3dB) (Hz) =1/(2*PI()*(Total Input Capacitance/1000000000000)*Rg)
        LF output roll-off due to Co and Rl  (-3dB) (Hz) =1/(2*PI()*(Anode Output Impedance con CAP+Rl)*(Co/1000000000))
        Half boost freq. Hz=(1/(2*PI()*Rk*(Ck/1000000)))*(1+((Rk*mu+1)/((2*(Ra+ra))+(0.5*(Rk*mu+1)))))^0.5
        """

window = Tk()
window.title('andmarti Load Line Plotter')
#window.state('zoomed')
start = mclass(window)
window.mainloop()
